"""Frozen final-layer ESM-2 residue-mean embeddings, one row per accession.

Project-local environment: .venv/Scripts/python.exe code/generate_esm2_embeddings.py
Run --pilot first; then run without --pilot. Interrupted runs require --resume.
No labels, classifiers, split generation, or similarity recomputation are used.
"""
import argparse
from contextlib import nullcontext
from datetime import datetime, timezone
import hashlib
import json
import logging
import os
from pathlib import Path
import platform
import sqlite3
import sys
import time

import numpy as np
import pandas as pd
import psutil
import torch
from filelock import FileLock
from tqdm import tqdm
import transformers
from transformers import AutoConfig, AutoModel, AutoTokenizer

ROOT = Path(__file__).resolve().parent.parent
MODEL = "facebook/esm2_t33_650M_UR50D"
REVISION = "08e4846e537177426273712802403f7ba8261b6c"
DIMENSION = 1280
MAX_LENGTH = 1022


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    temporary.replace(path)


def utc():
    return datetime.now(timezone.utc).isoformat()


def hardware():
    result = dict(python=platform.python_version(), platform=platform.platform(),
                  torch=torch.__version__, transformers=transformers.__version__,
                  numpy=np.__version__, cuda_available=torch.cuda.is_available(),
                  pytorch_cuda_version=torch.version.cuda,
                  cudnn_version=torch.backends.cudnn.version(),
                  cpu_logical_cores=os.cpu_count(), ram_bytes=psutil.virtual_memory().total,
                  available_ram_bytes=psutil.virtual_memory().available, gpus=[])
    for i in range(torch.cuda.device_count()):
        free, total = torch.cuda.mem_get_info(i)
        result["gpus"].append(dict(index=i, name=torch.cuda.get_device_name(i),
                                   free_vram_bytes=free, total_vram_bytes=total,
                                   compute_capability=list(torch.cuda.get_device_capability(i))))
    return result


def load_proteins(source):
    # Do not load EC3, taxonomy, split membership, or any other label column.
    frame = pd.read_csv(source, sep="\t", usecols=["accession", "sequence"],
                        dtype=str, keep_default_na=False)
    if frame.accession.eq("").any() or frame.accession.duplicated().any():
        raise ValueError("Missing or duplicated accession in cleaned dataset")
    if frame.sequence.eq("").any() or not frame.sequence.str.fullmatch("[A-Z]+").all():
        raise ValueError("Missing or invalid protein sequence; input is not altered")
    frame = frame.sort_values("accession", kind="stable").reset_index(drop=True)
    frame["length"] = frame.sequence.str.len()
    bad = frame[frame.length.gt(MAX_LENGTH)]
    if len(bad):
        raise ValueError(f"Overlength sequences; refusing truncation: {bad[['accession','length']].to_dict('records')}")
    frame["sequence_sha256"] = frame.sequence.map(lambda s: hashlib.sha256(s.encode("ascii")).hexdigest())
    frame.index.name = "row_index"
    return frame


def frozen_inputs(source):
    paths = {source.resolve()}
    paths.update((ROOT/"data").glob("*"))
    paths.update((ROOT/"splits_ec3").rglob("*.tsv"))
    paths.update((ROOT/"splits_ec3").rglob("settings.json"))
    for name in ["nearest_train_identity", "similarity_bins"]:
        directory = ROOT/"analysis"/name
        paths.update(p for p in directory.glob("*") if p.is_file())
        if name == "nearest_train_identity":
            for p in directory.glob("*/*.tsv"):
                if not p.name.endswith("_hits.tsv"):
                    paths.add(p)
    return {str(p.resolve()): digest(p) for p in sorted(paths) if p.is_file()}


def residue_mean(hidden, attention_mask, special_mask, lengths):
    mask = attention_mask.bool() & ~special_mask.bool()
    expected = torch.as_tensor(lengths, device=hidden.device)
    if not torch.equal(mask.sum(dim=1), expected):
        raise ValueError("Residue-token count differs from amino-acid length")
    pooled = (hidden.float() * mask.unsqueeze(-1)).sum(dim=1) / expected.unsqueeze(-1)
    if not torch.isfinite(pooled).all():
        raise ValueError("Nonfinite pooled embedding")
    return pooled


def embed(frame, rows, tokenizer, model, device, precision):
    sequences = frame.loc[rows, "sequence"].tolist()
    encoded = tokenizer(sequences, padding=True, truncation=False,
                        return_special_tokens_mask=True, return_tensors="pt")
    special = encoded.pop("special_tokens_mask")
    mask = encoded["attention_mask"].bool() & ~special.bool()
    for k, sequence in enumerate(sequences):
        ids = encoded["input_ids"][k][mask[k]].tolist()
        expected = tokenizer.convert_tokens_to_ids(list(sequence))
        if ids != expected or len(ids) != len(sequence) or tokenizer.unk_token_id in ids:
            raise ValueError(f"Tokenization or unknown-residue error: {frame.loc[rows[k],'accession']}")
        if int(encoded["attention_mask"][k].sum()) != len(sequence) + 2:
            raise ValueError("Unexpected special-token layout")
    if encoded["input_ids"].shape[1] > MAX_LENGTH + 2:
        raise ValueError("Tokenized input exceeds the declared limit")
    encoded = {k:v.to(device) for k,v in encoded.items()}
    context = torch.autocast("cuda", dtype=torch.float16) if precision == "float16" else nullcontext()
    with torch.inference_mode(), context:
        hidden = model(**encoded, output_hidden_states=False, output_attentions=False,
                       return_dict=True).last_hidden_state
        if hidden.shape[-1] != DIMENSION:
            raise ValueError("Unexpected final hidden dimension")
        pooled = residue_mean(hidden, encoded["attention_mask"], special.to(device),
                              [len(s) for s in sequences])
    result = pooled.cpu().numpy().astype(np.float32, copy=True)
    del hidden, pooled, encoded
    return result


def check_ledger(connection, matrix, frame):
    completed = set()
    for row, accession, sequence_hash, vector_hash in connection.execute(
            "SELECT row_index,accession,sequence_sha256,vector_sha256 FROM completed"):
        if row < 0 or row >= len(frame):
            raise ValueError("Out-of-range checkpoint row")
        if accession != frame.loc[row,"accession"] or sequence_hash != frame.loc[row,"sequence_sha256"]:
            raise ValueError("Checkpoint accession or sequence mismatch")
        vector = np.asarray(matrix[row], dtype=np.float32)
        if not np.isfinite(vector).all() or hashlib.sha256(vector.tobytes()).hexdigest() != vector_hash:
            raise ValueError(f"Completed embedding failed checksum: {accession}")
        completed.add(row)
    return completed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT/"data"/"proteins_ec3.tsv")
    parser.add_argument("--output", type=Path, default=ROOT/"analysis"/"esm2_embeddings")
    parser.add_argument("--cache-dir", type=Path, default=ROOT/".cache"/"huggingface")
    parser.add_argument("--revision", default=REVISION, help="Immutable public HF checkpoint revision")
    parser.add_argument("--device", choices=["auto","cuda","cpu"], default="auto")
    parser.add_argument("--precision", choices=["auto","float32","float16"], default="auto")
    parser.add_argument("--batch-size", type=int, default=None)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--pilot", action="store_true", help="Eight length-spanning proteins in a separate pilot output")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--restart", action="store_true", help="Start in a new timestamped sibling, preserving the old run")
    parser.add_argument("--recompute-samples", type=int, default=5)
    args = parser.parse_args()
    if args.resume and args.restart:
        parser.error("--resume and --restart are mutually exclusive")
    if args.threads < 1 or args.recompute_samples < 1 or args.batch_size is not None and args.batch_size < 1:
        parser.error("Thread, batch, and recomputation counts must be positive")
    torch.set_num_threads(args.threads)
    torch.manual_seed(355)
    np.random.seed(355)
    hw = hardware()
    print(json.dumps(hw,indent=2),flush=True)
    device = "cuda" if args.device == "auto" and hw["cuda_available"] else args.device
    device = "cpu" if device == "auto" else device
    if device == "cuda" and not hw["cuda_available"]:
        raise RuntimeError("CUDA requested but unavailable")
    precision = ("float16" if device == "cuda" else "float32") if args.precision == "auto" else args.precision
    if precision == "float16" and device != "cuda":
        raise ValueError("CPU fallback uses float32")
    batch_size = args.batch_size or (4 if device == "cuda" and hw["gpus"][0]["free_vram_bytes"] >= 12*2**30 else 1)
    out = args.output.resolve()
    if args.pilot:
        out = out.with_name(out.name + "_pilot")
    if args.restart:
        out = out.with_name(out.name + "_restart_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S"))
    if out.parent != (ROOT/"analysis").resolve() or not out.name.startswith("esm2_embeddings"):
        raise ValueError("Output must be a separate analysis/esm2_embeddings* directory")
    source = args.source.resolve()
    frame = load_proteins(source)
    all_lengths = frame.length.to_numpy()
    all_count = len(frame)
    if args.pilot:
        order = frame.sort_values(["length","accession"]).index.to_numpy()
        selected = order[np.linspace(0,len(order)-1,8,dtype=int)]
        frame = frame.loc[selected].sort_values("accession").reset_index(drop=True)
        frame.index.name = "row_index"
    if out.exists() and not args.resume:
        raise FileExistsError(f"Output exists: {out}; use --resume or --restart")
    if args.resume and not (out/"settings.json").exists():
        raise FileNotFoundError("Resume requires an existing configured run")
    out.mkdir(parents=True,exist_ok=True)
    logging.basicConfig(level=logging.INFO,format="%(asctime)s %(message)s",
                        handlers=[logging.FileHandler(out/"run.log",encoding="utf-8"),logging.StreamHandler(sys.stdout)])
    with FileLock(str(out/".run.lock"),timeout=0):
        cached = json.loads((out/"settings.json").read_text()) if args.resume else None
        protected = cached["frozen_input_sha256"] if cached else frozen_inputs(source)
        if any(digest(p) != h for p,h in protected.items()):
            raise ValueError("A frozen input differs from its saved hash")
        config = AutoConfig.from_pretrained(MODEL,revision=args.revision,cache_dir=args.cache_dir,token=False)
        if config.hidden_size != DIMENSION or config.num_hidden_layers != 33:
            raise ValueError("Checkpoint dimension or depth does not match requested ESM-2")
        resolved = getattr(config,"_commit_hash",None)
        if resolved != args.revision or len(resolved) != 40:
            raise ValueError("An immutable resolved checkpoint SHA is required")
        tokenizer = AutoTokenizer.from_pretrained(MODEL,revision=resolved,cache_dir=args.cache_dir,
                                                  use_fast=False,token=False)
        alphabet = set("".join(frame.sequence))
        if any(tokenizer.convert_tokens_to_ids(a) in tokenizer.all_special_ids for a in alphabet):
            raise ValueError("Input contains residues unsupported by checkpoint vocabulary")
        identity = dict(model=MODEL,resolved_revision=resolved,source_sha256=digest(source),
                        input_count=len(frame),full_cleaned_input_count=all_count,
                        pilot=args.pilot,hidden_dimension=DIMENSION,final_dtype="float32",
                        inference_precision=precision,device=device,
                        torch_version=torch.__version__,transformers_version=transformers.__version__,
                        pooling="last_hidden_state mean over attention_mask AND NOT special_tokens_mask",
                        row_order="lexicographic accession",source_columns=["accession","sequence"],
                        truncation=False,maximum_allowed_residues=MAX_LENGTH)
        if cached and cached["identity"] != identity:
            raise ValueError("Resume identity/model/runtime settings differ; use a separate restart output")
        metadata = frame[["accession","length","sequence_sha256"]].reset_index()
        metadata_path = out/"accessions.tsv"
        if cached:
            saved = pd.read_csv(metadata_path,sep="\t",dtype={"accession":str,"sequence_sha256":str})
            pd.testing.assert_frame_equal(saved,metadata)
        else:
            metadata.to_csv(metadata_path,sep="\t",index=False)
            settings = dict(identity=identity,hardware=hw,initial_batch_size=batch_size,
                            threads=args.threads,model_config=config.to_dict(),
                            frozen_input_sha256=protected,script_sha256=digest(Path(__file__)),
                            length_summary=dict(count=all_count,min=int(all_lengths.min()),
                            max=int(all_lengths.max()),mean=float(all_lengths.mean()),
                            median=float(np.median(all_lengths)),
                            quantiles={str(q):float(np.quantile(all_lengths,q)) for q in [.05,.25,.75,.95,.99]}),
                            created_utc=utc(),no_labels_used=True,model_frozen=True)
            atomic_json(out/"settings.json",settings)
        model = AutoModel.from_pretrained(MODEL,revision=resolved,cache_dir=args.cache_dir,token=False,
                    use_safetensors=True,torch_dtype=torch.float32,low_cpu_mem_usage=True,
                    add_pooling_layer=False,attn_implementation="eager")
        model.requires_grad_(False)
        model.eval()
        model.to(device)
        if model.training or any(p.requires_grad for p in model.parameters()):
            raise ValueError("Model is not frozen")
        matrix_path = out/"embeddings.npy"
        if matrix_path.exists():
            matrix = np.load(matrix_path,mmap_mode="r+")
            if matrix.shape != (len(frame),DIMENSION) or matrix.dtype != np.float32:
                raise ValueError("Embedding matrix header mismatch")
        else:
            matrix = np.lib.format.open_memmap(matrix_path,mode="w+",dtype="float32",shape=(len(frame),DIMENSION))
            for start in range(0,len(frame),256):
                matrix[start:start+256] = np.nan
            matrix.flush()
        connection = sqlite3.connect(out/"checkpoint.sqlite3")
        connection.execute("PRAGMA journal_mode=WAL")
        connection.execute("PRAGMA synchronous=FULL")
        connection.execute("CREATE TABLE IF NOT EXISTS completed (row_index INTEGER PRIMARY KEY, accession TEXT UNIQUE NOT NULL, sequence_sha256 TEXT NOT NULL, vector_sha256 TEXT NOT NULL)")
        done = check_ledger(connection,matrix,frame)
        pending = [int(r) for r in frame.sort_values(["length","accession"]).index if r not in done]
        logging.info("Model loaded. %d/%d validated rows; batch size %d; %s %s",len(done),len(frame),batch_size,device,precision)
        atomic_json(out/"progress.json",dict(status="running",completed=len(done),expected=len(frame),updated_utc=utc()))
        started = time.perf_counter()
        processed_residues = 0
        progress = tqdm(total=len(frame),initial=len(done),unit="protein",mininterval=10)
        try:
            for start in range(0,len(pending),batch_size):
                rows = pending[start:start+batch_size]
                vectors = embed(frame,rows,tokenizer,model,device,precision)
                receipts = []
                for row, vector in zip(rows,vectors):
                    matrix[row] = vector
                    receipts.append((row,frame.loc[row,"accession"],frame.loc[row,"sequence_sha256"],
                                     hashlib.sha256(vector.tobytes()).hexdigest()))
                matrix.flush()
                with matrix_path.open("r+b") as handle:
                    os.fsync(handle.fileno())
                with connection:
                    connection.executemany("INSERT INTO completed VALUES (?,?,?,?)",receipts)
                done.update(rows)
                processed_residues += int(frame.loc[rows,"length"].sum())
                progress.update(len(rows))
                elapsed = time.perf_counter()-started
                atomic_json(out/"progress.json",dict(status="running",completed=len(done),expected=len(frame),
                    elapsed_this_run_seconds=elapsed,processed_residues_this_run=processed_residues,
                    batch_size=batch_size,device=device,inference_precision=precision,updated_utc=utc()))
                if len(done) % 100 == 0 or args.pilot:
                    logging.info("Checkpoint %d/%d; %.2f seconds elapsed",len(done),len(frame),elapsed)
        except BaseException as error:
            atomic_json(out/"progress.json",dict(status="interrupted_or_failed",completed=len(done),expected=len(frame),
                                                exception_type=type(error).__name__,message=str(error),updated_utc=utc()))
            logging.exception("Job interrupted; completed rows are resumable")
            raise
        finally:
            progress.close()
        check_ledger(connection,matrix,frame)
        if len(done) != len(frame):
            raise ValueError("Incomplete accession coverage")
        duplicate_groups = connection.execute("SELECT vector_sha256,COUNT(*),COUNT(DISTINCT sequence_sha256) FROM completed GROUP BY vector_sha256 HAVING COUNT(*)>1").fetchall()
        if any(distinct != 1 for _,_,distinct in duplicate_groups):
            raise ValueError("Identical embedding rows for different sequences; possible indexing error")
        count = 0
        total = square_total = 0.0
        minimum, maximum = float("inf"),float("-inf")
        for start in range(0,len(frame),256):
            block = np.asarray(matrix[start:start+256],dtype=np.float64)
            if not np.isfinite(block).all():
                raise ValueError("NaN/Inf in completed embedding matrix")
            count += block.size
            total += float(block.sum())
            square_total += float(np.square(block).sum())
            minimum,maximum = min(minimum,float(block.min())),max(maximum,float(block.max()))
        rng = np.random.default_rng(355)
        samples = sorted(set(rng.choice(len(frame),min(args.recompute_samples,len(frame)),replace=False).tolist()
                             + [int(frame.length.idxmax())]))
        tolerance = dict(rtol=5e-3,atol=5e-3) if precision == "float16" else dict(rtol=1e-4,atol=2e-5)
        checks = []
        for row in samples:
            vector = embed(frame,[row],tokenizer,model,device,precision)[0]
            difference = float(np.max(np.abs(vector-matrix[row])))
            passed = bool(np.allclose(vector,matrix[row],**tolerance))
            checks.append(dict(row_index=row,accession=frame.loc[row,"accession"],max_absolute_difference=difference,passed=passed))
        unchanged = all(digest(p) == h for p,h in protected.items())
        validation = dict(all_passed=all(c["passed"] for c in checks) and unchanged,
            input_proteins=len(frame),embeddings_generated=len(done),unique_accessions=int(frame.accession.nunique()),
            every_accession_exactly_once=True,shape=list(matrix.shape),dtype=str(matrix.dtype),dimension_verified=DIMENSION,
            all_values_finite=True,failed_proteins=0,sequences_truncated=0,
            embedding_value_statistics=dict(min=minimum,max=maximum,mean=total/count,
                std=float(np.sqrt(max(0,square_total/count-(total/count)**2)))),
            duplicate_vector_groups_for_identical_sequences=len(duplicate_groups),unexpected_duplicate_vector_groups=0,
            recomputation_checks=checks,recomputation_tolerance=tolerance,frozen_inputs_unchanged=unchanged,
            inference_seconds_this_run=time.perf_counter()-started,processed_residues_this_run=processed_residues,
            matrix_bytes=matrix_path.stat().st_size,matrix_sha256=digest(matrix_path),
            no_labels_used=True,model_frozen=True,classifier_trained=False,finished_utc=utc())
        atomic_json(out/"validation.json",validation)
        if not validation["all_passed"]:
            raise ValueError("Recomputation or frozen-input validation failed")
        atomic_json(out/"progress.json",dict(status="complete",completed=len(done),expected=len(frame),updated_utc=utc()))
        connection.execute("PRAGMA wal_checkpoint(TRUNCATE)")
        connection.close()
        report = f"""# Frozen ESM-2 embeddings{' pilot' if args.pilot else ''}

Completed and validated {len(frame):,} accession-indexed embeddings from the cleaned dataset.
Model: `{MODEL}`; resolved revision `{resolved}`; final hidden dimension {DIMENSION}.
Matrix shape: `{matrix.shape}`; float32 NumPy .npy memory map; {matrix_path.stat().st_size:,} bytes.
Device: {device}; inference precision: {precision}; batch size this run: {batch_size}; CPU threads: {args.threads}.
Pooling: arithmetic mean of final-layer residue states, excluding BOS/CLS, EOS and padding using attention and special-token masks. No L2 normalization or sequence truncation.

Accession rows are sorted lexicographically; `accessions.tsv` includes row index, amino-acid length and sequence hash. Labels and taxonomy are not loaded or used for embedding generation. ESM-2 is frozen in evaluation and inference mode. No classifier was trained or evaluated. Embeddings are generated once from the cleaned dataset and can be joined by accession to any held-out experiment.

All values are finite. Every accession has exactly one checksummed completed row. Identical rows are allowed only for identical source sequences; unexpected duplicate vectors: zero. Independent single-protein recomputation checks: {len(checks)}, all within the tolerance in validation.json. Frozen dataset, split, nearest-identity and stratum inputs are unchanged.

Files: embeddings.npy, accessions.tsv, checkpoint.sqlite3, settings.json, progress.json, validation.json, run.log and report.md. The SQLite ledger commits a row only after the matrix is flushed; resume verifies completed row hashes and source mapping before skipping them. Do not use unfinished rows; progress.json must say complete and validation.json must pass. Exact source/runtime configuration is enforced on resume; changed numerical settings require a separate restart output.

Reproduce using the project-local Python environment. Run `code/generate_esm2_embeddings.py --pilot` before a full run; use `--resume` after interruption, or `--restart` to start in a new timestamped sibling without destroying the previous run. Public checkpoint downloads need no API key and use the local cache. CPU float32 is supported; optional CUDA mixed float16 inference still pools and saves in float32. Mixed precision and different hardware can produce small numerical differences. Length sorting changes execution order only, not matrix/accession ordering. No model checkpoint smaller than the requested 650M model is substituted.
"""
        (out/"report.md").write_text(report,encoding="utf-8")
        print(json.dumps(validation,indent=2),flush=True)


if __name__ == "__main__":
    main()
