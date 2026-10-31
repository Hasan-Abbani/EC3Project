# Frozen ESM-2 embedding execution and GPU handoff

## Scope

The local 8 GB Windows computer has CPU-only PyTorch and Intel integrated
graphics, with no CUDA device. Only a CPU pilot is authorized on that machine.
Do not run the full 146,677-protein job there. The full job is intended for a
separate CUDA-capable machine. No classifier, embedding fine-tuning, or
evaluation of test-set predictions is part of this stage.

## Fixed scientific method

- Public checkpoint: `facebook/esm2_t33_650M_UR50D`, pinned to revision
  `08e4846e537177426273712802403f7ba8261b6c`.
- The same `code/generate_esm2_embeddings.py` is used on CPU and CUDA.
- Source: `data/proteins_ec3.tsv`; only `accession` and `sequence` are loaded
  as model inputs. All 146,677 unique accessions are retained, including the
  two repeated sequences under distinct accessions.
- Sequences are used unchanged. Lengths must be 1--1,022 residues; empty,
  invalid, overlength, or unsupported-token inputs stop the run. No truncation.
- Frozen `AutoModel` encoder with pooling head disabled, `model.eval()`,
  `requires_grad_(False)`, and `torch.inference_mode()`.
- `last_hidden_state` is mean-pooled over amino-acid residue tokens only,
  using attention and special-token masks. CLS/BOS, EOS, and padding are
  excluded. Residue token counts and identities are checked against input.
- Dimension must be 1,280 in both checkpoint configuration and model output.
- No L2 normalization is applied. Saved vectors are always float32.
- Matrix rows follow lexicographic accession order. Scheduling by sequence
  length improves padding efficiency but never changes matrix row indexing.
- Each accession has one embedding in the global library. Later experiments
  join by accession; embeddings are not recomputed per held-out phylum.
- Existing data, split assignments, nearest identities, and similarity strata
  are hashed and remain unchanged.

## GPU resources and precision

Recommend **at least 16 GB GPU VRAM** and **16 GB host RAM** for a comfortable
starting point. An 8 GB GPU may handle batch size 1, but this has not been
measured here and is not a guaranteed memory bound. Actual peak memory depends
on length, batch size, PyTorch/CUDA build, and inference precision. A 24 GB GPU
provides more room for larger batches; it is not required by the method.

Start with `--batch-size 1`. Validate the CUDA pilot, which spans input lengths
and includes a 1,022-residue protein. If desired, test batch sizes 2 and 4 using
fresh pilot output names before selecting the full-run batch size. Retain VRAM
headroom for the longest padded batch. The script's default auto batch is 1,
or 4 when CUDA reports at least 12 GiB free; explicit size 1 is the most
conservative handoff command. No runtime OOM is silently ignored: completed
rows remain resumable, and the failing batch is not marked complete.

CPU uses float32 inference. With `--precision float16` on CUDA, model weights
remain float32 and eligible operations use CUDA float16 autocast. Pooling and
saved vectors are explicitly float32. Use `--precision float32` instead if
desired, but use the same precision for the CUDA pilot and full run. Precision,
resolved model revision, and runtime versions are saved. Mixed precision can
produce small differences from CPU float32; independently recomputed sample
rows must pass the documented numerical tolerances. Resume rejects changes to
device type, precision, model revision, source, and PyTorch/Transformers versions.

## Expected storage

The full matrix requires **750,986,368 bytes**, approximately **751 MB / 716 MiB**,
including the NumPy header, for shape `(146677, 1280)` and dtype float32.
Allow additional tens of MB for accession metadata and the SQLite checkpoint
ledger. Checkpoint weights download to approximately 2.6 GB, with small extra
tokenizer/config files. Cache behavior can require additional space on Windows
when symlinks are unavailable. Keep **15--20 GB free** for weights, download
temporaries, data, results, and environment overhead; a new CUDA software
installation or complete repository clone may need more. MMseqs2 temporary
databases and preflight artifacts are not required for generating embeddings.

## Environment on the CUDA machine

Recreate the Python environment; do not copy the Windows `.venv` to another
machine. Use a CUDA-enabled PyTorch 2.11 build compatible with the GPU driver,
plus `code/requirements_esm2.txt`. A CPU-only PyTorch installation does not
become CUDA-capable merely because a GPU is present. If cloning this repository,
obtain the actual cleaned TSV through Git LFS rather than its pointer file.
Keep the existing frozen analysis and split metadata for provenance checks.

From the project root, verify CUDA before running:

```bash
python -c "import torch; assert torch.cuda.is_available(), 'CUDA unavailable'; print(torch.__version__, torch.version.cuda, torch.cuda.get_device_name(0), torch.cuda.mem_get_info())"
python code/test_generate_esm2_embeddings.py
```

Here `python` means the CUDA machine's configured Python environment. Replace
it with `.venv/bin/python` or `.venv\Scripts\python.exe` as appropriate.

## Exact CUDA commands

These output names deliberately differ from the local CPU pilot, even if the
workspace is copied rather than cloned. The model and revision are explicit.

First run the CUDA pilot:

```bash
python code/generate_esm2_embeddings.py --pilot --device cuda --precision float16 --batch-size 1 --revision 08e4846e537177426273712802403f7ba8261b6c --output analysis/esm2_embeddings_cuda
```

Verify `analysis/esm2_embeddings_cuda_pilot/validation.json` has
`"all_passed": true`. Only then start the full job:

```bash
python code/generate_esm2_embeddings.py --device cuda --precision float16 --batch-size 1 --revision 08e4846e537177426273712802403f7ba8261b6c --output analysis/esm2_embeddings_cuda
```

Resume that same full job after interruption:

```bash
python code/generate_esm2_embeddings.py --resume --device cuda --precision float16 --batch-size 1 --revision 08e4846e537177426273712802403f7ba8261b6c --output analysis/esm2_embeddings_cuda
```

Do not resume CPU pilot vectors into a CUDA run. Both pilots and the full run
use the same method, but separate output directories prevent mixing numerical
settings or treating a pilot subset as the full embedding library.

## Completion, resume, and downstream loading

The full output must have `progress.json` with `status=complete` and
`validation.json` with `all_passed=true`, `input_proteins=146677`, and shape
`[146677, 1280]` before downstream use. A preallocated matrix alone does not
prove completion. `accessions.tsv` provides the explicit zero-based row index.

`checkpoint.sqlite3` commits each row only after the NumPy matrix is flushed
and synced. Resume checks source/row mapping and saved vector checksums before
skipping completed work. Existing outputs require `--resume`; `--restart`
creates a new timestamped sibling rather than deleting the previous run.
Final validation checks every row for finiteness, every completed vector hash,
unexpected duplicate vectors, matrix statistics, frozen-input hashes, and
random independent single-protein recomputations plus the longest protein.
Identical vectors for identical source sequences are legitimate.

Example loading after full validation:

```python
import numpy as np
import pandas as pd

matrix = np.load("analysis/esm2_embeddings_cuda/embeddings.npy", mmap_mode="r")
mapping = pd.read_csv("analysis/esm2_embeddings_cuda/accessions.tsv", sep="\t")
row_for_accession = mapping.set_index("accession")["row_index"]
# Join frozen experimental split records by accession; never assume split order.
```

Large caches and generated embedding output directories are ignored by Git.
Archive complete results and metadata intentionally after validation; do not
commit partial matrices as completed scientific artifacts.

## Validated local CPU pilot results

The local pilot completed and passed all validation checks. **The full
146,677-protein embedding job was not started on this machine.**

- Model: facebook/esm2_t33_650M_UR50D; revision
  08e4846e537177426273712802403f7ba8261b6c.
- Eight proteins spanning 34--1,022 residues, deliberately selected by length.
- Matrix: (8, 1280), float32 NumPy .npy, **41,088 bytes**.
- Device: CPU; float32 inference; batch size 1; four PyTorch threads.
- Initial eight-protein generation: **94.27 seconds**, excluding
  model download/loading and additional independent validation. This is not a
  measured GPU runtime estimate.
- Three independent accession-matched recomputations, including the longest
  protein, gave maximum absolute difference **0.0**.
- Offline resume loaded all 8 completed rows, regenerated no missing vectors,
  and retained the original matrix SHA256:
  7c0083e8f14e6e46f51dfbb3311b87704ef3fb042ef752bd463d32d247088de0.
- No failed proteins, truncation, NaN/Inf, or unexpected duplicate vectors.
- Embedding value statistics: min -8.44131184, max 0.39688271,
  mean -0.000672599, population standard deviation 0.198374249.
- All original scientific input hashes and the executed script hash verified.
  Special-token exclusion, mapping and actual 1,022-residue tokenization were
  independently audited. Four regression tests passed.
- Python 3.13.12, PyTorch 2.11.0+cpu,
  Transformers 4.57.6; CUDA unavailable.
- Public weights were downloaded once into the local cache (2,609,506,392 bytes).
  Initial network interruptions resumed successfully. Windows does not support
  symlinks in this cache configuration; the cache still works. The pinned
  Transformers version emits a harmless torch_dtype deprecation warning.
- Pilot artifacts under analysis/esm2_embeddings_pilot/: embeddings.npy,
  accessions.tsv, checkpoint.sqlite3, settings.json, progress.json,
  validation.json, pilot_handoff_validation.json, run.log, and report.md.

The full cleaned dataset was validated as 146,677 unique accessions and
146,675 distinct sequences, maximum length 1,022. Identical sequences under
different accessions remain distinct indexed records. No EC3 labels were
used to generate embeddings; no ESM-2 fine-tuning, classifier training or
test-set performance inspection was performed. The existing split and
similarity-stratum assignments and nearest-training identities are unchanged.
