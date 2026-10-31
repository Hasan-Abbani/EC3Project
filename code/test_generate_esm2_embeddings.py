"""Regression checks for residue-only pooling and resumable row validation."""
import hashlib
from pathlib import Path
import sqlite3
import tempfile
import unittest

import numpy as np
import pandas as pd
import torch

from generate_esm2_embeddings import residue_mean, load_proteins, check_ledger


class EmbeddingTests(unittest.TestCase):
    def test_pooling_excludes_bos_eos_and_padding(self):
        hidden = torch.tensor([[[100.,100.],[1.,3.],[3.,5.],[200.,200.],[999.,999.]]])
        attention = torch.tensor([[1,1,1,1,0]])
        special = torch.tensor([[1,0,0,1,1]])
        actual = residue_mean(hidden,attention,special,[2])
        torch.testing.assert_close(actual,torch.tensor([[2.,4.]]))

    def test_pooling_rejects_token_length_mismatch(self):
        with self.assertRaises(ValueError):
            residue_mean(torch.zeros(1,4,2),torch.ones(1,4),torch.tensor([[1,0,0,1]]),[3])

    def test_accession_order_label_independence_and_overlength(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"input.tsv"
            pd.DataFrame(dict(accession=["B","A"],sequence=["ACD","ME"],ec3=["unused","unused"])).to_csv(path,sep="\t",index=False)
            first = load_proteins(path)
            self.assertEqual(first.accession.tolist(),["A","B"])
            pd.DataFrame(dict(accession=["A","B"],sequence=["ME","ACD"],ec3=["changed","changed"])).to_csv(path,sep="\t",index=False)
            pd.testing.assert_frame_equal(first,load_proteins(path))
            pd.DataFrame(dict(accession=["A"],sequence=["A"*1023])).to_csv(path,sep="\t",index=False)
            with self.assertRaisesRegex(ValueError,"Overlength"):
                load_proteins(path)

    def test_resume_validates_vector_and_accession(self):
        frame = pd.DataFrame(dict(accession=["A","B"],sequence_sha256=["seqA","seqB"]))
        matrix = np.zeros((2,1280),dtype=np.float32)
        database = sqlite3.connect(":memory:")
        database.execute("CREATE TABLE completed (row_index INTEGER,accession TEXT,sequence_sha256 TEXT,vector_sha256 TEXT)")
        row_hash = hashlib.sha256(matrix[0].tobytes()).hexdigest()
        database.execute("INSERT INTO completed VALUES (0,'A','seqA',?)",(row_hash,))
        self.assertEqual(check_ledger(database,matrix,frame),{0})
        matrix[0,0] = 1
        with self.assertRaisesRegex(ValueError,"checksum"):
            check_ledger(database,matrix,frame)
        database.close()


if __name__ == "__main__":
    unittest.main()
