import unittest,sys,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from model import design,synthetic_sequence,read_fasta
class Tests(unittest.TestCase):
    def test_primer_binding_and_product(self):
        s=synthetic_sequence();r=design(s,1000)
        self.assertEqual(r['pairs_found'],1)
        l=r['left_primer'];rr=r['right_primer'];end=r['right_end']
        self.assertEqual(s[r['left_start']-1:r['left_start']-1+len(l)],l)
        self.assertEqual(s[end-len(rr):end].translate(str.maketrans('ACGT','TGCA'))[::-1],rr)
        self.assertEqual(end-r['left_start']+1,r['product_size'])
        self.assertTrue(600<=r['product_size']<=800)
        self.assertLess(r['left_start']+len(l)-1,r['target_start'])
        self.assertGreater(end-len(rr)+1,r['target_start']+r['target_length']-1)
    def test_boundaries_and_no_pair(self):
        for pos in [1,2000]:
            r=design(synthetic_sequence(),pos);self.assertGreaterEqual(r['window_start'],1);self.assertLessEqual(r['window_end'],2000)
        self.assertEqual(design('A'*100,50)['pairs_found'],0)
        for seq,pos in [('ACGN',2),('ACGT',0),('ACGT',5)]:
            with self.assertRaises(ValueError):design(seq,pos)
    def test_fasta(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.fa';p.write_text('>a description\nAC\nGT\n>b\nTT\n');self.assertEqual(read_fasta(p,'a'),'ACGT')
            with self.assertRaises(ValueError):read_fasta(p,'c')
if __name__=='__main__':unittest.main()
