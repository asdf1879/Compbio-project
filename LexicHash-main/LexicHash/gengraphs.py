from utils import metric_utils
from utils.metric_utils import get_gt_df, get_gt_df_nctc , get_overlaps, get_pred_df ,get_thetas, roc , roc_and_prc, prc, pr_vs_k
from utils.seq_utils import get_fasta, get_fastq, get_seqs, write_fasta
import numpy as np
import pyfastx

gt_path = 'NCTC1080_daligner_ground_truth.txt'

fasta_file_path = 'NCTC1080_reads.fasta.gz'
pred_path1 = 'aln_lh_500h.tsv'
seqs = pyfastx.Fastx(fasta_file_path)
seqs = [seq[1] for seq in seqs]
# print('seqs', seqs)
print(len(seqs))
seq_lens = np.array([len(seq) for seq in seqs])

print('seq_lens', seq_lens)

gt_df = get_gt_df( gt_path, seq_lens, 0 )

#store to csv
gt_df.to_csv('gt_df.csv', index=False)


print(gt_df.head())

pred_df = get_pred_df( pred_path1 , one_idx=True )

print(pred_df.head())

names = ["Lexichash"]

pred_dfs= [pred_df]

roc_and_prc(pred_dfs, gt_df, names, "NCTC1080", n_seq=len(seqs),rc=True)


