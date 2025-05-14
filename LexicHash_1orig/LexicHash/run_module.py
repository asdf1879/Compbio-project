
#!/home/gcgreen2/envs/py3.9venv/bin/python3
import sys
import os
import argparse
from os.path import join
import pandas as pd

from LexicHash import LexicHash, MinHash
from LexicHash.utils import utils

sys.argv = [
  "run_module.py",
  "--fasta", "data/NCTC1080/NCTC1080_reads.fasta.gz",
  "--out", "../LH_out",
  "--n_hash", "500",
  "--n_cpu", "8"
]


def check_args(args):
    assert args['method'] in ['lexichash', 'minhash'], 'method argument must be either "lexichash" or "minhash"'
    if args['method'] == 'minhash':
        assert args['k'] is not None, 'must set a k-value for minhash'
        assert args['k'] >= 7 and args['k'] <= 24, 'k-value must be between 7 and 24'
    assert args['n_cpu'] <= os.cpu_count(), f'only {os.cpu_count()} cpus available'

def add_args(args, stage_suffix):
    run_id = utils.get_run_id(**args)
    args['n_bits'] = 30
    args['aln_path'] = join(args['out_dir'], f'aln_{stage_suffix}.tsv')  # modified
    args['sketch_path'] = join(args['out_dir'], 'tmp', f'sketches_{run_id}.npy')
    args['args_path'] = join(args['out_dir'], 'tmp', f'args_{run_id}.npy')
    return args

def parse_args():
    parser = argparse.ArgumentParser(description='Sequence Similarity Estimation via Lexicographic Comparison of Hashes')
    parser.add_argument('--out', type=str, dest='out_dir', help='output directory', required=True)
    parser.add_argument('--fasta', type=str, dest='fasta_path', help='path to reads fasta', required=True)
    parser.add_argument('--method', type=str, help='alignment method (lexichash/minhash)', default='lexichash')
    parser.add_argument('--n_hash', type=int, default=500, help='number of hash functions to use')
    parser.add_argument('--k', type=int, help='k-value (for minhash only)', default=16)
    parser.add_argument('--min_k', type=int, help='min match length (for lexichash only)', default=14)
    parser.add_argument('--max_k', type=int, help='max match length (for lexichash only)', default=32)
    parser.add_argument('--no_rc', action='store_false', help='do not account for reverse-complements', dest='rc')
    parser.add_argument('--min_n_col', type=int, help='min # of minhash collisions (for minhash only)', default=1)
    parser.add_argument('--n_cpu', type=int, help='number of cpus to use', default=os.cpu_count())
    parser.add_argument('--alpha', type=float, help='alpha*n_seq pairs will be output', default=None)
    parser.set_defaults(rc=True)
    return vars(parser.parse_args())

def merge_prediction_files(paths):
    merged = {}
    for path in paths:
        df = pd.read_csv(path, sep="\t", header=None, names=["i1", "i2", "score", "direction"])
        for _, row in df.iterrows():
            key = (int(row.i1), int(row.i2), row.direction)
            score = int(row.score)
            if key not in merged or score > merged[key]:
                merged[key] = score
    merged_df = pd.DataFrame([(i1, i2, score, d) for (i1, i2, d), score in merged.items()],
                             columns=["i1", "i2", "score", "direction"])
    return merged_df

def main():
    args = parse_args()
    k_values = [(127, 256), (63, 128), (31, 64), (14, 32)]

    for min_k, max_k in k_values:
        print(f"Running LexicHash with min_k = {min_k} and max_k = {max_k}")
        args['min_k'] = min_k
        args['max_k'] = max_k
        stage_suffix = f'k{max_k}'
        args = add_args(args, stage_suffix)
        check_args(args)
        utils.setup_out_dir(**args)
        utils.print_header(**args)

        if args['method'] == 'lexichash':
            LexicHash.find_overlaps(**args)
        elif args['method'] == 'minhash':
            MinHash.find_overlaps(**args)

    # Merge results from all kmax stages -- added
    merged_paths = [join(args['out_dir'], f'aln_k{K}.tsv') for K in [256, 128, 64, 32]]
    merged_df = merge_prediction_files(merged_paths)
    merged_output_path = join(args['out_dir'], 'merged_aln.tsv')
    merged_df.to_csv(merged_output_path, sep="\t", index=False, header=False)
    print(f"[INFO] Merged output written to: {merged_output_path}")


# import sys
# import os


# from LexicHash.utils.metric_utils import get_gt_df_nctc, get_pred_df, roc_and_prc


# # Path to prediction file
# pred_path = r"C:\Users\panum\Downloads\LH_out\merged_aln.tsv"

# # Load your prediction
# pred_df = get_pred_df(pred_path)

# # Sequence lengths 
# seq_lens = pd.read_csv(r"C:\Users\panum\Downloads\seq_lens.txt", header=None)[0].values

# # Ground truth file 
# gt_df = get_gt_df_nctc(r"C:\Users\panum\Downloads\NCTC1080_daligner_ground_truth (1).txt", seq_lens)

# # Run metrics
# aurocs, fprs, tprs, auprcs, precs, recalls = roc_and_prc(
#     [pred_df],     # List of predictions
#     gt_df,         # Ground truth
#     ["LH_out"],    # Names for the legend
#     "MyDataset",   # Dataset name for title
#     len(seq_lens)  # Number of sequences



if __name__ == '__main__':
    main()

