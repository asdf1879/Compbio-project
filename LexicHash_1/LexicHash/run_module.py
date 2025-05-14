#!/home/gcgreen2/envs/py3.9venv/bin/python3
import sys
import os
import argparse
from os.path import join
import pandas as pd

from LexicHash import LexicHash, MinHash
from LexicHash.utils import utils


def check_args(args):
    assert args['method'] in ['lexichash', 'minhash'], 'method argument must be either "lexichash" or "minhash"'
    if args['method'] == 'minhash':
        assert args['k'] is not None, 'must set a k-value for minhash'
        assert 7 <= args['k'] <= 24, 'k-value must be between 7 and 24'
    assert args['n_cpu'] <= os.cpu_count(), f'only {os.cpu_count()} cpus available'


def add_args(args, stage_suffix):
    run_id = utils.get_run_id(**args)
    args['n_bits'] = 30
    args['aln_path'] = join(args['out_dir'], f'aln_{stage_suffix}.tsv')
    args['sketch_path'] = join(args['out_dir'], 'tmp', f'sketches_{run_id}.npy')
    args['args_path'] = join(args['out_dir'], 'tmp', f'args_{run_id}.npy')
    return args


def parse_args():
    parser = argparse.ArgumentParser(description='Sequence Similarity Estimation via Lexicographic Comparison of Hashes')
    parser.add_argument('--out', type=str, dest='out_dir', required=True, help='output directory')
    parser.add_argument('--fasta', type=str, dest='fasta_path', required=True, help='path to reads fasta')
    parser.add_argument('--method', type=str, default='lexichash', help='alignment method (lexichash/minhash)')
    parser.add_argument('--n_hash', type=int, default=500, help='number of hash functions to use')
    parser.add_argument('--k', type=int, default=16, help='k-value (for minhash only)')
    parser.add_argument('--min_k', type=int, default=14, help='min match length (for lexichash only)')
    parser.add_argument('--max_k', type=int, default=32, help='max match length (for lexichash only)')
    parser.add_argument('--no_rc', action='store_false', dest='rc', help='do not account for reverse-complements')
    parser.add_argument('--min_n_col', type=int, default=1, help='min # of minhash collisions (for minhash only)')
    parser.add_argument('--n_cpu', type=int, default=os.cpu_count(), help='number of cpus to use')
    parser.add_argument('--alpha', type=float, default=None, help='alpha*n_seq pairs will be output')
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

    # Mapping of (min_k, max_k) to (n_hash, alpha)
    stages = [
        ((63, 128), 700, 10.0),
        ((31, 64), 600, 100.0),
        ((14, 32), 600, 1000.0),
        ((13, 16), 400, 10.0)
    ]

    merged_paths = []

    for (min_k, max_k), n_hash, alpha in stages:
        print(f"\n[INFO] Running LexicHash with min_k = {min_k}, max_k = {max_k}, n_hash = {n_hash}, alpha = {alpha}")
        args['min_k'] = min_k
        args['max_k'] = max_k
        args['n_hash'] = n_hash
        args['alpha'] = alpha

        stage_suffix = f'k{max_k}'
        args = add_args(args, stage_suffix)
        check_args(args)
        utils.setup_out_dir(**args)
        utils.print_header(**args)

        if args['method'] == 'lexichash':
            LexicHash.find_overlaps(**args)
        elif args['method'] == 'minhash':
            MinHash.find_overlaps(**args)

        merged_paths.append(args['aln_path'])

    # Merge all output files into one
    merged_df = merge_prediction_files(merged_paths)
    merged_output_path = join(args['out_dir'], 'merged_aln.tsv')
    merged_df.to_csv(merged_output_path, sep="\t", index=False, header=False)
    print(f"\n[INFO] Merged output written to: {merged_output_path}")


if __name__ == '__main__':
    main()
