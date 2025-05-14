

import pandas as pd
from metric_utils import get_gt_df_nctc, get_pred_df, roc_and_prc

# Path to prediction file
pred_path = r"C:\Users\panum\Downloads\LH_out\merged_aln.tsv"

# Load your prediction
pred_df = get_pred_df(pred_path)

# Sequence lengths (dummy example: length of each read)
# You need the actual lengths in the correct order
seq_lens = pd.read_csv(r"C:\Users\panum\Downloads\seq_lens.txt", header=None)[0].values

# Ground truth file (must be provided or generated)
gt_df = get_gt_df_nctc(r"C:\Users\panum\Downloads\NCTC1080_daligner_ground_truth (1).txt", seq_lens)

# Run metrics
aurocs, fprs, tprs, auprcs, precs, recalls = roc_and_prc(
    [pred_df],     # List of predictions
    gt_df,         # Ground truth
    ["LH_out"],    # Names for the legend
    "MyDataset",   # Dataset name for title
    len(seq_lens)  # Number of sequences
)
