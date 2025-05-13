import re
import pyfastx
import numpy as np
import matplotlib.pyplot as plt
from math import comb

# --- Load sequences and compute lengths ---
seqs = list(pyfastx.Fastx("NCTC1080_reads.fasta.gz"))
n_seqs = len(seqs)
map_len = {i: len(seq[1]) for i, seq in enumerate(seqs)}

for i in range(n_seqs):
    print(f"seq {i}: {map_len[i]}")

# --- Load ground truth alignments ---
ground_truth_raw = {}
with open("NCTC1080_daligner_ground_truth.txt", "r") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 3:
            a, b, aln_len = int(parts[0]), int(parts[1]), float(parts[2])
            pair = (min(a, b), max(a, b))
            ground_truth_raw[pair] = aln_len

# --- Load predictions ---
predicted_map = {}
with open("pairwise_results.txt", "r") as f:
    for line in f:
        parts = re.split(r"\s+", line.strip())
        if len(parts) >= 4:
            q, t, l, s = parts[:4]
            q, t = int(q)-1, int(t)-1
            pair = (min(q, t), max(q, t))
            predicted_len = int(l)
            predicted_map[pair] = predicted_len

# --- Filter ground truth by Jaccard threshold ---
jaccard_threshold = 0.2
true_pairs = set()
for (i1, i2), a in ground_truth_raw.items():
    l1 = map_len.get(i1, 0)
    l2 = map_len.get(i2, 0)
    denom = l1 + l2 - a
    jaccard = a / denom if denom > 0 else 0
    if jaccard >= jaccard_threshold:
        true_pairs.add((i1, i2))

# --- Sort predicted scores descending ---
sorted_predicted = sorted(predicted_map.items(), key=lambda x: -x[1])
labels = [1 if pair in true_pairs else 0 for pair, _ in sorted_predicted]

# --- Cumulative counts ---
tp_cumsum = np.cumsum(labels)
fp_cumsum = np.cumsum([1 - l for l in labels])
total_true = len(true_pairs)
total_pairs = comb(n_seqs, 2)
total_neg = total_pairs - total_true

# --- Sample every 10 predictions ---
Ns = np.arange(10, len(labels) + 1, 10)
TP = tp_cumsum[Ns - 1]
FP = fp_cumsum[Ns - 1]
FN = total_true - TP
TN = total_neg - FP

# --- Compute metrics ---
precisions = TP / (TP + FP)
recalls = TP / (TP + FN)
fprs = FP / (FP + TN)
tprs = recalls  # Same as TPR

# --- Plot PR Curve ---
plt.figure(figsize=(8,6))
plt.plot(recalls, precisions, marker='o', linewidth=1.5, markersize=2)
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Precision-Recall Curve (Every 10 Predictions)")
plt.grid(True)
plt.tight_layout()
plt.savefig("precision_recall_every_10.png")
plt.show()

# --- Plot ROC Curve ---
plt.figure(figsize=(8,6))
plt.plot(fprs, tprs, marker='o', linewidth=1.5, markersize=2)
plt.xlabel("False Positive Rate (FPR)")
plt.ylabel("True Positive Rate (TPR)")
plt.title("ROC Curve (Every 10 Predictions)")
plt.grid(True)
plt.tight_layout()
plt.savefig("roc_every_10.png")
plt.show()
