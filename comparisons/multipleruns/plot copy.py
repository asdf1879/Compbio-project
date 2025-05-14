import numpy as np
import matplotlib.pyplot as plt
from math import comb
import pyfastx
import os
import csv
import re
from sklearn.metrics import auc

# --- Setup paths ---
base_dir = os.path.dirname(__file__)
fasta_path = os.path.join(base_dir, "NCTC1080_reads.fasta.gz")
gt_path = os.path.join(base_dir, "NCTC1080_daligner_ground_truth.txt")

files = [
  
    ("aln_lh_500hkmax32.tsv", "500hkmax32", "tab:blue"),
    ("merged_aln.tsv", "merged_aln", "tab:orange"),

]


# --- Load reads and sequence lengths ---
seqs = list(pyfastx.Fastx(fasta_path))
n_seqs = len(seqs)
map_len = {i: len(seq[1]) for i, seq in enumerate(seqs)}

# --- Load and filter ground truth ---
true_pairs = set()
with open(gt_path) as f:
    for line in f:
        parts = line.strip().split()
        a, b, aln = int(parts[0]), int(parts[1]), float(parts[2])
        l1, l2 = map_len[a], map_len[b]
        jaccard = aln / (l1 + l2 - aln)
        if jaccard >= 0.2:
            true_pairs.add((min(a, b), max(a, b)))

# --- Evaluate from CSV ---
def evaluate_csv(file_path):
    predicted = {}
    with open(file_path) as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            if len(row) >= 5:
                a, b = int(row[0]), int(row[1])
                score = float(row[4])
                predicted[(min(a, b), max(a, b))] = score
    return evaluate_common(predicted)

# --- Evaluate from TSV ---
def evaluate_tsv(file_path):
    predicted = {}
    with open(file_path) as f:
        for line in f:
            parts = re.split(r"\s+", line.strip())
            if len(parts) >= 3:
                a, b, score = int(parts[0])-1, int(parts[1])-1, float(parts[2])
                predicted[(min(a, b), max(a, b))] = score
    return evaluate_common(predicted)

# --- Common evaluation logic ---
def evaluate_common(predicted):
    sorted_preds = sorted(predicted.items(), key=lambda x: -x[1])
    labels = [1 if p in true_pairs else 0 for p, _ in sorted_preds]

    if len(labels) < 10:
        return np.array([]), np.array([]), np.array([]), np.array([]), 0.0, 0.0

    tp = np.cumsum(labels)
    fp = np.cumsum([1 - x for x in labels])
    total_true = len(true_pairs)
    total_pairs = comb(n_seqs, 2)
    total_neg = total_pairs - total_true

    Ns = np.arange(10, len(labels), 10)
    TP = tp[Ns]
    FP = fp[Ns]
    FN = total_true - TP
    TN = total_neg - FP

    recall = TP / total_true
    precision = TP / (TP + FP)
    tpr = TP / total_true
    fpr = FP / (FP + TN)

    pr_auc = auc(recall, precision)
    roc_auc = auc(fpr, tpr)

    return recall, precision, fpr, tpr, pr_auc, roc_auc

# --- Plot Precision-Recall Curve ---
plt.figure(figsize=(7, 5))
for fname, label, color in files:
    full_path = os.path.join(base_dir, fname)
    
    rec, prec, _, _, pr_auc, _ = evaluate_tsv(full_path)

    if len(rec) == 0:
        print(f"Skipping {label}: not enough data")
        continue

    plt.plot(rec, prec, label=f"{label} (AUC={pr_auc:.2f})", color=color)

plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Precision-Recall Curve")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(base_dir, "prc.png"))
plt.show()

# --- Plot ROC Curve ---
plt.figure(figsize=(7, 5))
for fname, label, color in files:
    full_path = os.path.join(base_dir, fname)
    if fname.endswith(".csv"):
        _, _, fpr, tpr, _, roc_auc = evaluate_csv(full_path)
    else:
        _, _, fpr, tpr, _, roc_auc = evaluate_tsv(full_path)

    if len(fpr) == 0:
        print(f"Skipping {label}: not enough data")
        continue

    plt.plot(fpr, tpr, label=f"{label} (AUC={roc_auc:.2f})", color=color)

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(base_dir, "roc.png"))
plt.show()
