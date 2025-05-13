import os
import re
import numpy as np
import matplotlib.pyplot as plt
from math import comb
import pyfastx
from sklearn.metrics import auc

# --- Load sequence lengths ---
seqs = list(pyfastx.Fastx("NCTC1080_reads.fasta.gz"))
n_seqs = len(seqs)
map_len = {i: len(seq[1]) for i, seq in enumerate(seqs)}

# --- Load ground truth ---
ground_truth_raw = {}
with open("NCTC1080_daligner_ground_truth.txt", "r") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 3:
            a, b, aln_len = int(parts[0]), int(parts[1]), float(parts[2])
            pair = (min(a, b), max(a, b))
            ground_truth_raw[pair] = aln_len

# --- Filter ground truth with Jaccard threshold ---
jaccard_threshold = 0.2
true_pairs = set()
for (i1, i2), a in ground_truth_raw.items():
    l1, l2 = map_len.get(i1, 0), map_len.get(i2, 0)
    denom = l1 + l2 - a
    jaccard = a / denom if denom > 0 else 0
    if jaccard >= jaccard_threshold:
        true_pairs.add((i1, i2))

# --- Setup groups ---
base_dir = "outputs"
group_patterns = {
    "kmax32": [],
    "kmax60": [],
    "top2k": []
}

for root, _, files in os.walk(base_dir):
    for fname in files:
        if not fname.endswith(".tsv"):
            continue
        full_path = os.path.join(root, fname)
        if "kmax32" in fname:
            group_patterns["kmax32"].append(full_path)
        elif "kmax60" in fname:
            group_patterns["kmax60"].append(full_path)
        elif "top2k" in fname.lower():
            group_patterns["top2k"].append(full_path)

# --- Plot PRC and ROC per group ---
def evaluate_and_plot(group_name, file_list):
    colors = plt.cm.get_cmap("tab10", len(file_list))  # Get distinct colors from the "tab10" colormap
    plt.figure(figsize=(8,6))
    for idx, path in enumerate(file_list):
        predicted_map = {}
        with open(path, "r") as f:
            for line in f:
                parts = re.split(r"\s+", line.strip())
                if len(parts) >= 4:
                    q, t, l = int(parts[0])-1, int(parts[1])-1, int(parts[2])
                    pair = (min(q, t), max(q, t))
                    predicted_map[pair] = l

        sorted_predicted = sorted(predicted_map.items(), key=lambda x: -x[1])
        labels = [1 if pair in true_pairs else 0 for pair, _ in sorted_predicted]

        tp_cumsum = np.cumsum(labels)
        fp_cumsum = np.cumsum([1 - l for l in labels])
        total_true = len(true_pairs)
        total_pairs = comb(n_seqs, 2)
        total_neg = total_pairs - total_true

        Ns = np.arange(10, len(labels) + 1, 10)
        TP = tp_cumsum[Ns - 1]
        FP = fp_cumsum[Ns - 1]
        FN = total_true - TP
        TN = total_neg - FP

        precisions = TP / (TP + FP)
        recalls = TP / (TP + FN)
        fprs = FP / (FP + TN)
        tprs = recalls

        # Calculate AUC-PR
        pr_auc = auc(recalls, precisions)
        
        label = os.path.basename(path).split(".")[0]
        plt.plot(recalls, precisions, label=f"{label} (AUC-PR: {pr_auc:.2f})", color=colors(idx))

    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(f"Precision-Recall Curve: {group_name}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{group_name}_prc.png")
    plt.show()

    # ROC
    plt.figure(figsize=(8,6))
    for idx, path in enumerate(file_list):
        predicted_map = {}
        with open(path, "r") as f:
            for line in f:
                parts = re.split(r"\s+", line.strip())
                if len(parts) >= 4:
                    q, t, l = int(parts[0])-1, int(parts[1])-1, int(parts[2])
                    pair = (min(q, t), max(q, t))
                    predicted_map[pair] = l

        sorted_predicted = sorted(predicted_map.items(), key=lambda x: -x[1])
        labels = [1 if pair in true_pairs else 0 for pair, _ in sorted_predicted]

        tp_cumsum = np.cumsum(labels)
        fp_cumsum = np.cumsum([1 - l for l in labels])
        total_true = len(true_pairs)
        total_pairs = comb(n_seqs, 2)
        total_neg = total_pairs - total_true

        Ns = np.arange(10, len(labels) + 1, 10)
        TP = tp_cumsum[Ns - 1]
        FP = fp_cumsum[Ns - 1]
        FN = total_true - TP
        TN = total_neg - FP

        fprs = FP / (FP + TN)
        tprs = TP / (TP + FN)

        # Calculate AUC-ROC
        roc_auc = auc(fprs, tprs)

        label = os.path.basename(path).split(".")[0]
        plt.plot(fprs, tprs, label=f"{label} (AUC-ROC: {roc_auc:.2f})", color=colors(idx))

    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC Curve: {group_name}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{group_name}_roc.png")
    plt.show()

# Run for each group
for group_name, paths in group_patterns.items():
    if paths:
        evaluate_and_plot(group_name, paths)
