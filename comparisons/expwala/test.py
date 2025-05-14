import matplotlib.pyplot as plt
import re
import pyfastx
import os
import glob
import csv
def percent_top100(predicted, ground_truth, k):
    sorted_ground_truth = sorted(ground_truth.items(), key=lambda x: x[1], reverse=True)
    ground_truth_pairs = {pair: score for pair, score in sorted_ground_truth[:k]}
    percent_overlap = 0

    pred = sorted(predicted.items(), key=lambda x: x[1], reverse=True)
    predicted = {pair: score for pair, score in pred[:k]}

    for key in predicted:
        if key in ground_truth_pairs:
            percent_overlap += 1

    return percent_overlap / k * 100

def above_mean(predicted_map, mean):
    cnt = 0
    for key in predicted_map:
        true_score, _ = predicted_map[key]
        if true_score > mean:
            cnt += 1
    return (cnt / len(predicted_map)) * 100

# Load sequences
seqs = list(pyfastx.Fastx("NCTC1080_reads.fasta.gz"))
seqs = [seq[1] for seq in seqs]

# Ground truth
ground_truth = {}
list_of_scores = []

with open("NCTC1080_daligner_ground_truth.txt", "r") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 3:
            a, b, aln_len = int(parts[0]), int(parts[1]), float(parts[2])
            key = (min(a, b), max(a, b))
            list_of_scores.append(aln_len)
            ground_truth[key] = aln_len

# Plot GT distribution once
plt.hist(list_of_scores, bins=100)
plt.xlabel('Alignment Score')
plt.ylabel('Frequency')
plt.title('Distribution of Ground Truth Alignment Scores')
plt.savefig('alignment_scores_distribution.png')
plt.close()

sorted_scores = sorted(list_of_scores)
mean_score = sum(sorted_scores) / len(sorted_scores)
median_score = sorted_scores[len(sorted_scores) // 2]
min_score = sorted_scores[0]
max_score = sorted_scores[-1]

print(f"Ground truth mean: {mean_score:.2f}")
print(f"Ground truth median: {median_score:.2f}")
print(f"Ground truth min: {min_score:.2f}")
print(f"Ground truth max: {max_score:.2f}")

# Process each predicted .tsv file

tsv_files = glob.glob("*.tsv")
os.makedirs("results", exist_ok=True)

for tsv_path in tsv_files:
    k = 2000
    matches = []
    not_in_ground_truth = []
    predicted_map = {}
    in_ground_truth = {}
    list_of_similar_predictions = []

    with open(tsv_path, "r") as f:
        for line in f:
            parts = re.split(r"\s+", line.strip())
            if len(parts) >= 4:
                q, t, l, s = parts[:4]
                q, t = int(q)-1, int(t)-1
                pair = (min(q, t), max(q, t))
                predicted_len = int(l)
                predicted_map[pair] = predicted_len

                if pair in ground_truth:
                    true_len = ground_truth[pair]
                    matches.append((q, t, predicted_len, true_len))
                    list_of_similar_predictions.append(true_len)
                    in_ground_truth[pair] = [true_len, predicted_len]
                else:
                    not_in_ground_truth.append((q, t, predicted_len))

    if not list_of_similar_predictions:
        print(f"{tsv_path}: No overlaps found.")
        continue

    mean_pred = sum(list_of_similar_predictions) / len(list_of_similar_predictions)
    median_pred = sorted(list_of_similar_predictions)[len(list_of_similar_predictions) // 2]
    min_pred = min(list_of_similar_predictions)
    max_pred = max(list_of_similar_predictions)

    base_name = os.path.basename(tsv_path).replace(".tsv", "")
    #put in csv format
    
    print(f"\n=== {base_name} ===")
    print(f"Total predicted pairs: {len(matches) + len(not_in_ground_truth)}")
    print(f"Found in ground truth: {len(matches)}")
    print(f"Not in ground truth: {len(not_in_ground_truth)}")
    print(f"Predicted match mean: {mean_pred:.2f}")
    print(f"Predicted match median: {median_pred:.2f}")
    print(f"Predicted match min: {min_pred:.2f}")
    print(f"Predicted match max: {max_pred:.2f}")
    k=100
    print(f"Percentage of top {k} predicted pairs in top {k} GT pairs: {percent_top100(predicted_map, ground_truth, k):.2f}%")
    k=1000
    print(f"Percentage of top {k} predicted pairs in top {k} GT pairs: {percent_top100(predicted_map, ground_truth, k):.2f}%")
    k=2000
    print(f"Percentage of top {k} predicted pairs in top {k} GT pairs: {percent_top100(predicted_map, ground_truth, k):.2f}%")
    
    print(f"Percentage of predicted matches above GT mean: {above_mean(in_ground_truth, mean_score):.2f}%")


    

    # Save histogram and scores
    plt.hist(list_of_similar_predictions, bins=100)
    plt.xlabel('Predicted Alignment Score')
    plt.ylabel('Frequency')
    plt.title(f'{base_name} Score Distribution')
    plt.savefig(f"results/{base_name}_score_dist.png")
    plt.close()

    with open(f"results/{base_name}_matches.txt", "w") as out:
        for q, t, pred_len, true_len in matches:
            out.write(f"{q}\t{t}\t{pred_len}\t{true_len}\n")



csv_path = "finalorder.csv"
matches = []
not_in_ground_truth = []
predicted_map = {}
in_ground_truth = {}
list_of_similar_predictions = []
rowcnt=0
maxrows=2000
with open(csv_path) as f:
    reader = csv.reader(f)
    next(reader)  # Skip header
    for row in reader:
        rowcnt+=1
        if rowcnt>maxrows:
            break
        if len(row) >= 5:
            a, b = int(row[0]), int(row[1])
            score = float(row[4])
            predicted_map[(min(a, b), max(a, b))] = score
            pair = (min(a, b), max(a, b))
            if pair in ground_truth:
                true_len = ground_truth[pair]
                matches.append((a, b, score, true_len))
                list_of_similar_predictions.append(true_len)
                in_ground_truth[pair] = [true_len, score]
            else:
                not_in_ground_truth.append((a, b, score))
    
if not list_of_similar_predictions:
    print(f"{csv_path}: No overlaps found.")
    


mean_pred = sum(list_of_similar_predictions) / len(list_of_similar_predictions)
median_pred = sorted(list_of_similar_predictions)[len(list_of_similar_predictions) // 2]
min_pred = min(list_of_similar_predictions)
max_pred = max(list_of_similar_predictions)

base_name = "finalorder"
#put in csv format

print(f"\n=== {base_name} ===")
print(f"Total predicted pairs: {len(matches) + len(not_in_ground_truth)}")
print(f"Found in ground truth: {len(matches)}")
print(f"Not in ground truth: {len(not_in_ground_truth)}")
print(f"Predicted match mean: {mean_pred:.2f}")
print(f"Predicted match median: {median_pred:.2f}")
print(f"Predicted match min: {min_pred:.2f}")
print(f"Predicted match max: {max_pred:.2f}")
k=100
print(f"Percentage of top {k} predicted pairs in top {k} GT pairs: {percent_top100(predicted_map, ground_truth, k):.2f}%")
k=1000
print(f"Percentage of top {k} predicted pairs in top {k} GT pairs: {percent_top100(predicted_map, ground_truth, k):.2f}%")
k=2000
print(f"Percentage of top {k} predicted pairs in top {k} GT pairs: {percent_top100(predicted_map, ground_truth, k):.2f}%")

print(f"Percentage of predicted matches above GT mean: {above_mean(in_ground_truth, mean_score):.2f}%")




# Save histogram and scores
plt.hist(list_of_similar_predictions, bins=100)
plt.xlabel('Predicted Alignment Score')
plt.ylabel('Frequency')
plt.title(f'{base_name} Score Distribution')
plt.savefig(f"results/{base_name}_score_dist.png")
plt.close()

with open(f"results/{base_name}_matches.txt", "w") as out:
    for q, t, pred_len, true_len in matches:
        out.write(f"{q}\t{t}\t{pred_len}\t{true_len}\n")