import re

# Load ground truth as a dictionary: (min, max) -> alignment size
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

#plot distribution of scores
import matplotlib.pyplot as plt
plt.hist(list_of_scores, bins=100)
plt.xlabel('Alignment Score')
plt.ylabel('Frequency')
plt.title('Distribution of Alignment Scores')
plt.savefig('alignment_scores_distribution.png')
plt.show()


sorted_scores = sorted(list_of_scores)
mean_score = sum(sorted_scores) / len(sorted_scores)
print(f"Mean alignment score: {mean_score:.2f}")
#//print min and max scores
min_score = sorted_scores[0]
max_score = sorted_scores[-1]
print(f"Min alignment score: {min_score:.2f}")
print(f"Max alignment score: {max_score:.2f}")
#//print median score
median_score = sorted_scores[len(sorted_scores) // 2]
print(f"Median alignment score: {median_score:.2f}")
#//print top 100
print("Top 100 scores:")
for i in range(100):
    print(sorted_scores[-i-1])


# Load pairwise results and check matches
matches = []
not_in_ground_truth = []

list_of_similar_predictions = []

with open("pairwiseresultspython.txt", "r") as f:
    for line in f:
        parts = re.split(r"\s+", line.strip())
        if len(parts) >= 4:
            q, t, s, l = parts[:4]
            q, t = int(q) , int(t)   # convert 1-based to 0-based
            pair = (min(q, t), max(q, t))
            predicted_len = int(l)

            if pair in ground_truth:
                true_len = ground_truth[pair]
                matches.append((q, t, predicted_len, true_len))
                list_of_similar_predictions.append(true_len)
            else:
                not_in_ground_truth.append((q, t, predicted_len))

meanofpredicted = sum(list_of_similar_predictions) / len(list_of_similar_predictions)
print(f"Mean predicted alignment score: {meanofpredicted:.2f}")
#//print min and max predicted scores
min_predicted = min(list_of_similar_predictions)
max_predicted = max(list_of_similar_predictions)
print(f"Min predicted alignment score: {min_predicted:.2f}")
print(f"Max predicted alignment score: {max_predicted:.2f}")
#//print median predicted score
median_predicted = sorted(list_of_similar_predictions)[len(list_of_similar_predictions) // 2]
print(f"Median predicted alignment score: {median_predicted:.2f}")

#plot distribution of predicted scores
plt.hist(list_of_similar_predictions, bins=100)
plt.xlabel('Predicted Alignment Score')
plt.ylabel('Frequency')
plt.title('Distribution of Predicted Alignment Scores')
plt.savefig('predicted_alignment_scores_distribution.png')
plt.show()



# Output results
print(f"Total predicted pairs: {len(matches) + len(not_in_ground_truth)}")
print(f"Found in ground truth: {len(matches)}")
print(f"Not in ground truth: {len(not_in_ground_truth)}\n")

print("Matches with true and predicted lengths:")
for q, t, pred_len, true_len in matches[:10]:  # show first 10
    print(f"{q}\t{t}\tPredicted: {pred_len}\tTrue: {true_len}")

# Optional: save to file
with open("matches_with_scores.txt", "w") as out:
    for q, t, pred_len, true_len in matches:
        out.write(f"{q}\t{t}\t{pred_len}\t{true_len}\n")
