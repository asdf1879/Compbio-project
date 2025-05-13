import re
import pandas as pd

# Load ground truth as a dictionary: (min, max) -> alignment size

#percentage of top 100 predicted pairs in top 100 groundtruth score pairs

def percent_top100(predicted, ground_truth):
    # Sort the ground truth and predicted pairs based on alignment scores in descending order
    
    sorted_ground_truth = sorted(ground_truth.items(), key=lambda x: x[1], reverse=True)
    
    #get the i1,i2 -> score pairs
    ground_truth_pairs = {pair: score for pair, score in sorted_ground_truth[:2000]}

    percent_overlap=0

    keys = list(predicted.keys())
    # print(keys[1])

    for key in predicted:
        # print(key)
        if key in ground_truth_pairs:
            percent_overlap += 1

    
    return percent_overlap/20

def topkexp(df,ground_truth,k):
    # Sort the ground truth and predicted pairs based on alignment scores in descending order
    sorted_ground_truth = sorted(ground_truth.items(), key=lambda x: x[1], reverse=True)
    
    #get the i1,i2 -> score pairs
    ground_truth_pairs = {pair: score for pair, score in sorted_ground_truth[:k]}

    percent_overlap=0
    cnt=0
    for row in df.iterrows():
        # print(row)
        # print(row[1])
        i1 = row[1]['id1']
        i2 = row[1]['id2']
        # print(i1,i2)
        cnt += 1
        if (i1,i2) in ground_truth_pairs:
            percent_overlap += 1
        if cnt == k:
            break
    # print(keys[1])

    

    print(f"number of pairs in top k: {cnt} is {percent_overlap/k*100:.2f}%")
    return percent_overlap

def above_mean(predicted_map,mean,len):
    #
    #get the i1,i2 -> score pairs
   
    cnt=0


    keys = list(predicted_map.keys())
    # print(keys[1])

    for key in predicted_map:
        # print(key)
        [x,y]= predicted_map[key]
        if x> mean:
            cnt += 1
    print(mean)
    print(cnt)
    
    return (cnt/len)*100


import pandas as pd
import matplotlib.pyplot as plt
import pyfastx

def scatterplot(groundtruthwal):
    # Convert dictionary to DataFrame and transpose
    df = pd.DataFrame(groundtruthwal).T
    df.columns = ['True Length', 'Predicted Length']
    plt.figure(figsize=(8, 6))
    # Predicted on x, True on y
    plt.scatter(df['Predicted Length'], df['True Length'], alpha=0.6)
    plt.xlabel('Predicted Length')
    plt.ylabel('True Length')
    plt.title('True vs Predicted Lengths')
    # Diagonal reference line (limited by x-axis)
    plt.plot([14, 32], [0, 20000], 'r--', label='Reference')
    plt.xlim(14, 33)
    plt.ylim(0, 20000)
    plt.legend()
    plt.tight_layout()
    plt.savefig('scatterplot.png')
    plt.close()



seqs = pyfastx.Fastx("NCTC1080_reads.fasta.gz")
seqs = list(pyfastx.Fastx("NCTC1080_reads.fasta.gz"))

seqs = [seq[1] for seq in seqs]



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
# plt.show()


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
# print("Top 100 scores:")
# for i in range(100):
#     print(sorted_scores[-i-1])


# Load pairwise results and check matches
matches = []
not_in_ground_truth = []

predicted_map = {}

in_ground_truth ={}

list_of_similar_predictions = []

with open("pairwise_results.txt", "r") as f:
    for line in f:
        parts = re.split(r"\s+", line.strip())
        if len(parts) >= 4:
            q, t, l, s = parts[:4]
            q, t = int(q)-1 , int(t)-1   # convert 1-based to 0-based
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
# plt.show()



# Output results
print(f"Total predicted pairs: {len(matches) + len(not_in_ground_truth)}")
print(f"Found in ground truth: {len(matches)}")
print(f"Not in ground truth: {len(not_in_ground_truth)}\n")

# print("Matches with true and predicted lengths:")
# for q, t, pred_len, true_len in matches[:10]:  # show first 10
#     print(f"{q}\t{t}\tPredicted: {pred_len}\tTrue: {true_len}")

# Optional: save to file
with open("matches_with_scores.txt", "w") as out:
    for q, t, pred_len, true_len in matches:
        out.write(f"{q}\t{t}\t{pred_len}\t{true_len}\n")


percent_overlap = percent_top100(predicted_map, ground_truth)
print(f"Percentage of top 2000 predicted pairs in top 2000 ground truth pairs: {percent_overlap:.2f}%")



print(f"Percentage of predicted pairs above mean: {above_mean(in_ground_truth,mean_score,len(in_ground_truth)):.2f}%")

scatterplot(in_ground_truth)



df =pd.read_csv('finalorder.csv')


k=1000
topkindf = df.head(k)

print(topkindf.head())
(topkexp(topkindf,ground_truth,k))

map_exp ={}

for i in range(len(df)):
    row = df.iloc[i]
    i1 = row['id1']
    i2 = row['id2']
    pred_len = row['score']
    map_exp[(i1,i2)] = pred_len




from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from scipy.stats import spearmanr, pearsonr

# Assume these are your dictionaries
# predicted_scores = { (i1, i2): val1, ... }
# ground_truth = { (i1, i2): val2, ... }

# Find intersection keys (if needed)
common_keys = set(map_exp.keys()) & set(ground_truth.keys())

# Extract aligned values
pred_vals = [map_exp[k] for k in common_keys]
gt_vals = [ground_truth[k] for k in common_keys]

# Optional: Normalize to same scale
from sklearn.preprocessing import MinMaxScaler, StandardScaler, RobustScaler, MaxAbsScaler
import numpy as np

scaler = MinMaxScaler()
pred_scaled = scaler.fit_transform(np.array(pred_vals).reshape(-1, 1)).flatten()
gt_scaled = scaler.fit_transform(np.array(gt_vals).reshape(-1, 1)).flatten()

# Compute metrics
mse = mean_squared_error(gt_scaled, pred_scaled)
mae = mean_absolute_error(gt_scaled, pred_scaled)
r2 = r2_score(gt_scaled, pred_scaled)
spearman_corr, _ = spearmanr(gt_scaled, pred_scaled)
pearson_corr, _ = pearsonr(gt_scaled, pred_scaled)

# Print results
print(f"MSE: {mse:.4f}")
print(f"MAE: {mae:.4f}")
print(f"R² Score: {r2:.4f}")
print(f"Spearman Correlation: {spearman_corr:.4f}")
print(f"Pearson Correlation: {pearson_corr:.4f}")