# Compbio-project


Run using 

mycode directory has the c++ code
Instructions to compile and run-

g++ -o2  -fopenmp lexic_hash.cpp -o lexichash.exe -O3 -march=native -funroll-loops 
lexichash.exe --out output --fasta NCTC1080_reads.fasta --n_hash 100 --alpha 1 --n_cpu 8  

Output in pairwiseresults.txt


Python code run using 

python run_module.py --out ../LH_out  --fasta  data/NCTC1080/NCTC1080_reads.fasta.gz --n_hash 500 --alpha 1


Mask generation and exponential decay aggregation present in the experiment_with_lexichash directory. Comment out/change the regions in the lexichash.py of the directory to run with different masks or run the code by uncommenting the code before write_overlaps. Final result for exponential decay aggregation will be stored in csv.


Quoting general instructions from the actual tools readme

LexicHash
Sequence Similarity Estimation via Lexicographic Comparison of Hashes
For help running LexicHash/MinHash, run the following in LexicHash directory:

python3 run_module.py -h

Example on NCTC1080 dataset:

python3 run_module.py --out ../LH_out --fasta data/NCTC1080/NCTC1080_reads.fasta.gz --n_hash 100 --max_k 32

To install requirements, run

pip install -r requirements.txt

For now, LexicHash cannot deal with sequences with N's.

The groundtruth alignmnent file has the following column format:

First read index (zero-based)
Second read index (zero-based)
Alignment size in base-pairs
Second read alignment orientation (+/-)
The output alignmnent file has the following column format:

First read index (zero-based)
Second read index (zero-based)
Alignment score
Second read alignment orientation (+/-)
