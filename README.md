# Compbio-project


Run using 

\test accuracy>python test.py
lexichash.exe --out output --fasta NCTC1080_reads.fasta --n_hash 100 --alpha 1 --n_cpu 15  
g++ -o2  -fopenmp lexic_hash.cpp -o lexichash.exe -O3 -march=native -funroll-loops 
