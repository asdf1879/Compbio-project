# Compbio-project

Run using:

The `mycode` directory has the C++ code.  
Instructions to compile and run:

```bash
g++ -o2 -fopenmp lexic_hash.cpp -o lexichash.exe -O3 -march=native -funroll-loops
./lexichash.exe --out output --fasta NCTC1080_reads.fasta --n_hash 100 --alpha 1 --n_cpu 8
```

Output will be in `pairwiseresults.txt`.

---

Python code can be run using:

```bash
python run_module.py --out ../LH_out --fasta data/NCTC1080/NCTC1080_reads.fasta.gz --n_hash 500 --alpha 1
```

---

Mask generation and exponential decay aggregation are present in the `experiment_with_lexichash` directory.  
Comment out or change the relevant regions in `lexichash.py` in that directory to run with different masks, or uncomment the code before `write_overlaps`.  
Final results for exponential decay aggregation will be stored as CSV files.

---



Testing has been done on outputs and ground truth, present in various directories.  
Simply run the Python evaluation scripts after placing the results in the appropriate files.

---

Multiscale LexicHash is present in the 'Lexic_Hash1orig' directory and a modified version is present in 'LexicHash1'

Run using:

```bash
python run_module.py --out ../LH_out --fasta data/NCTC1080/NCTC1080_reads.fasta.gz --no_rc
```


---

## Quoting general instructions from the original LexicHash README

**LexicHash**  
*Sequence Similarity Estimation via Lexicographic Comparison of Hashes*

For help running LexicHash/MinHash, run the following in the LexicHash directory:

```bash
python3 run_module.py -h
```

Example on the NCTC1080 dataset:

```bash
python3 run_module.py --out ../LH_out --fasta data/NCTC1080/NCTC1080_reads.fasta.gz --n_hash 100 --max_k 32
```

To install the required dependencies:

```bash
pip install -r requirements.txt
```

> ⚠️ Note: LexicHash currently cannot process sequences containing `N`s.

---

### Ground Truth Alignment File Format:

```
[First read index (zero-based)]  
[Second read index (zero-based)]  
[Alignment size in base-pairs]  
[Second read alignment orientation (+/-)]
```

### Output Alignment File Format:

```
[First read index (zero-based)]  
[Second read index (zero-based)]  
[Alignment score]  
[Second read alignment orientation (+/-)]
```
