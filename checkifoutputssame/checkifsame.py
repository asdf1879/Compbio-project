#read from aln_lh_500h.tsv and pairwise_results.txt and make sure that each line in pairwise_results.txt exists in aln_lh_500h.tsv
import pandas as pd
import os
import sys
import argparse
from os.path import join

lineslexic = []

filename1 = 'aln_lh_500h.tsv'
filename2 = 'pairwise_results.txt'

def read_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    return lines
def compare_files(lines1, lines2):

    for line in lines1:
        if line not in lines2:
            print(f"Line not found in second file: {line.strip()}")
            
        
    print("All lines in the first file are present in the second file.")

def main():
    lines1= read_file(filename1)
    lines2= read_file(filename2)
    compare_files(lines1, lines2)

if __name__ == "__main__":
    main()