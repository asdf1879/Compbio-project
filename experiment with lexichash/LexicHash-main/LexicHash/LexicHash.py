import numpy as np
import sys
from os.path import exists
from multiprocessing import Pool, cpu_count
from time import perf_counter, process_time
import pyfastx
from LexicHash.utils import utils, seq_utils

BASE_TO_INT = {'A':0, 'T':1, 'G':2, 'C':3, 'a':0, 't':1, 'g':2, 'c':3}

def find_overlaps(**args):
    """
    Runs the LexicHash pipeline on a set of long-read sequencing reads
    to find similarities between pairs of reads.

    Args (in args dictionary):
        fasta_path: Path of the reads in fasta format
        out: Path of the output directory to write alignment file
        n_hash: Number of hash functions to consider (int)
        min_k: Minimum match length to consider for output (int)
        max_k: Maximum match length to consider for output 
        rc: Whether or not to consider reverse-complement reads (boolean)
        alpha: Outputs alpha*n_seq total pairs
        
    Globals:
        sketches: Sketches of reads. Each row corresponds to a single read. 
            Each column corresponds to a single mask (i.e. hash function)
        masks: Each mask indicates lexicographic order for sketching procedure.
            `masks` has shape (n_hash, MAX_K)
        MIN_K: Minimum match-length to consider for output
        MAX_K: Maximum match-length to consider for output
        RC: Whether to consider reverse-complement reads in whole pipeline (typically True)
    """
    global MAX_K 
    MAX_K= args['max_k']
    start_program = perf_counter()
    
    # OBTAIN READS (no RAM)
    seqs = pyfastx.Fastx(args['fasta_path'])
    
    # SKETCHING
    utils.print_banner('SKETCHING READS')
    seqs = list(pyfastx.Fastx(args['fasta_path']))  # Force full read into list of (id, seq)

    masks = get_masks(args['n_hash'], args['max_k'])
    start = perf_counter()
    sketches, n_seq = sketching(seqs, masks, **args)
    print(f'# sequences: {n_seq}')
    utils.print_clocktime(start, perf_counter(), 'sketching')

    # PAIRWISE COMPARISON
    utils.print_banner('PERFORMING PAIRWISE COMPARISON')
    start = perf_counter()
    pair_match_lens = pairwise_comparison(sketches, n_seq, **args)

    utils.print_clocktime(start, perf_counter(), 'pairwise comparison')

    # # # REFINE MATCH PAIRS
    sorted_df = refine_match_pairs(pair_match_lens,sketches)
    sorted_df.to_csv('finalorder.csv', index=False)
    
    
    # WRITE RESULTS
    write_overlaps(pair_match_lens, **args)
    utils.print_clocktime(start_program, perf_counter(), 'full process')
    
    
def write_overlaps(pair_match_lens, aln_path, **args):
    '''
    Writes the overlap pairs to output file. 
    
    Args:
        pair_match_lens: Dict of pair:match-length. Pair is a tuple of the form (id1,id2,+/-).
        aln_path: Output tsv filepath. Each row is a pair with columns (id1,id2,+/-,match-length)
    '''
    print(f'# overlaps: {len(pair_match_lens)}')
    with open(aln_path, 'w') as fh:
        for pair,match_len in pair_match_lens.items():
            id1,id2,is_rc = pair
            line = '\t'.join([str(id1+1), str(id2+1), str(match_len), is_rc]) # 1 indexed
            fh.write(line+'\n')
                
#################################################
#                SKETCHING                      #
#################################################

    
def sketching(seqs, masks, n_hash, n_cpu, rc, max_k, min_k, **args):
    '''
    Use multiprocessing to compute the sketches for all sequences. 
    
    Returns:
        sketches: Sketches of reads. Each row corresponds to a single read. 
            Each column corresponds to a single mask (i.e. hash function).
            Includes reverse-complements unless no_rc is set
        n_seq: Number of sequences
    '''
    chunksize = 10
    with Pool(processes=n_cpu, initializer=init_worker_sketching, initargs=(masks,rc,max_k,min_k)) as pool:
        all_sketches = pool.map(get_seq_sketch, seqs, chunksize)
    sketches, sketches_rc = list(map(list, zip(*all_sketches))) # list of two-tuples to two lists
    n_seq = len(sketches)
    if max_k <= 32:
        sketches = np.array(sketches, dtype=np.uint64)
        if rc: sketches = np.concatenate((sketches, np.array(sketches_rc,dtype=np.uint64)))
    else:
        sketches = np.array(sketches, dtype='object')
        if rc: sketches = np.concatenate((sketches, np.array(sketches_rc, dtype='object')))
    return sketches, n_seq


def init_worker_sketching(masks, rc, max_k, min_k):
    global shared_masks, RC, MAX_K, MIN_K
    shared_masks = masks
    MAX_K = max_k
    MIN_K = min_k
    RC = rc
    
    
def get_seq_sketch(seq): 
    '''
    Function called by multiprocessing.
    Computes the sketch of a single sequence.
    
    Args:
        seq: Sequence to sketch
        
    Returns:
        sketch: Sketch of the read. Each entry corresponds to a single mask (i.e. hash function)
        sketch_rc: Sketch of reverse-complement read
    '''
    seq = seq[1] # pyfastx.Fastx iters as (name,seq,comment)
    seq_int_repr = np.array([BASE_TO_INT[b] for b in seq])
    sketch = [lexicographic_first(seq_int_repr, mask) for mask in shared_masks]
    if RC:
        seq_rc = seq_utils.revcomp(seq)
        seq_int_repr = np.array([BASE_TO_INT[b] for b in seq_rc])
        sketch_rc = [lexicographic_first(seq_int_repr, mask) for mask in shared_masks]
    else:
        sketch_rc = None
    return sketch, sketch_rc
    
    
def lexicographic_first(seq,mask):
    '''
    Finds the lexicographic first substring of input sequence given 
        lexicographic order indicated by the input mask. 
    First finds candidate locations in the sequence which correspond to
        the starting locations of the longest substrings exactly matching 
        the prefix of the mask.
    Next extends these candidate substrings to find the lexicograhpic first
        one according the mask.
    
    Args:
        seq: Sequence to consider
        mask: Lexicographic order to consider
        
    Returns:
        min-hash: Rank (hash value) of lexicographic first substring.
    '''
    def get_candidate_locs():
        j = 0 
        cur_idxs = np.arange(len(seq)-MIN_K)
        while j < MAX_K:
            next_idxs = np.where(seq[cur_idxs] == mask[j])[0]
            if len(next_idxs)==0:
                return cur_idxs-j, j
            elif len(next_idxs)==1:
                return [cur_idxs[next_idxs[0]]-j], j+1
            cur_idxs = cur_idxs[next_idxs]+1
            if cur_idxs[-1] == len(seq): 
                cur_idxs = cur_idxs[:-1]  
            j += 1
        return cur_idxs-MAX_K, MAX_K
    
    def extend_candidates(candidate_locs, n_matching):
        j = n_matching
        while len(candidate_locs) > 1 and j < MAX_K:
            best = 4 # will always be overwritten immediately
            next_candidates = []
            for loc in candidate_locs:
                if loc+j < len(seq):
                    bm = (seq[loc+j]-mask[j]) % 4
                    if bm < best:
                        best = bm
                        next_candidates = [loc]
                    elif bm == best:
                        next_candidates.append(loc)
            if len(next_candidates) == 0: break
            candidate_locs = next_candidates
            j += 1    
        return candidate_locs[0]
    
    def hash_val(lex_first_idx):
        val = 0
        for b in seq[lex_first_idx:lex_first_idx+MAX_K]:
            val *= 4 # 2 bits per base, and it's 2 bases
            val += b.item()
        base_extend = max(0, MAX_K-(len(seq)-lex_first_idx))
        if base_extend > 0: # used if substring is at end of sequence (corner case)
            for b in mask[-base_extend:]:
                val *= 4 # 2 bits per base
                val += b.item() ^ 3 # extend with lexicographically last rank
        return int(val)

    candidate_locs, n_matching = get_candidate_locs() # n_matching is length exactly matching the mask
    lex_first_idx = extend_candidates(candidate_locs, n_matching)
    return hash_val(lex_first_idx)


def get_masks(n_masks, max_k): 
    '''
    Create randomized masks.
    
    Args:
        n_masks: Number of masks to create. Corresponds to input argument n_hash.
    Returns:
        masks: Numpy array (n_masks, max match length). Entries are integers in {0,1,2,3}.
    ''' 
    masks = np.random.randint(4, size=(n_masks,max_k))
    print(masks)
    return masks



def get_masks_uniform_blend(n_masks, max_k, sequences):
    '''
    Mix of random and data-derived masks.
    sequences: list of tuples (id, sequence_str)
    '''
    import numpy as np

   
    
    random_masks = np.random.randint(4, size=(n_masks // 2, max_k))

    # Convert sequences to just the strings
    seq_strings = [s[1] for s in sequences]

    data_masks = []
    while len(data_masks) < (n_masks - len(random_masks)):
        seq = np.random.choice(seq_strings)
        if len(seq) < max_k:
            continue
        start = np.random.randint(0, len(seq) - max_k + 1)
        mask_str = seq[start:start+max_k]
        mask_ints = [BASE_TO_INT[b] for b in mask_str]
        data_masks.append(mask_ints)

    data_masks = np.array(data_masks)
    masks = np.vstack([random_masks, data_masks])
    print(len(masks))
    return masks

def get_masks_freq_weighted(n_masks, max_k, sequences):
    '''
    Sample rare k-mers more often for mask creation.
    sequences: list of tuples (id, sequence_str)
    '''
    import numpy as np
    from collections import Counter

    # Convert sequences to just the strings
    seq_strings = [s[1] for s in sequences]

    kmers = []
    for seq in seq_strings:
        if len(seq) < max_k:
            continue
        # Extract all k-mers
        for i in range(len(seq) - max_k + 1):
            kmers.append(seq[i:i+max_k])

    # Count k-mer frequencies
    kmer_counts = Counter(kmers)

    # Assign inverse frequency weights
    weights = np.array([1 / (kmer_counts[k] + 1) for k in kmers])
    probs = weights / weights.sum()

    # Sample k-mers according to probabilities
    sampled = np.random.choice(kmers, size=n_masks, p=probs)

    masks = np.array([[BASE_TO_INT[c] for c in k] for k in sampled])
    print(masks)
    return masks
def get_masks_cluster_gc(n_masks, max_k, sequences):
    '''
    Generate masks using GC-content based clustering.
    sequences: list of tuples (id, sequence_str)
    '''
    import numpy as np

    

    low_gc, high_gc = [], []

    # Split sequences into low-GC and high-GC
    for _, seq in sequences:
        if not seq:  # skip empty
            continue
        gc_content = (seq.count('G') + seq.count('C')) / len(seq)
        if gc_content > 0.5:
            high_gc.append(seq)
        else:
            low_gc.append(seq)

    masks = []
    while len(masks) < n_masks:
        # Choose low_gc or high_gc list randomly
        pool = low_gc if np.random.rand() < 0.5 else high_gc
        if not pool:  # if the pool is empty
            continue
        seq = np.random.choice(pool)
        if len(seq) < max_k:
            continue
        start = np.random.randint(0, len(seq) - max_k + 1)
        mask_str = seq[start:start + max_k]
        mask = [BASE_TO_INT[b] for b in mask_str]
        masks.append(mask)

    return np.array(masks)


def get_masks_mutated_data(n_masks, max_k, sequences, mutation_rate=0.1):
    '''
    Sample real substrings and add random mutations.
    sequences: list of tuples (id, sequence_str)
    '''
    import numpy as np

    # Convert sequences to just the strings
    seq_strings = [s[1] for s in sequences]

    masks = []
    while len(masks) < n_masks:
        seq = np.random.choice(seq_strings)
        if len(seq) < max_k:
            continue
        start = np.random.randint(0, len(seq) - max_k + 1)
        kmer = list(seq[start:start+max_k])

        # Mutate each base with probability 'mutation_rate'
        for i in range(max_k):
            if np.random.rand() < mutation_rate:
                original = kmer[i]
                options = [b for b in 'ACGT' if b != original]
                kmer[i] = np.random.choice(options)

        mask_ints = [BASE_TO_INT[b] for b in kmer]
        masks.append(mask_ints)

    masks = np.array(masks)
    print(masks)
    return masks

def get_masks_progressive_refresh(n_masks, max_k, sequences, buffer_size=300):
    '''
    Start with random masks, include a few sampled from a rolling buffer of seen sequences.
    sequences: list of tuples (id, sequence_str)
    '''
    import numpy as np

    # Convert sequences to just the strings
    seq_strings = [s[1] for s in sequences]

    # Build buffer from first few sequences
    buffer_kmers = []
    for seq in seq_strings[:buffer_size]:
        if len(seq) < max_k:
            continue
        start = np.random.randint(0, len(seq) - max_k + 1)
        buffer_kmers.append(seq[start:start+max_k])

    n_buffer = int(n_masks * 0.3)
    n_random = n_masks - n_buffer

    random_masks = np.random.randint(4, size=(n_random, max_k))
    buffer_masks = np.array([
        [BASE_TO_INT[b] for b in kmer]
        for kmer in np.random.choice(buffer_kmers, size=n_buffer)
    ])

    masks = np.vstack([random_masks, buffer_masks])
    print(masks)
    return masks

#################################################
#            PAIRWISE COMPARISON                #
#################################################

def pairwise_comparison(sketches, n_seq, n_hash, n_cpu, alpha, max_k, min_k, **args): 
    """
    Perform the pairwise comparison component of the LexicHash pipeline.

    Args:
        n_seq: Number of sequences
        n_hash: Number of hash functions used.
        n_cpu: Number of cpus to use
        alpha: Will aggregate alpha*n_seq pairs
        
    Returns:
        pair_match_lens: Dict of pair:match-length. Pair is a tuple of the form (id1,id2,+/-).
    """
    n_pairs = n_seq * alpha if alpha is not None else 10000000000000000000000
    print(f'number of pairs to output: {n_pairs}')
    all_matching_sets = prefix_tree_multiproc(sketches, n_hash, n_cpu, max_k, min_k)
#     min_k = get_k_thresh(all_matching_sets) if min_k is None else min_k
    pair_match_lens = bottom_up_pair_aggregation(all_matching_sets, n_seq, n_pairs, max_k, min_k)

    
    return pair_match_lens

    
def prefix_tree_multiproc(sketches, n_hash, n_cpu, max_k, min_k):
    '''
    Use multiprocessing to compute the pairwise comparisons. 
    Each process considers a single hash function.
    
    Returns:
        sketches: Sketches of reads. Each row corresponds to a single read. 
            Each column corresponds to a single mask (i.e. hash function)
        sketches_rc: Sketches of reverse-complement reads
    '''
    args = (i for i in range(n_hash))
    chunksize = int(np.ceil(n_hash/n_cpu/4))
    with Pool(processes=n_cpu, initializer=init_worker_prefix_tree, initargs=(sketches,max_k,min_k)) as pool:
        all_matching_sets = pool.map(get_matching_sets, args, chunksize)
    return all_matching_sets


def init_worker_prefix_tree(sketches, max_k, min_k):
    global shared_sketches, MAX_K, MIN_K, RC, get_chars
    shared_sketches = sketches
    MAX_K = max_k
    MIN_K = min_k
    if max_k <= 32:
        get_chars = lambda x,k: (x >> 2*(MAX_K-k-1)) & 3
    else:
        and3_ = lambda y: int(y)&3
        get_chars = lambda x,k: map(and3_, x // 4**(MAX_K-k-1))

    
def get_matching_sets(sketch_idx):
    '''
    Function called by multiprocessing.
    Efficiently computes the pairwise comparison using a structure similar to a prefix-tree.
    Starts at root of the prefix tree, which is represented as the set of all sequence indices.
    Next, partition the root set based on the groups of substrings with matching first base.
    Continue partitioning until a depth of MAX_K or until no non-singleton sets remain.
    The partition at depth k is a list of sets, where the number of sets corresponds to the number of
        subtrees at depth k in the prefix tree, and each set corresponds to the sequence indices 
        in that subtree.
    
    Args:
        sketch_idx: Index in sketch to consider (i.e. hash function index)
        
    Returns:
        matching_sets: Dict of match-length: list of sets of sequence indices
            with match of corresponding length 
    '''
    subtrees = [list(range(len(shared_sketches)))]
    matching_sets = {}
    for k in range(MAX_K): # current position in all substrings k (i.e. match-length - 1)
        next_subtrees = []
        for seq_idxs in subtrees: 
            partition = {0:[],1:[],2:[],3:[]}
#             print(shared_sketches[seq_idxs, sketch_idx], 2*(MAX_K-k-1))
            chars = get_chars(shared_sketches[seq_idxs, sketch_idx], k)
            for char,seq_idx in zip(chars,seq_idxs):
                partition[char].append(seq_idx)
            partition = [p for p in partition.values() if len(p)>1]
            next_subtrees.extend(partition)
        if len(next_subtrees) == 0:
            return matching_sets
        subtrees = next_subtrees
        if k+1 >= MIN_K: # k+1 because a match at index 0 is a length-1 match, e.g.
            matching_sets[k+1] = subtrees.copy()
    return matching_sets
                   

def bottom_up_pair_aggregation(all_matching_sets, n_seq, n_pairs, max_k, min_k):
    '''
    Aggregate pairs of similar sequences from the ``bottom" up. Starts with maximum match-length (MAX_K),
        and for each set in the list of sets at that level, it adds all possible pairs of sequence indices
        to the dictionary.
        
    Args:
        matching_sets_comb: Single giant dict of match-length:list of all sets of sequence indices
    
    Returns:
        pair_match_lens: Dict of pairs:match-length. Pair is a tuple of the form (id1,id2,+/-).
    '''
    # print(f'number of pairs to output: {n_pairs}')
    # print(f'n_seq: {n_seq}')
# go from combined dict to dict of seq idx pairs (i1,i2) --> (match_len, +/-)
    pair_match_lens = {}
    for k in range(max_k, min_k-1, -1): # start from bottom to give max match length to each pair
        for matching_sets in all_matching_sets:
            for matching_set in matching_sets.get(k, []):
                for i1 in matching_set: # seq index 1
                    for i2 in matching_set: # seq index 2
                        _i1,_i2 = i1%n_seq, i2%n_seq
                        if _i2 > _i1:
                            is_rc = '-' if (i1>=n_seq and i2<n_seq) or (i1<n_seq and i2>=n_seq) else '+'
                            key = (_i1,_i2,is_rc)
                            if key not in pair_match_lens:
                                pair_match_lens[key] = k
                                
                                if len(pair_match_lens) == int(n_pairs):
                                    return pair_match_lens
    return pair_match_lens



def prefix_match_len_bases(sketch1, sketch2):
    """
    Calculate the prefix match length in bases between two sequences.
    
    Args:
        s1: First sequence
        s2: Second sequence
    
    Returns:
        match_len: The prefix match length in bases
    """
    matchlen=0

    bitrepr1 = []
    bitrepr2 = []

    for i in range(MAX_K):
        b1 = (sketch1 >> 2*(MAX_K-i-1)) & 3
        b2 = (sketch2 >> 2*(MAX_K-i-1)) & 3
        if b1 == b2:
            matchlen += 1
        else:
            break
        # bitrepr1.append(b1)
        # bitrepr2.append(b2)
        
    return matchlen
    # print(f"sketch1: {sketch1}, sketch2: {sketch2}, matchlen: {matchlen}")
    # print(bitrepr1)
    # print(bitrepr2)

def yscore(ind1,ind2,sketches):
    """
    Calculate the y-score for a pair of indices.
    
    Args:
        ind1: First index
        ind2: Second index
    
    Returns:
        y_score: The y-score for the pair of indices
    """
    #get all sketches of ind1 and ind2
    shared_sketches=sketches    
    sketch1 = shared_sketches[ind1]
    sketch2 = shared_sketches[ind2]

    #get the y score which is list of prefix match lengths
    y_score = []

    for i in range(len(sketch1)):
        #get the prefix match length
        # print(f"sketch1: {sketch1[i]}, sketch2: {sketch2[i]}")
        y_score.append(prefix_match_len_bases(sketch1[i], sketch2[i]))

    #return the y score
    return y_score

def expdecayscore(y_score, decay):

    revsorted_y_score = sorted(y_score, reverse=True)
    # print(revsorted_y_score)
    n = len(revsorted_y_score)
    
    finalscore=0
    for i in range(max(n,100)):
        finalscore += revsorted_y_score[i] * (1-decay)**i

    return finalscore
import pandas as pd

def refine_match_pairs(pair_match_lens, sketches):
    '''
    Refine the match pairs to only include the top n_pairs.

    Args:
        pair_match_lens: Dict of pairs:match-length. Pair is a tuple of the form (id1,id2,+/-).
        sketches: Data used to compute scores

    Returns:
        sorted_df: DataFrame sorted by score
    '''
    final_key_val = {}
    records = []

    for key in pair_match_lens:
        ind1, ind2, is_rc = key
        val = pair_match_lens[key]
        y_score = yscore(ind1, ind2, sketches)
        score = expdecayscore(y_score, 0.1)
        final_key_val[key] = score
        records.append({'id1':ind1, 'id2':ind2, 'match_len': val, 'is_rc': is_rc, 'score': score})

    df = pd.DataFrame.from_records(records)
    sorted_df = df.sort_values(by='score', ascending=False).reset_index(drop=True)
    
    return sorted_df

# TO USE, YOU MUST CHANGE bottom_up_pair_aggregation SO THAT YOU CAN SET THE MIN_K, THEN GET 1 FULL SAMPLE TREE, THEN CHANGE THE get_k_thresh ACCORDINGLY 
def get_k_thresh(all_matching_sets, max_quantile=0.99):
    '''
    This function finds a minimum match-length to consider based on the output of the prefix-tree step.
    Calculates the match-length PMF by determining the number of pairs at each length. 
    The minimum match-length is set to the inflection point of the log-PMF. 
    This method has proven to be reliable.
    '''
    def get_all_n_pairs_arr():
        k_vals = np.arange(MAX_K+1)
        all_n_pairs_arr = np.zeros(MAX_K+1)
        for matching_sets in all_matching_sets:
            cur_n_pairs_arr = get_n_pairs_arr(matching_sets)
            for k,n_pairs in zip(k_vals,cur_n_pairs_arr):
                all_n_pairs_arr[k] += n_pairs
        return k_vals, all_n_pairs_arr
    
    def get_n_pairs_arr(matching_sets):
        n_pairs_arr = np.zeros(MAX_K+1)
        for k in np.sort(list(matching_sets))[::-1]:
            n_pairs = 0
            for s in matching_sets[k]:
                n_pairs += len(s)*(len(s)-1)//2
            n_pairs_arr[k] = n_pairs
            if k <= MAX_K: 
                n_pairs_arr[k] -= sum(n_pairs_arr[k+1:])
        return n_pairs_arr
  
    def get_extreme_value_cdf(n_pairs_arr):
        n_hash = len(sketches[0])
        cdf = np.cumsum(n_pairs_arr)/sum(n_pairs_arr)
        ev_cdf = cdf**n_hash
        return ev_cdf
    
    def get_maximum_k_thresh(cdf, max_quantile=0.9):
        idx = np.flatnonzero(np.diff(cdf>max_quantile))[0]
        return k_vals[idx]

    def get_inflection_thresh(k_vals, cdf):
        log_cdf = np.log(np.gradient(cdf)+1)
        log_d2 = np.gradient(np.gradient(log_cdf))
        smooth_log_d2 = gaussian_filter1d(log_d2, 1)
        smooth_log_d2 = np.insert(smooth_log_d2, 0, 1000)
        infls = np.flatnonzero(np.diff(np.sign(smooth_log_d2)))
        k_thresh = k_vals[infls[1]]+1
        return k_thresh
    
    k_vals, n_pairs_arr = get_all_n_pairs_arr()
    cdf = get_extreme_value_cdf(n_pairs_arr)
    k_thresh = get_inflection_thresh(k_vals, cdf)
    maximum_k_thresh = get_maximum_k_thresh(cdf)
    k_thresh = min(k_thresh, maximum_k_thresh)
    return k_thresh

