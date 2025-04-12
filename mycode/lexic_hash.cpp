#include <bits/stdc++.h>
#include <omp.h> // OpenMP
// #include <seqan3/io/sequence_file/input.hpp>

using namespace std;
int min_k = 14;
int max_k = 32;

int n_hash = 500;
int n_cpu = 4;
int MIN_K = 14;
int MAX_K = 32;
bool RC = false;
vector<vector<int>> shared_masks;
unordered_map<char, int> BASE_TO_INT = {
    {'A', 0}, {'T', 1}, {'G', 2}, {'C', 3}, {'a', 0}, {'t', 1}, {'g', 2}, {'c', 3}};

vector<string> read_fasta(const string &fasta_path, size_t max_k) {
    vector<string> seqs;
    ifstream infile(fasta_path);
    if (!infile) {
        cerr << "Error opening file: " << fasta_path << "\n";
        return seqs;
    }

    string line, current_seq;
    while (getline(infile, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!current_seq.empty()) {
                if (current_seq.size() >= max_k)
                    seqs.push_back(current_seq);
                current_seq.clear();
            }
        } else {
            current_seq += line;
        }
    }

    // Add last sequence
    if (!current_seq.empty() && current_seq.size() >= max_k)
        seqs.push_back(current_seq);

    return seqs;
}


vector<vector<int>> get_masks(int n_masks, int max_k)
{
    vector<vector<int>> masks(n_masks, vector<int>(max_k));
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(0, 3);

    for (int i = 0; i < n_masks; ++i)
    {
        for (int j = 0; j < max_k; ++j)
        {
            masks[i][j] = dist(gen);
        }
    }
    return masks;
}
string reverse_complement(const string &seq)
{
    string rc;
    rc.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it)
    {
        switch (*it)
        {
        case 'A':
        case 'a':
            rc += 'T';
            break;
        case 'T':
        case 't':
            rc += 'A';
            break;
        case 'G':
        case 'g':
            rc += 'C';
            break;
        case 'C':
        case 'c':
            rc += 'G';
            break;
        default:
            rc += 'N'; // unknown
        }
    }
    return rc;
}
vector<int> to_int_seq(const string &seq)
{
    vector<int> int_seq(seq.size());
    for (size_t i = 0; i < seq.size(); ++i)
        int_seq[i] = BASE_TO_INT[seq[i]];
    return int_seq;
}

pair<vector<int>, int> get_candidate_locs(const vector<int> &seq, const vector<int> &mask)
{
    int j = 0;
    vector<int> cur_idxs;
    for (int i = 0; i < (int)(seq.size()) - MIN_K; ++i)
        cur_idxs.push_back(i);
    // cout<<"curidxs ";
    // for(int i=0;i<cur_idxs.size();i++){
    //     cout<<cur_idxs[i]<<" ";
    // }
    // cout<<endl;

    while (j < MAX_K)
    {
        vector<int> next_idxs;
        for (int idx = 0; idx < cur_idxs.size(); ++idx)
        {
            if (cur_idxs[idx] < seq.size() && seq[cur_idxs[idx]] == mask[j])
                next_idxs.push_back(idx);
        }
        // cout << "nextids ";
        // for (int i = 0; i < next_idxs.size(); i++)
        // {
        //     cout << next_idxs[i] << " ";
        // }
        // cout << endl;

        if (next_idxs.empty())
        {
            vector<int> result;
            for (int idx : cur_idxs)
                result.push_back(idx - j);
            return {result, j};
        }
        else if (next_idxs.size() == 1)
        {
            return {{cur_idxs[next_idxs[0]] - j}, j + 1};
        }

        vector<int> new_cur_idxs;
        for (int idx : next_idxs)
        {
            if (cur_idxs[idx] + 1 < seq.size())
                new_cur_idxs.push_back(cur_idxs[idx] + 1);
        }

        cur_idxs = new_cur_idxs;
        ++j;
    }

    vector<int> result;
    for (int idx : cur_idxs)
        result.push_back(idx - MAX_K);
    return {result, MAX_K};
}

int extend_candidates(const vector<int> &seq, const vector<int> &mask, vector<int> &candidate_locs, int n_matching)
{
    int j = n_matching;

    while (candidate_locs.size() > 1 && j < MAX_K)
    {
        int best = 4;
        vector<int> next_candidates;

        for (int loc : candidate_locs)
        {
            if (loc + j < seq.size())
            {
                int bm = (seq[loc + j] - mask[j] + 4) % 4;
                if (bm < best)
                {
                    best = bm;
                    next_candidates = {loc};
                }
                else if (bm == best)
                {
                    next_candidates.push_back(loc);
                }
            }
        }

        if (next_candidates.empty())
            break;
        candidate_locs = next_candidates;
        ++j;
    }

    return candidate_locs[0];
}

uint64_t hash_val(const vector<int> &seq, const vector<int> &mask, int lex_first_idx)
{
    uint64_t val = 0;
    for (int i = lex_first_idx; i < lex_first_idx + MAX_K && i < seq.size(); ++i)
    {
        val = val * 4 + (seq[i]);
    }

    int base_extend = max(0, MAX_K - (int)(seq.size() - lex_first_idx));
    if (base_extend > 0)
    {
        for (int i = (int)(mask.size()) - base_extend; i < (int)(mask.size()); ++i)
        {
            val = val * 4 + (mask[i] ^ 3);
        }
    }

    return (int)(val);
}

uint64_t lexicographic_first(const vector<int> &seq, const vector<int> &mask)
{
    // cout<<"reached inside lexicographic_first\n";
    auto [candidate_locs, n_matching] = get_candidate_locs(seq, mask);
    int lex_first_idx = extend_candidates(seq, mask, candidate_locs, n_matching);
    return hash_val(seq, mask, lex_first_idx);
}

string convert_hash_to_string(uint64_t hash, vector<int> mask)
{
    string result;
    int i = mask.size() - 1;
    while (hash > 0)
    {
        int digit = hash % 4;
        if (digit == 0)
            result += mask[i];
        else
            result += '_';

        hash /= 4;
        i--;
    }
    reverse(result.begin(), result.end());
    return result;
}

uint64_t naivelexicfirst(const vector<int> &seq, const vector<int> &mask)
{
    uint64_t min_hash = UINT64_MAX;
    uint64_t rolling_hash = 0;
    uint64_t hash = 0;
    uint64_t finans = 0;

    int l = mask.size();
    string minstring;
    // cout << "l: " << l << endl;
    string s = "";
    for (int i = 0; i < l; i++)
    {
        rolling_hash = rolling_hash * 4 + (mask[i] ^ seq[i]);
        hash = hash * 4 + seq[i];
        if (seq[i] == 0)
            s += 'A';
        else if (seq[i] == 1)
            s += 'T';
        else if (seq[i] == 2)
            s += 'G';
        else if (seq[i] == 3)
            s += 'C';
    }

    // cout << "rolling_hash: " << rolling_hash << " " << convert_hash_to_string(rolling_hash, mask) << endl;
    min_hash = min(min_hash, rolling_hash);

    for (int i = 1; i <= seq.size() - l; i++)
    {
        rolling_hash = 0;
        hash = 0;
        // cout<< "len of s " << s.length() << endl;
        s = "";
        for (int j = 0; j < l; j++)
        {
            rolling_hash = rolling_hash * 4 + (mask[j] ^ seq[i + j]);
            hash = hash * 4 + seq[i + j];
            if (seq[i + j] == 0)
                s += 'A';
            else if (seq[i + j] == 1)
                s += 'T';
            else if (seq[j + i] == 2)
                s += 'G';
            else if (seq[i + j] == 3)
                s += 'C';
        }
        if (rolling_hash < min_hash)
        {
            min_hash = rolling_hash;
            minstring = s;
            finans = hash;

            // cout << "rolling_hash: " << rolling_hash << " " << convert_hash_to_string(rolling_hash, mask) << endl;
        }
    }
    // cout << "minstring: " << minstring << endl;
    return finans;
}
/*
int main()

{
    string dna_sequence = "ACGTACGTGCTACATAGCTAGTCAG"; // Example sequence
    vector<int> int_seq = to_int_seq(dna_sequence);

    // Generate one random mask
    vector<vector<int>> masks = get_masks(1, 6);
    vector<int> mask = masks[0];

    for (int i = 0; i < mask.size(); ++i)
    {
        cout << mask[i] << " ";
    }
    cout << endl;
    for (int i = 0; i < mask.size(); i++)
    {

        if (mask[i] == 0)
            cout << 'A';
        else if (mask[i] == 1)
            cout << 'T';
        else if (mask[i] == 2)
            cout << 'G';
        else if (mask[i] == 3)
            cout << 'C';
    }
    cout << endl;
    // Compute lexicographic hash for original sequence
    int hash_original = lexicographic_first(int_seq, mask);

    int hash_naive = naivelexicfirst(int_seq, mask);

    // Compute for reverse complement
    string rc_sequence = reverse_complement(dna_sequence);
    vector<int> int_rc_seq = to_int_seq(rc_sequence);
    // int hash_rc = lexicographic_first(int_rc_seq, mask);

    cout << "Original sequence: " << dna_sequence << endl;
    cout << "Reverse complement: " << rc_sequence << endl;
    cout << "Hash (original): " << hash_original << endl;
    cout << "Hash (naive): " << hash_naive << endl;
    cout << "Hash (naive): " << convert_hash_to_string(hash_naive, mask) << endl;

    // cout << "Hash (rev comp): " << hash_rc << endl;
    // cout << "Min Hash (canonical): " << min(hash_original, hash_rc) << endl;

    return 0;
}*/

pair<vector<uint64_t>, vector<uint64_t>> get_seq_sketch(const string &seq,const vector<vector<int>> &masks)
{
    auto seq_int = to_int_seq(seq);
    vector<uint64_t> sketch, sketch_rc;

    // cout<<"reached inside get_seq_sketch\n";

    for (const auto &mask : masks)
    {
        sketch.push_back(lexicographic_first(seq_int, mask)); // Assuming one hash value per mask
    }
    if (RC)
    {
        auto rc_seq = reverse_complement(seq);
        auto seq_int_rc = to_int_seq(rc_seq);
        for (const auto &mask : masks)
        {
            sketch_rc.push_back(lexicographic_first(seq_int_rc, mask));
        }
    }
    return {sketch, sketch_rc};
}

pair<vector<vector<uint64_t>>, int> sketching(const vector<string> &seqs,const vector<vector<int>> &masks,
          int n_threads, bool rc, int max_k, int min_k)
{

    // shared_masks = masks;
    RC = rc;

    vector<vector<uint64_t>> sketches(seqs.size());
    vector<vector<uint64_t>> sketches_rc;

    int n_seq = seqs.size();
    cout << "sketching " << n_seq << " sequences\n";
    sketches.resize(n_seq);
    if (rc)
        sketches_rc.resize(n_seq);

    // write seq number to file

    // ofstream seq_num_file("seq_num.txt");
    // if (!seq_num_file.is_open())
    // {
    //     cerr << "Error opening file for writing sequence numbers." << endl;
    //     return {};
    // }

#pragma omp parallel for num_threads(n_threads) schedule(static)
    for (int i = 0; i < n_seq; ++i)
    {
        // cout<<"sketching sequence "<<i<<endl;
        // seq_num_file << i << endl;
        if (seqs[i].size() < max_k)
        {
            // cout<<"sketching sequence "<<i<<" with size less than max_k "<<max_k<<endl;
            continue;
        }
        auto [sketch, sketch_rc] = get_seq_sketch(seqs[i], masks);
        sketches[i] = sketch;
        if (rc)
            sketches_rc[i] = sketch_rc;
    }

    if (rc)
    {
        // Concatenate original and RC sketches
        sketches.insert(sketches.end(), sketches_rc.begin(), sketches_rc.end());
    }

    int cnt = 0;
    for (int i = 0; i < sketches.size(); i++)
    {
        if (sketches[i].size() == 0)
        {
            cnt++;
        }
    }
    cout << "sketches with 0 size: " << cnt << endl;

    return {sketches, rc ? 2 * n_seq : n_seq};
}

/*
int main() {
    int n_cpu = 4;
    bool rc = false;
    int max_k = 5, min_k = 3;

    // Dummy sequences (name, sequence)
    vector<pair<string, string>> seqs = {
        {"read1", "ACGTACGT"},
        {"read2", "TGCATGCA"},
        {"read3", "GGCCATTA"}
    };

    // Dummy masks (like hash functions)
    vector<string> seqs_str;
    for (const auto& seq : seqs) {
        seqs_str.push_back(seq.second);
    }

    vector<vector<int>> masks = get_masks(10, max_k);
    for(auto mask : masks){
        for(auto x : mask)
        {
            if(x==0)
                cout<<'A';
            else if(x==1)
                cout<<'T';
            else if(x==2)
                cout<<'G';
            else if(x==3)
                cout<<'C';
        }
        cout<<endl;
    }


    auto [sketches, n_seq] = sketching(seqs_str, masks, n_cpu, rc, max_k, min_k);

    // Output
    cout << "Total sequences: " << n_seq << "\n";
    for (int i = 0; i < sketches.size(); ++i) {
        cout << "Sketch " << i + 1 << ": ";
        for (auto x : sketches[i]) cout << x << " ";
        cout << "\n";
    }

    return 0;
}


*/

// using vector<vector<uint64_t>> = vector<vector<uint64_t>>;
// using MatchSets = map<int, vector<vector<int>>>; // match length -> list of groups
// using vector<map<int, vector<vector<int>>>> = vector<map<int, vector<vector<int>>>>;
// using pair<int, pair<int, char>> = pair<int, pair<int, char>>; // (_i1, (_i2, strand))

// Forward declarations
vector<map<int, vector<vector<int>>>> prefix_tree_multiproc(const vector<vector<uint64_t>> &sketches, int n_hash, int max_k, int min_k, int n_cpu);
map<int, vector<vector<int>>> get_matching_sets(const vector<vector<uint64_t>> &sketches, int sketch_idx, int max_k, int min_k);
map<pair<int, pair<int, char>>, int> bottom_up_pair_aggregation(
    const vector<map<int, vector<vector<int>>>> &all_matching_sets,
    int n_seq,
    int n_pairs,
    int max_k,
    int min_k);

// Main pairwise comparison function
map<pair<int, pair<int, char>>, int> pairwise_comparison(
    const vector<vector<uint64_t>> &sketches,
    int n_seq,
    int n_hash,
    int n_cpu,
    float alpha,
    int max_k,
    int min_k)
{
    int n_pairs = (alpha > 0) ? static_cast<int>(n_seq * alpha) : INT_MAX;
    cout<<"n_seq: " << n_seq<<"alpha "<<alpha << endl;
    cout<<"n_pairs: " << n_pairs << endl;
    vector<map<int, vector<vector<int>>>> all_matching_sets = prefix_tree_multiproc(sketches, n_hash, max_k, min_k, n_cpu);
    return bottom_up_pair_aggregation(all_matching_sets, n_seq, n_pairs, max_k, min_k);
}

vector<map<int, vector<vector<int>>>> prefix_tree_multiproc(const vector<vector<uint64_t>> &sketches, int n_hash, int max_k, int min_k, int n_cpu)
{
    vector<map<int, vector<vector<int>>>> all_sets(n_hash); 
    vector<thread> threads;

    
    int current_idx = 0;
    mutex idx_mutex;

    auto worker = [&]()
    {
        while (true)
        {
            int i;
            {
                lock_guard<mutex> lock(idx_mutex);
                if (current_idx >= n_hash)
                    break;
                i = current_idx++;
            }
            all_sets[i] = get_matching_sets(sketches, i, max_k, min_k);
        }
    };

    for (int t = 0; t < n_cpu; ++t)
    {
        threads.emplace_back(worker);
    }

    for (auto &t : threads)
    {
        t.join();
    }

    return all_sets;
}


map<int, vector<vector<int>>> get_matching_sets(const vector<vector<uint64_t>> &sketches, int sketch_idx, int max_k, int min_k)
{
    int n = sketches.size();
    vector<vector<int>> subtrees = {vector<int>(n)};
    iota(subtrees[0].begin(), subtrees[0].end(), 0); 

    map<int, vector<vector<int>>> matching_sets;

    for (int k = 0; k < max_k; ++k)
    {
        vector<vector<int>> next_subtrees;

        for (const auto &seq_idxs : subtrees)
        {
            vector<int> partitions[4]; // DNA bases: 0, 1, 2, 3
            for (int idx : seq_idxs)
            {
                if (sketches[idx].size() <= sketch_idx)
                    continue; // Skip if no sketch
                uint64_t val = sketches[idx][sketch_idx];
                int ch = (val >> (2 * (max_k - k - 1))) & 3;
                partitions[ch].push_back(idx);
            }
            for (int i = 0; i < 4; ++i)
            {
                if (partitions[i].size() > 1)
                {
                    next_subtrees.push_back(partitions[i]);
                }
            }
        }

        if (next_subtrees.empty())
            break;
        subtrees = next_subtrees;

        if (k + 1 >= min_k)
        {
            matching_sets[k + 1] = subtrees;
        }
    }

    return matching_sets;
}

// Aggregate all pairwise max-match-lengths with strand awareness
map<pair<int, pair<int, char>>, int> bottom_up_pair_aggregation(
    const vector<map<int, vector<vector<int>>>> &all_matching_sets,
    int n_seq,
    int n_pairs,
    int max_k,
    int min_k)
{

    map<pair<int, pair<int, char>>, int> pair_match_lens;

    for (int k = max_k; k >= min_k; --k)
    {
        for (const auto &matching_sets : all_matching_sets)
        {
            auto it = matching_sets.find(k);
            if (it == matching_sets.end())
                continue;

            for (const auto &group : it->second)
            {
                for (int i1 : group)
                {
                    for (int i2 : group)
                    {
                        int _i1 = i1 % n_seq;
                        int _i2 = i2 % n_seq;
                        if (_i2 > _i1)
                        {
                            char strand = '+';
                            if((i1>=n_seq && i2<n_seq) || (i1<n_seq && i2>=n_seq))
                            {
                                strand = '-';
                            }
                            pair<int, pair<int, char>> key = {_i1, {_i2, strand}};
                            if (pair_match_lens.find(key) == pair_match_lens.end())
                            {
                                pair_match_lens[key] = k;
                                if (pair_match_lens.size() == static_cast<size_t>(n_pairs))
                                {
                                    return pair_match_lens;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return pair_match_lens;
}

int main(int argc, char *argv[])
{
    string out_dir;
    string fasta_path;
    string method = "lexichash";
    float alpha = 0.1;
    // Simple manual parsing
    for (int i = 1; i < argc; ++i)
    {
        string arg = argv[i];
        if (arg == "--out" && i + 1 < argc)
        {
            out_dir = argv[++i];
        }
        else if (arg == "--fasta" && i + 1 < argc)
        {
            fasta_path = argv[++i];
        }
        else if (arg == "--method" && i + 1 < argc)
        {
            method = argv[++i];
        }
        // Add more args as needed
        else if (arg == "--n_cpu" && i + 1 < argc)
        {
            n_cpu = stoi(argv[++i]);
        }
        else if (arg == "--min_k" && i + 1 < argc)
        {
            min_k = stoi(argv[++i]);
        }
        else if (arg == "--max_k" && i + 1 < argc)
        {
            max_k = stoi(argv[++i]);
        }
        else if (arg == "--n_hash" && i + 1 < argc)
        {
            n_hash = stoi(argv[++i]);
        }
        else if (arg =="--alpha" && i + 1 < argc)
        {
            alpha = stof(argv[++i]);
        }
        else if (arg == "--no_rc")
        {
            RC = false;
        }
        
        
        else if (arg == "--help")
        {
            cout << "Usage: " << argv[0] << " --out <dir> --fasta <file> [--method <lexichash|minhash>]\n";
            return 0;
        }
        // Add more args as needed
    }

    // Check required
    if (out_dir.empty() || fasta_path.empty())
    {
        cerr << "Usage: " << argv[0] << " --out <dir> --fasta <file> [--method <lexichash|minhash>]\n";
        return 1;
    }

    // Validation
    if (method != "lexichash" && method != "minhash")
    {
        cerr << "Error: method must be 'lexichash' or 'minhash'\n";
        return 1;
    }

    // Output to confirm
    cout << "out_dir = " << out_dir << "\n";
    cout << "fasta_path = " << fasta_path << "\n";
    cout << "method = " << method << "\n";
    cout << "n_cpu = " << n_cpu << "\n";
    cout << "min_k = " << min_k << "\n";
    cout << "max_k = " << max_k << "\n";
    cout << "n_hash = " << n_hash << "\n";
    cout << "RC = " << (RC ? "true" : "false") << "\n";
    cout<< "alpha = " << alpha << "\n";
    // Start timing

    auto start_program = chrono::high_resolution_clock::now();

    // Read FASTA file
    vector<string> seqs;
    
    seqs= read_fasta(fasta_path, max_k);

    cout << "Read " << seqs.size() << " sequences from " << fasta_path << "\n";
    // Generate masks

    cout << "Generating masks...\n\n\n";

    
    vector<vector<int>> masks = get_masks(n_hash, max_k);
    shared_masks = masks;
    //print masks to a file
    ofstream mask_file(out_dir + "/masks.txt");
    if (!mask_file.is_open())
    {
        cerr << "Error: Could not open output file " << out_dir + "/masks.txt" << "\n";
        return 1;
    }
    for (const auto &mask : masks)
    {
        for (int i = 0; i < mask.size(); i++)
        {
            if (mask[i] == 0)
                mask_file << 'A';
            else if (mask[i] == 1)
                mask_file << 'T';
            else if (mask[i] == 2)
                mask_file << 'G';
            else if (mask[i] == 3)
                mask_file << 'C';
        }
        mask_file << endl;
    }


    cout << "Generated " << masks.size() << " masks\n";
    auto start = chrono::high_resolution_clock::now();
    auto [sketches, n_seq] = sketching(seqs, masks, n_cpu, RC, max_k, min_k);

    cout << "# sequences: " << n_seq << endl;

    //    utils.print_clocktime(start, perf_counter(), 'sketching')
    cout << "Sketching took: " << chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start).count() << " ms\n";

    cout << "Pairwise comparison...\n\n\n";
    // Pairwise comparison
    auto start_pairwise = chrono::high_resolution_clock::now();

    auto pairwise_results = pairwise_comparison(sketches, n_seq, n_hash, n_cpu, alpha, max_k, min_k);
    cout << "Pairwise comparison took: " << chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start_pairwise).count() << " ms\n";

    // Output results
    ofstream out_file(out_dir + "/pairwise_results.txt");
    if (!out_file.is_open())
    {
        cerr << "Error: Could not open output file " << out_dir + "/pairwise_results.txt" << "\n";
        return 1;
    }
    for (const auto &pair : pairwise_results)
    {
        out_file << pair.first.first << "\t" << pair.first.second.first << "\t" << pair.first.second.second << "\t" << pair.second << "\n";
    }
    out_file.close();

    cout << "Pairwise results written to " << out_dir + "/pairwise_results.txt" << "\n";
    // End timing
    auto end_program = chrono::high_resolution_clock::now();

    cout << "Total program took: " << chrono::duration_cast<chrono::milliseconds>(end_program - start_program).count() << " ms\n";
}