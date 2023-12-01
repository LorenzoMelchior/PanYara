#include "pibf.h"


//using namespace seqan3::literals;


void generate_kmers(const std::string &filename, unsigned int k, const std::function<void(const std::string&)>& callback) {
    std::ifstream input(filename);
    if (!input.is_open()) {
        throw std::runtime_error("Cannot open file");
    }

    std::deque<char> window;
    char c;
    bool is_sequence_line = false;
    while (input.get(c)) {
        if (c == '\n') {  // End of line, check next line
            is_sequence_line = true;
            continue;
        }
        
        if (is_sequence_line) {
            window.push_back(c);
            if (window.size() == k) {  // When window is full
                if(std::count(window.begin(), window.end(), 'N') == 0) { // No 'N' in the k-mer
                    callback(std::string(window.begin(), window.end()));
                }
                window.pop_front();
            }
        }
    }
}

void generate_vcf_kmers(const std::string &ref_filename, 
                        const std::string &vcf_filename, 
                        unsigned int k, 
                        const std::function<void(const std::vector<std::string>&)>& callback) 
{
    std::ifstream ref_input(ref_filename);
    if (!ref_input.is_open()) {
        throw std::runtime_error("Cannot open reference file");
    }
    
    std::vector<char> ref_sequence;
    char c;
    bool is_sequence_line = false;
    while (ref_input.get(c)) {
        if (c == '\n') {  // End of line, check next line
            is_sequence_line = true;
            continue;
        }
        
        if (is_sequence_line) {
            ref_sequence.push_back(c);
        }
    }
    
    std::size_t view_pos_unsigned;  // Unsigned version of view.pos
    for (auto view : ivio::vcf::reader{{vcf_filename}}) {
        view_pos_unsigned = static_cast<std::size_t>(view.pos);
        if(view_pos_unsigned < k - 1) {
            continue;  // Skip if the k-mer starting at the position would be out of bounds
        }

        for (size_t i = 0; i < view.alts.size(); ++i) {
            std::vector<std::string> kmers;
            const auto& alt = view.alts[i];
            const auto& ref = view.ref[i];
            for(unsigned int offset = 0; offset < k && view_pos_unsigned - offset + k <= ref_sequence.size(); ++offset) {
                std::string kmer(ref_sequence.begin() + view.pos - offset, ref_sequence.begin() + view.pos - offset + k);
                std::string destroyed_kmer = kmer;
                kmer[offset] = alt;  // Replace the corresponding position in the k-mer with the alternate nucleotide
                destroyed_kmer[offset] = ref;
                kmers.push_back(kmer);
                kmers.push_back(destroyed_kmer);
                callback(kmers);
            }
        }
        
    }
}

size_t size_of_bit_array_estimation(double p, int n) {
    // Computes the size of a BloomFilter BitArray for a certain given probability $p of false positives
    // and the number of items in the filter $n. 
    // Source for formula: https://hur.st/bloomfilter/
    double result = std::ceil((n * std::log(p)) / std::log(1 / (std::pow(2, std::log(2)))));
    if(result < 0 || result > std::numeric_limits<size_t>::max()) {
        throw std::out_of_range("The result exceeds the range of unsigned int (size_t).");
    }
    return static_cast<size_t>(result);
}


std::vector<IndexRange> merge_ranges(std::vector<IndexRange> ranges) 
{
    // Merges overlapping ranges in a vector of ranges.
    std::vector<IndexRange> merged_ranges;
    std::sort(ranges.begin(), ranges.end(), [](const IndexRange& lhs, const IndexRange& rhs) { return lhs.begin < rhs.begin; });
    for (const auto& range : ranges) {
        if (merged_ranges.empty() || merged_ranges.back().end < range.begin) {
            merged_ranges.push_back(range);
        } else {
            merged_ranges.back().end = std::max(merged_ranges.back().end, range.end);
        }
    }
    return merged_ranges;
}


void write_vcf_fastas(const std::string &vcf_sequence, 
                         const std::string &output_folder, 
                         const std::string &filename,
                         const IndexRange &range) 
{
    /// Here there only should be written a fasta file for a given vcf-sequence (fasta+vcf applied as string) and the range. Maybe no func needed? 
    std::string sequence = vcf_sequence.substr(range.begin, range.end - range.begin);

    std::ofstream file(output_folder + "/" + filename + ".fasta", std::ios::out | std::ios::binary);
    if (!file) {
        throw std::runtime_error("Failed to open the file.");
    }

    file << ">" << filename << "\n";
    file << sequence << "\n";
    file.close();
}

void execute_yara_indexer(const std::string& folder, const std::string& fasta_file)
{
    std::string yara_filter = "bin/yara_indexer " + fasta_file + " -o " + folder + "/REF.index";
    std::cout << yara_filter << std::endl;
    int cmd = system(yara_filter.c_str());

    if (cmd != 0) {
        throw std::runtime_error("Failed to execute yara_indexer.");
    }
    
}

void execute_yara_mapper(const std::string& output_bam, const std::string& index_dir, const std::string& reads_fasta)
{
    std::string yara_args = "bin/yara_mapper " + index_dir + "/REF.index " + reads_fasta + " -o " + output_bam;
    int cmd = system(yara_args.c_str());

    if (cmd != 0) {
        throw std::runtime_error("Failed to execute yara_mapper.");
    }
}

void execute_dream_yara_indexer(const std::string& folder, const std::string& fasta_file) 
{
    std::string yara_filter = "bin/dream_yara_build_filter --threads 8 --kmer-size 18 --filter-type bloom --bloom-size 16 --num-hash 3 --output-file " + folder + "/IBF.filter " + fasta_file;
    std::string yara_indexer = "bin/dream_yara_indexer --threads 8 --output-prefix " + folder + "/ " + fasta_file;
    int cmd1 = system(yara_filter.c_str());
    int cmd2 = system(yara_indexer.c_str());

    if (cmd1 != 0 || cmd2 != 0) {
        throw std::runtime_error("Failed to execute yara_indexer.");
    }

}

void execute_dream_yara_mapper(const std::string& output_bam, const std::string& index_dir, const std::string& reads_fasta) 
{
    std::string yara_args = "bin/dream_yara_mapper -t 8 -ft bloom -e 0.03 -fi " + index_dir + "IBF.filter -o " + output_bam + " " + index_dir + " " + reads_fasta;
    int cmd = system(yara_args.c_str());

    if (cmd != 0) {
        throw std::runtime_error("Failed to execute yara_mapper.");
    }
}