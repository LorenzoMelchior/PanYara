#include <iostream>
#include <filesystem>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <fstream>
#include <string>
#include <deque>
#include <functional>
#include <cereal/archives/binary.hpp>
#include <ivio/ivio.h>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <chrono>


struct pibf_results 
{
    bool allVariants = false;
    seqan3::interleaved_bloom_filter<seqan3::uncompressed>::membership_agent_type::binning_bitvector positives;
    seqan3::interleaved_bloom_filter<seqan3::uncompressed>::membership_agent_type::binning_bitvector negatives;
};

struct IndexRange 
{
        std::size_t begin;
        std::size_t end;
};

void generate_kmers(const std::string &filename, unsigned int k, const std::function<void(const std::string&)>& callback);
void generate_vcf_kmers(const std::string &ref_filename, 
                        const std::string &vcf_filename, 
                        unsigned int k, 
                        const std::function<void(const std::vector<std::string>&)>& callback1);
size_t size_of_bit_array_estimation(double p, int n);
std::vector<IndexRange> merge_ranges(std::vector<IndexRange> ranges);
void write_vcf_fastas(const std::string &ref_filename, 
                      const std::string &vcf_filename, 
                      const std::string &output_folder, 
                      const IndexRange &ranges);
void execute_yara_indexer(const std::string& folder, const std::string& fasta_file);
void execute_yara_mapper(const std::string& output_bam, const std::string& index_dir, const std::string& reads_fasta);
void execute_dream_yara_indexer(const std::string& folder, const std::string& fasta_file);
void execute_dream_yara_mapper(const std::string& output_bam, const std::string& index_dir, const std::string& reads_fasta);


class pangenomic_hibf 
{

    private:
        std::size_t k;
        std::size_t number_of_hash_funcs;
        std::size_t bit_array_size;
        std::size_t number_of_bins;
        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> all_variants;
        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> positives;
        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> negatives;
        

    public:
        pangenomic_hibf(std::size_t number_of_hash_funcs, std::size_t bit_array_size, std::size_t number_of_vcf_files)
                : all_variants{seqan3::bin_count{1u}, seqan3::bin_size{bit_array_size}},
                   positives{seqan3::bin_count{number_of_vcf_files}, seqan3::bin_size{bit_array_size/200}}, 
                   negatives{seqan3::bin_count{number_of_vcf_files}, seqan3::bin_size{bit_array_size/200}} {
            this->k = 21u;
            this->number_of_hash_funcs = number_of_hash_funcs;
            this->bit_array_size = bit_array_size;

        }

        std::size_t get_k() {
            return this->k;
        }

        void feed_reference(const std::string &reference){
            generate_kmers(reference, this->k,[this](const std::string& kmer) {
                this->all_variants.emplace(std::hash<std::string>{}(kmer), seqan3::bin_index{0u});
            });
        };

        void feed_vcf(const std::string &reference, const std::string &vcf, std::size_t index) {
            generate_vcf_kmers(reference, vcf, this->k,[this, index](const std::vector<std::string>& kmers) {
                this->all_variants.emplace(std::hash<std::string>{}(kmers[0]), seqan3::bin_index{0u});
                this->positives.emplace(std::hash<std::string>{}(kmers[0]), seqan3::bin_index{index});
                this->negatives.emplace(std::hash<std::string>{}(kmers[1]), seqan3::bin_index{index});
            });

        }

        void save_to_disk(const std::string &path) {
            std::filesystem::create_directory(path + "/pibf");
            std::ofstream os{path + "/pibf/pibf.bin", std::ios::binary};
            cereal::BinaryOutputArchive oarchive{os};
            oarchive(this->all_variants, this->positives, this->negatives);
        }

        void load_from_disk(const std::string &path) {
            std::ifstream is{path + "/pibf/pibf.bin", std::ios::binary};
            cereal::BinaryInputArchive iarchive{is};
            iarchive(this->all_variants, this->positives, this->negatives);
        }

        inline const pibf_results query_kmer(const std::string &kmer) {

            pibf_results results;

            auto av_agent = this->all_variants.membership_agent();
            auto av_result = av_agent.bulk_contains(std::hash<std::string>{}(kmer));

           if (!av_result[0]) {
                // all variants does not contain kmer -> return 0 for all bins (all_variants (1) + positives / negatives)
                results.allVariants = false;
                return results;
            } else {
                results.allVariants = true;
            }

            auto p_agent = this->positives.membership_agent();
            auto p_result = p_agent.bulk_contains(std::hash<std::string>{}(kmer));
            results.positives = p_result;

            auto n_agent = this->negatives.membership_agent();            
            auto n_result = n_agent.bulk_contains(std::hash<std::string>{}(kmer));
            results.negatives = n_result;

            return results;
        }
};
