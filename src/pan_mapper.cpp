#include "pibf.h"


int main(int argc, char const ** argv)
{
    // Args: Output.bam, Index_dir, Reads.fasta, kmer_size
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <output.bam> <index_dir> <reads.fasta> <kmer_size>" << std::endl;
        return 1;
    }

    std::string output_bam = argv[1];
    std::string index_dir = argv[2];
    std::string reads_fasta = argv[3];
    std::size_t kmer_size = std::stoi(argv[4]);

    std::size_t vcf_count = 0;
    for (const auto& entry : std::filesystem::directory_iterator(index_dir)) {
    // Check if the entry is a directory and has the specified prefix _vcf
        if (std::filesystem::is_directory(entry) && entry.path().filename().string().rfind("vcf_", 0) == 0) {
            vcf_count++;
        }
    }

    pangenomic_hibf pibf{4u, 3u, 8*8192u, vcf_count};
    pibf.load_from_disk(index_dir);

    

    
    
    auto reader = ivio::fasta::reader{{.input = reads_fasta, .compressed = false}};
    std::size_t is_in_var[vcf_count] = {};
    std::size_t is_in_ref[vcf_count] = {};

    for (auto record_view : reader) {

        std::string read(record_view.seq);
        std::size_t k_mer_count = 0;

        for (std::size_t i = 0; i <= read.length() - kmer_size; ++i) {
            std::string kmer = read.substr(i, kmer_size);

            pibf_results results = pibf.query_kmer(kmer);

            if (results.allVariants) 
            {
                for (std::size_t i = 0; i < vcf_count; i++) 
                {
                    if (results.positives[i] == results.negatives[i]) 
                    {
                        is_in_var[i]++;
                        is_in_ref[i]++;
                    }
                    else if (results.negatives[i]) {
                        is_in_ref[i]++;
                    }
                    else if (results.positives[i]) {
                        is_in_var[i]++;
                    }
                }
            }
            k_mer_count++;
        }    

        for (std::size_t i = 0; i < vcf_count; i++) 
        {
            std::size_t need_ref = 0;
            if (is_in_ref[i] >=  is_in_var[i]) 
            {
                if (!need_ref) {
                    need_ref = 1;
                    std::cout << 1 << std::endl;
                    execute_yara_mapper(index_dir + "/reference/" + output_bam, index_dir + "/reference/", reads_fasta);
                }
                std::filesystem::copy_file(index_dir + "/reference/" + output_bam, index_dir + "/vcf_" + std::to_string(i+1) + "/" + output_bam, std::filesystem::copy_options::overwrite_existing);
            }
            else if (is_in_var[i] >= k_mer_count*0.5) 
            {
                for (const auto& entry : std::filesystem::directory_iterator(index_dir + "/vcf_" + std::to_string(i+1))) 
                {
                    std::cout << 2 << std::endl;
                    execute_yara_mapper(entry.path().string() + "/" + output_bam, entry.path().string() + "/", reads_fasta);
                }
            }

        }



    }




    return 0;
}