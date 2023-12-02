#include "pibf.h"
#include <argparse.hpp>

int main(int argc, char const ** argv)
{

    /////////////////// Parse arguments ////////////////////////

    argparse::ArgumentParser program("pan_yara_mapper");

    program.add_argument("-r", "--reads")
        .help("Read file location")
        .required();

    program.add_argument("-i", "--index_folder")
        .help("Index folder location")
        .default_value("index");

    program.add_argument("--dream")
        .help("Use Dream-Yara instead of Yara")
        .default_value(false)
        .implicit_value(true);

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cout << err.what() << std::endl;
        std::cout << program;
        exit(1);
    }

    std::string reads_fasta = program.get<std::string>("--reads");
    std::string index_dir = program.get<std::string>("--index_folder");
    const bool dream = program.get<bool>("--dream");


    /////////////////// Read VCF files ////////////////////////

    std::size_t vcf_count = 0;
    for (const auto& entry : std::filesystem::directory_iterator(index_dir)) {
    // Check if the entry is a directory and has the specified prefix _vcf
        if (std::filesystem::is_directory(entry) && entry.path().filename().string().rfind("vcf_", 0) == 0) {
            vcf_count++;
        }
    }

    const std::string output_bam = "results.bam";

    /////////////////// Load PIBF ////////////////////////

    pangenomic_hibf pibf{3u, 8*8192u, vcf_count};
    std::size_t kmer_size = pibf.get_k();
    pibf.load_from_disk(index_dir);

    ///////////////// Execute Mapping //////////////////////
    
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
            bool need_ref = false;
            if (is_in_ref[i] >=  is_in_var[i]) 
            {
                if (!need_ref) {
                    need_ref = true;
                    if (dream) {
                        execute_dream_yara_mapper(index_dir + "/reference/" + output_bam, index_dir + "/reference/", reads_fasta);
                    } else {
                        execute_yara_mapper(index_dir + "/reference/" + output_bam, index_dir + "/reference/", reads_fasta);
                    }
                }
                std::filesystem::copy_file(index_dir + "/reference/" + output_bam, index_dir + "/vcf_" + std::to_string(i+1) + "/" + output_bam, std::filesystem::copy_options::overwrite_existing);
            }
            else if (is_in_var[i] >= k_mer_count*0.5) 
            {
                for (const auto& entry : std::filesystem::directory_iterator(index_dir + "/vcf_" + std::to_string(i+1))) 
                {
                    if (dream) {
                        execute_dream_yara_mapper(entry.path().string() + "/" + output_bam, entry.path().string() + "/", reads_fasta);
                    } else {
                        execute_yara_mapper(entry.path().string() + "/" + output_bam, entry.path().string() + "/", reads_fasta);
                    }
                }
            }

        }



    }




    return 0;
}