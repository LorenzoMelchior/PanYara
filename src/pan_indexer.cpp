#include <stdexcept>
#include <cstdlib>

#include "pibf.h"



int main(int argc, char const ** argv)
{

    ///////////////// Read reference and VCF files ////////////////////////

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <reference.fasta> <vcf_folder_path> <read_length>" << std::endl;
        return 1;
    }

    std::string reference_filename = argv[1];

    std::string folder_path = argv[2];
    std::vector<std::string> vcf_files;

    std::size_t read_length = std::stoi(argv[3]);

    try {
        for (const auto& entry : std::filesystem::directory_iterator(folder_path)) 
        {
            if (entry.is_regular_file() && entry.path().extension() == ".vcf") 
            {
                std::string file_path = entry.path().string();
                vcf_files.push_back(file_path);
            }
        }
    } catch (const std::filesystem::filesystem_error& e) 
    {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
        return 1;
    }


    ////////////////// Create indexing folder structure /////////////////////

    std::string index_folder_path = "index";

    if (argc > 4) {
        index_folder_path = argv[4];
    }

    if (std::filesystem::exists(index_folder_path)) 
    {
        std::cout << "Index folder already exists. Do you want to override it (y/n)?" << std::endl;
        char answer;
        std::cin >> answer;
        if (answer == 'y') {
            std::filesystem::remove_all(index_folder_path);
        } else {
            std::cout << "Aborting." << std::endl;
            return 1;
        }
    }

    std::filesystem::create_directory(index_folder_path);

    std::filesystem::create_directory(index_folder_path + "/reference");

    for (std::size_t i = 1; i < vcf_files.size(); i++) 
    {
        std::filesystem::create_directory(index_folder_path + "/vcf_" + std::to_string(i));
    }

    ///////////////// Index the PIBF ////////////////////////



    pangenomic_hibf pibf{4u, 3u, 8*8192u, vcf_files.size()};

    pibf.feed_reference(reference_filename);

    for (std::size_t index = 0; index < vcf_files.size(); ++index) 
    {
        pibf.feed_vcf(reference_filename, vcf_files[index], index);
    }

    pibf.save_to_disk(index_folder_path);


    /////////////////// Index Dream-Yara ////////////////////////


    auto ref_reader = ivio::fasta::reader{{.input = reference_filename}};
    std::string reference;

    for (auto record_view : ref_reader) {
        reference += record_view.seq;
    }

    // MUst convert filename from char to string!
    std::string reference_folder_path = index_folder_path + "/reference";
    std::string reference_filename_string = reference_filename;
    execute_yara_indexer(reference_folder_path, reference_filename_string);


    std::size_t i = 1;
    for (const std::string& vcf: vcf_files)
    {
        std::string vcf_sequence = reference;
        std::vector<IndexRange> ranges;
        for (auto view : ivio::vcf::reader{{vcf}}) 
        {
            std::size_t position = view.pos;
            std::size_t start_position;

            vcf_sequence.replace(position, 1, view.alts);

            if (position > read_length-1) {
                start_position = position - (read_length-1);
            } else {
                start_position = 0;
            }
            
            
            std::size_t end_position = position + read_length;
            IndexRange range = {start_position, end_position};
            ranges.push_back(range);
        }
        std::vector<IndexRange> merged_ranges = merge_ranges(ranges);


        std::filesystem::create_directory(index_folder_path + "/vcf_" + std::to_string(i));
        for (IndexRange range: merged_ranges)
        {
            std::string folder_path = index_folder_path + "/vcf_" + std::to_string(i) + "/" + std::to_string(range.begin) + "_" + std::to_string(range.end);
            std::filesystem::create_directory(folder_path);

            std::string filename = "vcf_" + std::to_string(i) + "_" + std::to_string(range.begin) + "_" + std::to_string(range.end) + ".fasta";
            write_vcf_fastas(vcf_sequence, folder_path, filename, range);


            std::string fasta_input = folder_path + "/*.fasta";
            execute_yara_indexer(folder_path, fasta_input);

        }
        i++;
    }



    return 0;
}