#include <bitset>
#include <experimental/filesystem>
#include <iostream>
#include <limits>
#include <vector>

#include "../../src/dna.hpp"
#include "../../src/types/CombinePattern.hpp"
#include "../../src/types/IOConfig.hpp"
#include "../../src/utilities.hpp"

#include "gtest/gtest.h"

namespace fs = std::experimental::filesystem;

using namespace priset;

class io_config_test_f : public ::testing::Test {

private:
    std::string fp = __FILE__;
    std::string rp = fp.substr(0, fp.find_last_of('/'));

public:
    fs::path lib_dir = fs::path(rp) / "../library/3041";
    fs::path work_dir = fs::path(rp) / "../work/3041";
    fs::path index_dir = work_dir / "index";
    IOConfig io_cfg;

protected:
    void SetUp() override {
        io_cfg.assign(lib_dir, work_dir, 1);
    }
};

class io_config_test2_f : public ::testing::Test {

private:
    std::string fp = __FILE__;
    std::string rp = fp.substr(0, fp.find_last_of('/'));

public:
    fs::path lib_dir = fs::path(rp) / "../library/one_seq";
    fs::path work_dir = fs::path(rp) / "../work/one_seq";
    fs::path index_dir = work_dir / "index";
    IOConfig io_cfg;

protected:
    void SetUp() override {
        io_cfg.assign(lib_dir, work_dir, 0);
    }
};


TEST(io_config_test, constructor)
{
    IOConfig io_cfg1;
    IOConfig io_cfg2();
    IOConfig io_cfg3{};
    IOConfig io_cfg4{io_cfg3};
    IOConfig io_cfg5 = io_cfg1;
}

TEST(io_config_test, assign)
{
    std::string fp = __FILE__;
    std::string rp = fp.substr(0, fp.find_last_of('/'));
    fs::path lib_dir = fs::path(rp) / "../library/3041";
    fs::path work_dir = fs::path(rp) / "../work/3041";
    IOConfig io_cfg;
    io_cfg.assign(lib_dir, work_dir, 1);
}

// FM_idx_flag = true and idx_dir non-existent
TEST_F(io_config_test_f, FM_idx_flag_1_wo_idx_dir)
{
    EXPECT_EQ(1, io_cfg.get_FM_idx_flag());
}

// FM_idx_flag = true and idx_dir existent => del idx_dir
TEST_F(io_config_test_f, FM_idx_flag_1_w_idx_dir)
{
    char cmd[50];
    sprintf(cmd, "mkdir -p %s", index_dir.c_str());
    system(cmd);
    EXPECT_EQ(1, io_cfg.get_FM_idx_flag());
}

// FM_idx_flag = false and idx_dir non-existent => exit(-1)
TEST_F(io_config_test_f, FM_idx_flag_0_wo_idx_dir)
{
    IOConfig io_cfg_0;
    EXPECT_DEATH(io_cfg_0.assign(lib_dir, work_dir, 0), "");
}

// idx_only = true and skip_idx = false and idx_dir existent
TEST_F(io_config_test_f, FM_idx_flag_0_w_idx_dir)
{
    char cmd[50];
    sprintf(cmd, "mkdir -p %s", index_dir.c_str());
    system(cmd);
    IOConfig io_cfg_0;
    io_cfg_0.assign(lib_dir, work_dir, 0);
    EXPECT_EQ(0, io_cfg_0.get_FM_idx_flag());
}

TEST_F(io_config_test_f, get_acc_file)
{
    fs::path acc_file = io_cfg.get_acc_file();
    EXPECT_EQ("root_3041.acc", acc_file.filename().string());
    EXPECT_TRUE(fs::exists(acc_file));
}

TEST_F(io_config_test_f, get_app_template)
{
    fs::path app_template = io_cfg.get_app_template();
    EXPECT_EQ("app_template.R", app_template.filename().string());
    EXPECT_TRUE(fs::exists(app_template));
}

TEST_F(io_config_test_f, get_fasta_file)
{
    fs::path fasta_file = io_cfg.get_fasta_file();
    EXPECT_EQ("root_3041.fasta", fasta_file.filename().string());
    EXPECT_TRUE(fs::exists(fasta_file));
}

TEST_F(io_config_test2_f, get_index_dir)
{
    fs::path index_dir = io_cfg.get_index_dir();
    std::string index_dir_str = index_dir.string();
    size_t pos = index_dir_str.find_last_of('/');
    EXPECT_EQ("index", index_dir_str.substr(pos + 1, index_dir_str.size() - pos + 1));
    EXPECT_TRUE(fs::exists(index_dir));
}

TEST_F(io_config_test2_f, get_index_txt_path)
{
    std::string index_txt_path_str = io_cfg.get_index_txt_path().string();
    size_t pos = index_txt_path_str.find_last_of('/');
    EXPECT_EQ("index.txt", index_txt_path_str.substr(pos + 1, index_txt_path_str.size() - pos + 1));
    EXPECT_TRUE(fs::exists(index_txt_path_str + ".concat"));
    EXPECT_TRUE(fs::exists(index_txt_path_str + ".limits"));
}

TEST_F(io_config_test_f, get_mapping_dir)
{
    fs::path mapping_dir = fs::canonical(work_dir / "mapping");
    EXPECT_EQ(mapping_dir.string(), io_cfg.get_mapping_dir().string());
    EXPECT_TRUE(fs::exists(io_cfg.get_mapping_dir()));
}

TEST_F(io_config_test_f, get_primer_info_file)
{
    fs::path primer_info_file = io_cfg.get_primer_info_file();
    std::cout << "primer_info_file = " << primer_info_file << std::endl;
    std::string primer_info_file_str = primer_info_file.string();
    size_t pos = primer_info_file_str.find_last_of('/');
    EXPECT_EQ("primer_info.csv", primer_info_file_str.substr(pos + 1, primer_info_file_str.size() - pos + 1));
    EXPECT_TRUE(fs::exists(fs::path(primer_info_file_str.substr(0, pos))));
}

TEST_F(io_config_test_f, get_result_file)
{
    fs::path result_file = io_cfg.get_result_file();
    std::string result_file_str = result_file.string();
    size_t pos = result_file_str.find_last_of('/');
    EXPECT_EQ("result.csv", result_file_str.substr(pos + 1, result_file_str.size() - pos + 1));
    EXPECT_TRUE(fs::exists(result_file.parent_path()));
}

TEST_F(io_config_test_f, get_script_file)
{
    fs::path script_file = io_cfg.get_script_file();
    EXPECT_EQ("app.R", script_file.filename().string());
    EXPECT_TRUE(fs::exists(script_file.parent_path()));
}

TEST_F(io_config_test_f, get_script_runner)
{
    fs::path script_runner = io_cfg.get_script_runner();
    EXPECT_EQ("app_run.R", script_runner.filename().string());
    EXPECT_TRUE(fs::exists(script_runner.parent_path()));
}

TEST_F(io_config_test_f, get_tax_file)
{
    fs::path tax_file = io_cfg.get_tax_file();
    EXPECT_EQ("root_3041.tax", tax_file.filename().string());
    EXPECT_TRUE(fs::exists(tax_file));
}

TEST_F(io_config_test_f, get_work_dir)
{
    fs::path work_dir = io_cfg.get_work_dir();
    EXPECT_EQ(fs::canonical(work_dir), io_cfg.get_work_dir());
    EXPECT_TRUE(fs::exists(work_dir));
}

TEST_F(io_config_test_f, get_library_size)
{
    EXPECT_EQ(245UL, io_cfg.get_library_size());
}


TEST_F(io_config_test_f, get_species_count)
{
    EXPECT_EQ(9051ULL, io_cfg.get_species_count());
}

TEST_F(io_config_test_f, is_species)
{
    EXPECT_EQ(false, io_cfg.is_species(3166));
    EXPECT_EQ(false, io_cfg.is_species(1764312));
    EXPECT_EQ(true, io_cfg.is_species(1615899));
    EXPECT_EQ(true, io_cfg.is_species(1478115));
}

TEST_F(io_config_test_f, get_acc_by_seqNo)
{
    std::vector<std::pair<TSeqNo, Accession>> const accs =
        {{1, "X05806.1"}, {2, "X55877.1"}, {245, "L41348.1"}};
    for (std::pair<TSeqNo, Accession> p : accs)
        EXPECT_EQ(p.second, io_cfg.get_acc_by_seqNo(p.first));
}

TEST_F(io_config_test_f, get_taxid_by_accession)
{
    // check presence of first, random and last rows of root_3041.acc
    std::vector<std::pair<Taxid, std::vector<Accession>>> rows = {
        {3190, {"X69220"}},
        {3074, {"D11346", "DI448914", "DI486630", "DI486668", "HW821071", "M68064", "M68064.1", "M76669", "X55349", "X55349.1"}},
        {3046, {"HW823297", "M23531", "M23531.1", "M84320", "M84320.1", "X15280"}}
    };
    for (std::pair<Taxid, std::vector<Accession>> row : rows)
    {
        for (Accession acc : row.second)
        {
            EXPECT_EQ(row.first, io_cfg.get_taxid_by_accession(acc));
        }
    }
}

TEST_F(io_config_test_f, get_taxid_by_seqNo)
{
    // check presence of first, random and last sequence from FASTA file
    std::vector<std::tuple<Accession, TSeqNo, Taxid>> data = {
        {"X05806", 1, 35845},
        {"X68928.1", 128, 28459},
        {"L41348.1", 245, 3059}
    };
    for (std::tuple<Accession, TSeqNo, Taxid> item : data)
        EXPECT_EQ(std::get<2>(item), io_cfg.get_taxid_by_seqNo(std::get<1>(item)));
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
