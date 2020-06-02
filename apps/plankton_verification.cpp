// #include <array>
// #include <chrono>
// #include <cstdlib>
// #include <ctime>
// #include <iostream>
// #include <experimental/filesystem>
// #include <fstream>
// #include <numeric>
// #include <regex>
// #include <sstream>
// #include <string>
// #include <sys/wait.h>
// #include <unistd.h>
// #include <vector>
//
// #include "../src/argument_parser.hpp"
// #include "../src/filter.hpp"
// #include "../src/fm.hpp"
// #include "../src/types/all.hpp"
// #include "../src/utilities.hpp"
//
// namespace fs = std::experimental::filesystem;
//
// // g++ ../PriSeT/apps/plankton_verification.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -Wno-unknown-pragmas -lstdc++fs -DNDEBUG -O3 -I~ -L~/lib -ldivsufsort -o verif
// // ./verif $taxid $lib_dir $work_dir
//
// using namespace priset;
//
// std::string dna_decoder(uint64_t, uint64_t const);
//
//
// struct TPrimerKey
// {
//     uint64_t fwd;
//     uint64_t rev;
//     TPrimerKey(uint64_t fwd_, uint64_t rev_) : fwd(fwd_), rev(rev_) {}
//     bool operator==(TPrimerKey const & lhs) const
//     {
//         return (fwd == lhs.fwd && rev == lhs.rev);
//     }
// };
//
// // hash defining in global namespace is ill-formed, needs to be declared in same namespace like std::hash!
// //namespace std
// //{ // doesn't work with struct hash, overwriting the one in std::hash
// struct hash_pp
// {
//     template<typename TPrimerKey>
//     std::size_t operator()(TPrimerKey const & key) const
//     {
//         return std::hash<std::string>()(std::to_string(key.fwd) + std::to_string(key.rev));
//     }
// };
//
// // Verify that known primers pass the chemical filter
// void chemical_filter_test(std::unordered_map<TPrimerKey, std::string, hash_pp> & pairs_known)
// {
//     TKmerID fwd, rev;
//     for (auto p : pairs_known)
//     {
//         fwd = p.first.fwd;
//         chemical_filter_single_pass(fwd);
//         rev = p.first.rev;
//         chemical_filter_single_pass(rev);
//         std::cout << p.second << ": ";
//         if ((PREFIX_SELECTOR & fwd) && (PREFIX_SELECTOR & rev))
//             std::cout << " OK\n";
//         else if (!(PREFIX_SELECTOR & fwd) && (PREFIX_SELECTOR & rev))
//             std::cout << " FWD DID NOT PASS FILTER!\n";
//         else if ((PREFIX_SELECTOR & fwd) && !(PREFIX_SELECTOR & rev))
//             std::cout << " REV DID NOT PASS FILTER!\n";
//         else
//             std::cout << " FWD + REV DID NOT PASS FILTER!\n";
//     }
// }
//
// // Note: reverse primers are noted in 5' to 3' direction for the forward sting
// // (no reverse complement like needed for real PCR!), because we search on same
// // string directly in combiner. So, we have to reverse complement a published primer
// // eg. noted as GCTTAATTTGACTCAACACGGG
// // For later display the reverse sequence would be transformed.
// template<typename TPrimerKey>
// void load_primers_known(std::unordered_map<TPrimerKey, std::string, hash_pp> & pairs_known)
// {
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GGACAGAAAGACCCTATG"), dna_encoder_with_lbit("GGATAACAGGCTGATCT")}] = "23S 0 0";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GGACAAAAAGACCCTATG"), dna_encoder_with_lbit("GGATAACAGGCTGATCT")}] = "23S 1 0";
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("TGTTGCAGTTAAAAAGCTCGT"), dna_encoder_with_lbit("GTTGGGGGTGCTAGTATTCA")}] = "CERC";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("TGGCCTATCTTGTTGGTCTGT"), dna_encoder_with_lbit("GTTGCCTTGTCAGGTTGATTC")}] = "CHLORO";
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATTCC"), dna_encoder_with_lbit("TAGGGGATCAAAGACGATCAGA")}] = "CV 0 0";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATTCC"), dna_encoder_with_lbit("TAGGGGATCAAAGACGACCAGA")}] = "CV 0 1";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATTCC"), dna_encoder_with_lbit("TAGGGGATCAAAGACAATCAGA")}] = "CV 0 2";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATTCC"), dna_encoder_with_lbit("TAGGGGATCAAAGACAACCAGA")}] = "CV 0 3";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATTCC"), dna_encoder_with_lbit("TAAGGGATCAAAGACGATCAGA")}] = "CV 0 4";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATTCC"), dna_encoder_with_lbit("TAAGGGATCAAAGACGACCAGA")}] = "CV 0 5";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATTCC"), dna_encoder_with_lbit("TAAGGGATCAAAGACAATCAGA")}] = "CV 0 6";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATTCC"), dna_encoder_with_lbit("TAAGGGATCAAAGACAACCAGA")}] = "CV 0 7";
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATACC"), dna_encoder_with_lbit("TAAGGGATCAAAGACAACCAGA")}] = "CV 1 0";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATACC"), dna_encoder_with_lbit("TAGGGGATCAAAGACGACCAGA")}] = "CV 1 1";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATACC"), dna_encoder_with_lbit("TAGGGGATCAAAGACAATCAGA")}] = "CV 1 2";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATACC"), dna_encoder_with_lbit("TAGGGGATCAAAGACAACCAGA")}] = "CV 1 3";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATACC"), dna_encoder_with_lbit("TAAGGGATCAAAGACGATCAGA")}] = "CV 1 4";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATACC"), dna_encoder_with_lbit("TAAGGGATCAAAGACGACCAGA")}] = "CV 1 5";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATACC"), dna_encoder_with_lbit("TAAGGGATCAAAGACAATCAGA")}] = "CV 1 6";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATACC"), dna_encoder_with_lbit("TAAGGGATCAAAGACAACCAGA")}] = "CV 1 7";
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATTCC"), dna_encoder_with_lbit("TAGGGGATCAAAGACGATCAGA")}] = "CV 2 0";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATTCC"), dna_encoder_with_lbit("TAGGGGATCAAAGACGACCAGA")}] = "CV 2 1";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATTCC"), dna_encoder_with_lbit("TAGGGGATCAAAGACAATCAGA")}] = "CV 2 2";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATTCC"), dna_encoder_with_lbit("TAGGGGATCAAAGACAATCAGA")}] = "CV 2 3";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATTCC"), dna_encoder_with_lbit("TAAGGGATCAAAGACGATCAGA")}] = "CV 2 4";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATTCC"), dna_encoder_with_lbit("TAAGGGATCAAAGACGACCAGA")}] = "CV 2 5";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATTCC"), dna_encoder_with_lbit("TAAGGGATCAAAGACAATCAGA")}] = "CV 2 6";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATTCC"), dna_encoder_with_lbit("TAAGGGATCAAAGACAACCAGA")}] = "CV 2 7";
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATACC"), dna_encoder_with_lbit("TAGGGGATCAAAGACGATCAGA")}] = "CV 3 0";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATACC"), dna_encoder_with_lbit("TAGGGGATCAAAGACGACCAGA")}] = "CV 3 1";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATACC"), dna_encoder_with_lbit("TAGGGGATCAAAGACAATCAGA")}] = "CV 3 2";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATACC"), dna_encoder_with_lbit("TAGGGGATCAAAGACAACCAGA")}] = "CV 3 3";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATACC"), dna_encoder_with_lbit("TAAGGGATCAAAGACGATCAGA")}] = "CV 3 4";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATACC"), dna_encoder_with_lbit("TAAGGGATCAAAGACGACCAGA")}] = "CV 3 5";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATACC"), dna_encoder_with_lbit("TAAGGGATCAAAGACAATCAGA")}] = "CV 3 6";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATACC"), dna_encoder_with_lbit("TAAGGGATCAAAGACAACCAGA")}] = "CV 3 7";
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("ATTCCAGCTCCAATAGCG"), dna_encoder_with_lbit("GATTAGATACCATCGTAGTC")}] = "DIAZ";
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GGGGACGGGTGAAATAGGATG"), dna_encoder_with_lbit("GGAGTCTGCGGCTCAATTTG")}] = "DIM 0 0";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GGGGACAGGTGAAATAGGATG"), dna_encoder_with_lbit("GGAGTCTGCGGCTCAATTTG")}] = "DIM 1 0";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("AGGGACGGGTGAAATAGGATG"), dna_encoder_with_lbit("GGAGTCTGCGGCTCAATTTG")}] = "DIM 2 0";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("AGGGACAGGTGAAATAGGATG"), dna_encoder_with_lbit("GGAGTCTGCGGCTCAATTTG")}] = "DIM 3 0";
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCGGTAATTCCAGCTCCAATAG"), dna_encoder_with_lbit("TATTCGTATTCCATTGTCAGAG")}] = "DIV4";
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CAGCAGCCGCGGTAATTCC"), dna_encoder_with_lbit("GCTTAATTTGACTCAACACGGG")}] = "EUK14";
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCTGCGGTAATTCC"), dna_encoder_with_lbit("TTGATCAAGAACGAAAGT")}] = "EUK15 0 0";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCTGCGGTAATTCC"), dna_encoder_with_lbit("TCGATCAAGAACGAAAGT")}] = "EUK15 0 1";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCTGCGGTAATTCC"), dna_encoder_with_lbit("TTAATCAAGAACGAAAGT")}] = "EUK15 0 2";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCTGCGGTAATTCC"), dna_encoder_with_lbit("TCAATCAAGAACGAAAGT")}] = "EUK15 0 3";
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATTCC"), dna_encoder_with_lbit("TTGATCAAGAACGAAAGT")}] = "EUK15 1 0";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATTCC"), dna_encoder_with_lbit("TCGATCAAGAACGAAAGT")}] = "EUK15 2 1";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATTCC"), dna_encoder_with_lbit("TTAATCAAGAACGAAAGT")}] = "EUK15 3 2";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCAGCCGCGGTAATTCC"), dna_encoder_with_lbit("TCAATCAAGAACGAAAGT")}] = "EUK15 4 3";
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCTGCGGTAATTCC"), dna_encoder_with_lbit("TTGATCAAGAACGAAAGT")}] = "EUK15 2 0";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCTGCGGTAATTCC"), dna_encoder_with_lbit("TCGATCAAGAACGAAAGT")}] = "EUK15 2 1";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCTGCGGTAATTCC"), dna_encoder_with_lbit("TTAATCAAGAACGAAAGT")}] = "EUK15 2 2";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCTGCGGTAATTCC"), dna_encoder_with_lbit("TCAATCAAGAACGAAAGT")}] = "EUK15 2 3";
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATTCC"), dna_encoder_with_lbit("TTGATCAAGAACGAAAGT")}] = "EUK15 3 0";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATTCC"), dna_encoder_with_lbit("TCGATCAAGAACGAAAGT")}] = "EUK15 3 1";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATTCC"), dna_encoder_with_lbit("TTAATCAAGAACGAAAGT")}] = "EUK15 3 2";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CCAGCACCCGCGGTAATTCC"), dna_encoder_with_lbit("TCAATCAAGAACGAAAGT")}] = "EUK15 3 3";
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCCGCGGTAATTCCAGCTC"), dna_encoder_with_lbit("TTAATCAAGGGCGAAGG")}] = "EUKA 0 0";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCCGCGGTAATTCCAGCTC"), dna_encoder_with_lbit("TCAATCAAGGGCGAAGG")}] = "EUKA 0 1";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCCGCGGTAATTCCAGCTC"), dna_encoder_with_lbit("TTAATCAAGGGCGAAAG")}] = "EUKA 0 2";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCCGCGGTAATTCCAGCTC"), dna_encoder_with_lbit("TCAATCAAGGGCGAAAG")}] = "EUKA 0 3";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCCGCGGTAATTCCAGCTC"), dna_encoder_with_lbit("TTAATCAAGGACGAAGG")}] = "EUKA 0 4";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCCGCGGTAATTCCAGCTC"), dna_encoder_with_lbit("TCAATCAAGGACGAAGG")}] = "EUKA 0 5";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCCGCGGTAATTCCAGCTC"), dna_encoder_with_lbit("TTAATCAAGGACGAAAG")}] = "EUKA 0 6";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCCGCGGTAATTCCAGCTC"), dna_encoder_with_lbit("TCAATCAAGGACGAAAG")}] = "EUKA 0 7";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCCGCGGTAATTCCAGCTC"), dna_encoder_with_lbit("TTAATCAAGAGCGAAGG")}] = "EUKA 0 8";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCCGCGGTAATTCCAGCTC"), dna_encoder_with_lbit("TCAATCAAGAGCGAAGG")}] = "EUKA 0 9";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCCGCGGTAATTCCAGCTC"), dna_encoder_with_lbit("TTAATCAAGAGCGAAAG")}] = "EUKA 0 10";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCCGCGGTAATTCCAGCTC"), dna_encoder_with_lbit("TCAATCAAGAGCGAAAG")}] = "EUKA 0 11";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCCGCGGTAATTCCAGCTC"), dna_encoder_with_lbit("TTAATCAAGAACGAAGG")}] = "EUKA 0 12";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCCGCGGTAATTCCAGCTC"), dna_encoder_with_lbit("TCAATCAAGAACGAAGG")}] = "EUKA 0 13";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCCGCGGTAATTCCAGCTC"), dna_encoder_with_lbit("TTAATCAAGAACGAAAG")}] = "EUKA 0 14";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCCGCGGTAATTCCAGCTC"), dna_encoder_with_lbit("TCAATCAAGAACGAAAG")}] = "EUKA 0 15";
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("AACCTGGTTGATCCTGCCAGT"), dna_encoder_with_lbit("GTAGGTGAACCTGCAGAAGGATCA")}] = "FAD4";
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("GCTTGTCTCAAAGATTAAGCC"), dna_encoder_with_lbit("TCCAAGGAAGGCAGCAGGC")}] = "nSSU";
//
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CGCGGTAATTCCAGCTTC"), dna_encoder_with_lbit("GTCAGAGGTGAAATTCTTGGAT")}] = "SSU 0 0";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CGCGGTAATTCCAGCTTC"), dna_encoder_with_lbit("GTCAGAGGTGAAATTCTTGAAT")}] = "SSU 0 1";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CGCGGTAATTCCAGCTCC"), dna_encoder_with_lbit("GTCAGAGGTGAAATTCTTGGAT")}] = "SSU 1 0";
//     pairs_known[TPrimerKey{dna_encoder_with_lbit("CGCGGTAATTCCAGCTCC"), dna_encoder_with_lbit("GTCAGAGGTGAAATTCTTGAAT")}] = "SSU 1 1";
// }
//
//
// // same as combine, but accepts primer pair set to be verified
// template<typename PairList, typename TPrimerKey>
// void combine2(TReferences const & references, TKmerIDs const & kmerIDs, PairList & pairs, std::unordered_map<TPrimerKey, std::string, hash_pp> & pairs_known, std::unordered_set<std::string> & verified)
// {
//     pairs.clear();
//     for (uint64_t seqNo = 0; seqNo < references.size(); ++seqNo)
//     {
//         std::cout << "STATUS: current reference i = " << seqNo + 1 << "/" << references.size() << std::endl;
//         sdsl::bit_vector reference;
//         sdsl::util::assign(reference, references.at(seqNo));
//         sdsl::rank_support_v5<1, 1> r1s(&references.at(seqNo)); // replace after bugfix with
//         sdsl::select_support_mcl<1> s1s(&reference);
//
//         for (uint64_t r_fwd = 1; r_fwd < r1s.rank(reference.size()); ++r_fwd)
//         {
//             uint64_t idx_fwd = s1s.select(r_fwd);  // text position of r-th k-mer
//             TKmerID const kmerID_fwd = kmerIDs[seqNo][r_fwd - 1];
//             if (!(kmerID_fwd >> CODE_SIZE))
//             {
//                 std::cout << "ERROR: k length pattern is zero\n";
//                 exit(-1);
//             }
//
//             uint64_t w_begin = idx_fwd + PRIMER_MIN_LEN + TRANSCRIPT_MIN_LEN;
//             uint64_t w_end = std::min(reference.size(), idx_fwd + PRIMER_MAX_LEN + TRANSCRIPT_MAX_LEN + 1); // + 1: rank excludes upper bound
//
//             for (uint64_t r_rev = r1s.rank(w_begin) + 1; r_rev <= r1s.rank(w_end); ++r_rev)
//             {
//                 CombinePattern<TKmerID, TKmerLength> cp;
//                 uint64_t mask_fwd = ONE_LSHIFT_63;
//                 while ((((mask_fwd - 1) << 1) & kmerID_fwd) >> 54)
//                 {
//                     if ((mask_fwd & kmerID_fwd)) //&& filter_CG_clamp(kmerID_fwd, '+', mask_fwd) && filter_WWW_tail(kmerID_fwd, '+', mask_fwd))
//                     {
//                         TKmerID const kmerID_rev = kmerIDs.at(seqNo).at(r_rev - 1);
//                         uint64_t mask_rev = ONE_LSHIFT_63;
//                         while ((((mask_rev - 1) << 1) & kmerID_rev) >> 54)
//                         {
//                             if ((mask_rev & kmerID_rev)) // && filter_CG_clamp(kmerID_rev, '-') && filter_WWW_tail(kmerID_fwd, '+', mask_fwd))
//                             {
//                                 if (dTm(kmerID_fwd, mask_fwd, kmerID_rev, mask_rev) <= PRIMER_DTM)
//                                 {
//                                     TPrimerKey key{get_code(kmerID_fwd, mask_fwd) | mask_fwd, get_code(kmerID_rev, mask_rev) | mask_rev};
//                                     if (pairs_known.find(key) != pairs_known.end())
//                                     {
//                                         std::cout << "INFO: primer pair <" << pairs_known[key] << "> found for refID = " << seqNo << " at " << s1s.select(r_fwd) << " and " << s1s.select(r_rev) << std::endl;
//                                         verified.insert(pairs_known[key]);
//                                     }
//                                 }
//                             }
//                             mask_rev >>= 1; // does not affect search window, since starting position is fixed
//                         } // length mask_rev
//                     }
//                     mask_fwd >>= 1;
//                 } // length mask_fwd
//
//                 if (cp.is_set())
//                     pairs.push_back(Pair<CombinePattern<TKmerID, TKmerLength>>{seqNo, r_fwd, r_rev, cp});
//             } // kmerID rev
//         } // kmerID fwd
//     }
// }
//
//
// int main(int argc, char ** argv)
// {
//     // test chemically filter on primers independently
//     std::unordered_map<TPrimerKey, std::string, hash_pp> pairs_known;
//     load_primers_known<TPrimerKey>(pairs_known);
//     chemical_filter_test(pairs_known);
//
//     // exit(0);
//     if (argc != 4)
//     {
//         std::cout << "Give taxid, and paths to lib and work dirs.\n";
//         exit(-1);
//     }
//     // setup su{argv[2], argv[3]};
//     std::string taxid = argv[1];
//     unsigned const priset_argc = 6;
//     char * const priset_argv[priset_argc] = {"priset", "-l", argv[2], "-w", argv[3], "-s"};
//     for (unsigned i = 0; i < priset_argc; ++i) std::cout << priset_argv[i] << " ";
//     std::cout << std::endl;
//
//     // collect number of kmers or kmer pairs left after relevant processing steps
//     TKmerCounts kmerCounts{0, 0, 0, 0};
//
//     // set path prefixes for library files
//     IOConfig io_cfg{};
//
//     // get instance to primer sequence settings
//     PrimerConfig primer_cfg{};
//
//     // parse options and init io and primer configurators
//     options opt(priset_argc, priset_argv, primer_cfg, io_cfg);
//
//     std::chrono::time_point<std::chrono::system_clock> start, finish;
//
//     TKLocations locations;
//     TDirectoryInformation directoryInformation;
//     TSequenceNames sequenceNames;
//     TSequenceLengths sequenceLengths;
//
//     // compute k-mer mappings
//     fm_map(io_cfg, primer_cfg, locations);
//     auto kmer_cnt = std::accumulate(locations.begin(), locations.end(), 0, [](unsigned ctr, auto & location){return ctr + location.second.first.size();});
//     std::cout << "INFO: kmers init = " << kmer_cnt << std::endl;
//     if (!kmer_cnt)
//         exit(0);
//     TReferences references;
//     TKmerIDs kmerIDs;
//     TSeqNoMap seqNoMap;
//     transform_and_filter(io_cfg, primer_cfg, locations, references, seqNoMap, kmerIDs, &kmerCounts);
//     std::cout << "INFO: kmers after filter1 & transform = " << get_num_kmers(kmerIDs) << std::endl;
//
//     // TODO: delete locations
//     using PairList = PairList<Pair<CombinePattern<TKmerID, TKmerLength>>>;
//     PairList pairs;
//
//     std::unordered_set<std::string> verified;
//
//     combine2<PairList, TPrimerKey>(references, kmerIDs, pairs, pairs_known, verified);
//     std::cout << "INFO: pairs after combiner = " << get_num_pairs<PairList>(pairs) << std::endl;
//
//     std::cout << "Verified primers for current clade: \n";
//     if (!verified.size())
//         std::cout << 0 << std::endl;
//     else
//     {
//         for (std::string pp : verified)
//             std::cout << pp << std::endl;
//     }
//
//     return 0;
// }
