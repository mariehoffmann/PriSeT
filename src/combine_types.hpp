#pragma once

#include "vector"

#include "types.hpp"
#include "utilities.hpp"

namespace priset
{

template<typename uint_type>
std::string bits2str(uint_type i);

// Store enumerated length combinations of two kmers in two 64 bit unsigned integers.
// Enumerations follow lexicographical ordering. Since a kmerID may store up to
// 10 different kmer lengths, we have 100 possible kmer combinations. One bit for
// each kmer combination is reserved in the mask in big endian fashion, s.t. mask
// can be seen as the concatenation of 2x64 bits.
// 0: 0 with 0, i.e. length pattern l_min of kmerID1 combined with l_min of kmerID2
// 1: 0 with 1
// x: x/10 with x%10
template<typename TKmerID, typename TKmerLength>
struct TCombinePattern
{
private:
    // TODO: use union type for mask when doing SIMD vectorization
    std::uint64_t data[2]; // or std::bitset<100> ? or union with SIMD 128

public:

    // return true if at least one combination bit is set.
    inline bool is_set()
    {
        return (data[0] | data[1]);
    }
    // set a kmer combination by its lengths given the maximal length difference
    inline void set(TKmerID const mask1, TKmerID const mask2) noexcept
    {
        auto idx = __builtin_clzl(mask1) * LEN_MASK_SIZE + __builtin_clzl(mask2); // in [0:l_max^2[
    //    std::cout << "computed index for bit set in cp: " << idx << std::endl;
        data[idx >> 6] += 1 << (WORD_SIZE - 1 - (idx % WORD_SIZE));

        //std::cout << "DEBUG: set data at position " << (idx >> 6) << std::endl;
    }

    // unset bit if length combination doesn't pass a filter anymore
    inline void unset(TKmerLength const k1, TKmerLength const k2, TKmerLength const k_min) noexcept
    {
        uint16_t const l1 = k1 - k_min;
        uint16_t const l2 = k2 - k_min;
        data[(l1 * LEN_MASK_SIZE + l2) >> 6] -= 1 << ((WORD_SIZE - 1) - ((l1 * LEN_MASK_SIZE + l2) % WORD_SIZE));
    }

    // return all enumerated length combinations translated into kmer lengths
    void get_combinations(std::vector<std::pair<TKmerLength, TKmerLength>> & combinations)
    {
        combinations.clear();
        for (uint8_t i = 0; i < 2; ++i)
        {
            uint64_t data_part = data[i];
            while (data_part)
            {
                auto c = i * WORD_SIZE + __builtin_clzl(data_part); // combination index
                if (64 == c)
                    break;
                std::pair<TKmerLength, TKmerLength> pair{c / LEN_MASK_SIZE + PRIMER_MIN_LEN, (c % LEN_MASK_SIZE) + PRIMER_MIN_LEN};
                combinations.push_back(pair);
                data_part &= ((1 << (WORD_SIZE - c - 1)) - 1); // delete leading bit
            }
        }
    }
};

/*
 * A pair of encoded kmers (kmer IDs) given by their indices in the kmer ID
 * container associated to a reference bit vector. Since kmer IDs may encode up
 * to |[primer_max_length : primer_min_length]| kmers, the combine pattern stores
 * which kmers are combined. See encoding scheme here: TCombinePattern.
 */
template<typename TCombinePattern>
struct TPair
{
    TPair() = default;
    TPair(uint64_t reference_, uint64_t r_fwd_, uint64_t r_rev_, TCombinePattern cp_) :
        reference(reference_), r_fwd(r_fwd_), r_rev(r_rev_), cp(cp_) {}
    //std::tuple<uint64_t, uint64_t, TCombinePattern<TKmerID, TKmerLength>> TPair;
    ~TPair() = default;

    // Reference identifier (equivalent to position in corpus).
    uint64_t reference;
    // Rank of forward kmer.
    uint64_t r_fwd;
    // Rank of reverse kmer.
    uint64_t r_rev;
    // Length combinations of kmers as bit mask.
    TCombinePattern cp;
};

template<typename TPair>
using TPairList = std::vector<TPair>;
} // priset
