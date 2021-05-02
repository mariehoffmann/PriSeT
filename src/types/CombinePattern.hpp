// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <ozymandiaz147 AT gmail.com>
//          Manual: https://github.com/mariehoffmann/PriSeT

#pragma once

#include <bitset>
#include <string>
#include <vector>

#include "PrimerConfig.hpp"

namespace priset
{

// Store enumerated length combinations of two kmers in two 64 bit unsigned integers.
// Enumerations follow lexicographical ordering. Since a kmerID may store up to
// 10 different kmer lengths, we have 100 possible kmer combinations. One bit for
// each kmer combination is reserved in the mask in big endian fashion, s.t. mask
// can be seen as the concatenation of 2x64 bits.
struct CombinePattern
{
private:
    // Stores which lengths of two KMerIDs are combined. At most 10 lengths are
    // encoded in a single KMerID, resulting in 10^2 possible sequence combinations.
    // A set bit at position i corresponds to the length combination
    //              l1 = i/10 + KAPPA_MIN
    //              l2 = (i % 10) + KAPPA_MIN
    std::bitset<100> data;

public:

    // The offset type for data.
    using TOffset = uint8_t;

    // return true if at least one combination bit is set.
    inline bool is_set()
    {
        return data.any();
    }

    // Set a k-mer combination given the length masks in TKMerID prefix format.
    inline void set(uint64_t const prefix1, uint64_t const prefix2) noexcept
    {
        assert(__builtin_popcountll(prefix1) == 1 && __builtin_popcountll(prefix2) == 1);
        auto idx = __builtin_clzl(prefix1) * PREFIX_SIZE + __builtin_clzl(prefix2);
        data.set(idx);
    }

    // Unset bit for a specific length combination given both TKMerID prefixes.
    inline void reset(uint64_t const prefix1, uint64_t const prefix2)
    {
        assert(__builtin_popcountll(prefix1) == 1 && __builtin_popcountll(prefix2) == 1);
        auto idx = __builtin_clzl(prefix1) * PREFIX_SIZE + __builtin_clzl(prefix2);
        data.reset(idx);
    }

    // Get the number of combinations stored in data.
    uint64_t size()
    {
        return data.count();
    }

    // Return all enumerated length combinations translated into kmer length offsets.
    //The true kmer length can be retrieved by adding KAPPA_MIN.
    void get_combinations(std::vector<std::pair<TOffset, TOffset>> & combinations)
    {
        combinations.clear();
        for (uint8_t i = 0; i < PREFIX_SIZE * PREFIX_SIZE; ++i)
        {
            if (data[i])
            {
                std::pair<TOffset, TOffset> pair{i / PREFIX_SIZE, (i % PREFIX_SIZE)};
                combinations.push_back(pair);
            }
        }
    }

    // Return set combination bit for
    constexpr bool operator[](size_t pos) const
    {
        return data[pos];
    }

    // Wrapper for bitset::none(), returns true if no bit is set, else false.
    bool none() const noexcept
    {
        return data.none();
    }

    std::string to_string() const noexcept
    {
        return data.to_string();
    }
};

}   // namespace priset
