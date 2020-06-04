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
// 0: 0 with 0, i.e. length pattern l_min of kmerID1 combined with l_min of kmerID2
// 1: 0 with 1
// x: x/10 with x%10
template<typename TKmerID, typename TKmerLength>
struct CombinePattern
{
private:
    // TODO: use union type for mask when doing SIMD vectorization
    std::bitset<100> data;  //? or union with SIMD 128

public:

    using TOffset = uint8_t;

    // return true if at least one combination bit is set.
    inline bool is_set()
    {
        return data.any();
    }
    // set a kmer combination by its lengths given the maximal length difference
    inline void set(uint64_t const prefix1, uint64_t const prefix2) noexcept
    {
        auto idx = __builtin_clzl(prefix1) * PREFIX_SIZE + __builtin_clzl(prefix2); // in [0:l_max^2[
        data.set(idx);
    }

    // unset bit if length combination doesn't pass a filter anymore
    // To be reset bit is expressed as offset w.r.t. PRIMER_MIN_LEN.
    inline void reset(TOffset const k_offset1, TOffset const k_offset2)
    {
        assert(k_offset1 <= PREFIX_SIZE && k_offset2 <= PREFIX_SIZE);
        data.reset(k_offset1 * PREFIX_SIZE + k_offset2);
    }

    // The number of combinations stored in data.
    uint64_t size()
    {
        //std::cout << "Enter size in cp: popcount_0 = " << __builtin_popcountll(data[0]) << ", popcount_1 = " << __builtin_popcountll(data[1]) << std::endl;
        return __builtin_popcountll(data[0]) + __builtin_popcountll(data[1]);
    }

    // Return all enumerated length combinations translated into kmer length offsets, i.e.
    // the true kmer length can be retrieved by adding PRIMER_MIN_LEN.
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

    constexpr bool operator[](std::size_t pos) const
    {
        return data[pos];
    }

    // Wrapper for bitset::none(), returns true if no bit is set, else false.
    constexpr bool none() const noexcept
    {
        return data.none();
    }

    std::string to_string() const noexcept
    {
        return data.to_string();
    }
};

}   // namespace priset
