#pragma once

#include <numeric>
#include <vector>

#include "../common.hpp"
#include "types/all.hpp"

namespace priset
{

// https://www.linuxquestions.org/questions/programming-9/can-i-predefine-constant-array-in-c-823350/
/*
 * Bit map for (A|T)^3 patterns:
 * 0b000000 = 0b000011 = 0b001100 = 0b001111 = 0b110000 = 0b110011 = 0b111100 = 0b111111 = 1
 */
inline uint8_t const* WWW_tails()
{

    static uint8_t const tails[64] = {  1, 0, 0, 1, 0, 0, 0, 0,    // 0
                                        0, 0, 0, 0, 1, 0, 0, 1,    // 8
                                        0, 0, 0, 0, 0, 0, 0, 0,    // 16
                                        0, 0, 0, 0, 0, 0, 0, 0,    // 24
                                        0, 0, 0, 0, 0, 0, 0, 0,    // 32
                                        0, 0, 0, 0, 0, 0, 0, 0,    // 40
                                        1, 0, 0, 1, 0, 0, 0, 0,    // 48
                                        0, 0, 0, 0, 1, 0, 0, 1};    // 56

    /* not implemented yet in gnu gcc compiler mac distribution: designation
    static uint8_t const tails[64] =
    {
        [0b000000] = 1, [0b000011] = 1, [0b001100] = 1, [0b001111] = 1,
        [0b110000] = 1, [0b110011] = 1, [0b111100] = 1, [0b111111] = 1
    };
    */
    return tails;
}

/*
 * Filter for tails with a (A|T)^3 pattern.
 */
extern inline bool filter_WWW_tail(TKmerID const kmerID, char const sense,
    uint64_t const mask = 0)
{
    if (sense != '+' && sense != '-')
        throw std::invalid_argument("Expected sense '-' or '+'!");
    auto [prefix, code] = split_kmerID(kmerID);
    uint64_t encoded_len = WORD_SIZE - __builtin_clzll(code) - 1;
    if (sense == '-') // reverse complement!
    {
        code >>= encoded_len - 6; // shift right to have prefix of length 10
        return !WWW_tails()[code & 0b111111];
    }
    if (mask)
    {
        uint64_t target_len = (KAPPA_MIN + __builtin_clzl(mask)) << 1;
        code >>= encoded_len - target_len; // delete prefix and trim to target length
    }
    return !WWW_tails()[code & 0b111111];
}

}
