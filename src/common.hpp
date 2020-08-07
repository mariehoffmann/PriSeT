
namespace priset
{
// Split prefix and code given a kmerID.
#define split_kmerID(kmerID) std::pair<uint64_t, uint64_t>{(uint64_t)kmerID & PREFIX_SELECTOR, (uint64_t)kmerID & ~PREFIX_SELECTOR}

// The longest encoded k-mer length expressed in 2-bit format, i.e. 16 bp are 32 bits.
#define encoded_length(kmerID) WORD_SIZE - 1 - __builtin_clzl(kmerID & ~PREFIX_SELECTOR)

// The longest encoded k-mer length in prefix mask expressed in 2-bit format, i.e. 16 bp are 32 bits.
#define encoded_length_mask(mask) ((KAPPA_MIN + __builtin_clzl(mask)) << 1)


// Reset length bit in prefix encoded in 2-bit format.
// remains rest of kmerID unmodified, ensure that the prefix bit is set, otherwise
// it will be inverted (0 -> 1)
#define reset_length(kmerID, l) kmerID ^ (1ULL << (WORD_SIZE - 1 - (l >> 1) + KAPPA_MIN))


}  // namespace priset
