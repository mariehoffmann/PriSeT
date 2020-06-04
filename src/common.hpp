
namespace priset
{
// Split prefix and code given a kmerID.
#define split_kmerID(kmerID) std::pair<uint64_t, uint64_t>{(uint64_t)kmerID & PREFIX_SELECTOR, (uint64_t)kmerID & ~PREFIX_SELECTOR}

// The longest encoded k-mer length expressed in 2 bit format, i.e. 16 bp are 32 bits.
#define encoded_length(kmerID) WORD_SIZE - 1 - __builtin_clzl(kmerID & ~PREFIX_SELECTOR);

}  // namespace priset
