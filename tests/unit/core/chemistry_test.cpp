// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================

#include <gtest/gtest.h>

#include <priset/core/chemistry.hpp>
#include <priset/types/dna.hpp>

using namespace priset;

TEST(chemistry_test, primer_melt_wallace)
{
    using float_type = float;
    using sequence_type = std::vector<priset::dna>;
    float_type EPS = .05;
    sequence_type sequence(16);
    // AAAACCCCGGGGTTTT
    // 0123456789012345
    std::fill(sequence.begin(), sequence.begin()+4, dna::A);
    std::fill(sequence.begin()+4, sequence.begin()+8, dna::C);
    std::fill(sequence.begin()+8, sequence.begin()+12, dna::G);
    std::fill(sequence.begin()+12, sequence.end(), dna::T);
    // 2*cnt_AT + 4*cnt_CG;
    // primer = AAAA
    float_type Tm = chemistry::primer_melt_wallace<sequence_type, float_type>(primer, primer.begin(), primer.begin()+4);
    EXPECT_TRUE(Tm >= 8. - EPS && Tm <= 8. + EPS);

    // primer = AACC
    Tm = chemistry::primer_melt_wallace<sequence_type, float_type>(primer, primer.begin()+2, primer.begin()+6);
    EXPECT_TRUE(Tm >= 10. - EPS && Tm <= 10. + EPS);

    // primer = CCGG
    Tm = chemistry::primer_melt_wallace<sequence_type, float_type>(primer, primer.begin()+6, primer.begin()+10);
    EXPECT_TRUE(Tm >= 16. - EPS && Tm <= 16. + EPS);

    // primer = GTTTT
    Tm = chemistry::primer_melt_wallace<sequence_type, float_type>(primer, primer.begin()+11, primer.end());
    EXPECT_TRUE(Tm >= 12. - EPS && Tm <= 12. + EPS);

    // primer = ''
    Tm = chemistry::primer_melt_wallace<sequence_type, float_type>(primer, primer.begin(), primer.begin());
    EXPECT_TRUE(Tm >= 0. - EPS && Tm <= 0. + EPS);

    // primer = AAAACCCCGGGGTTTT
    Tm = chemistry::primer_melt_wallace<sequence_type, float_type>(primer, primer.begin(), primer.begin());
    EXPECT_TRUE(Tm >= 48. - EPS && Tm <= 48. + EPS);
}
