// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================

/*!\file
 * \brief Primer Settings.
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 */

#include <unordered_map>

namespace priset
{
//!\brief DNA codes as enums.
enum dna {A, C, G, T, B, CGT, D, AGT, H, ACT, K, GT, M, AC, N, ACGT, R, AG, S, \
    CG, V, ACG, W, AT, Y, CT};

//!\brief String representation of DNA enums.
std::array<std::string, 26> dna2str = {"A", "C", "G", "T", "B", "CGT", "D", \
    "AGT", "H", "ACT", "K", "GT", "M", "AC", "N", "ACGT", "R", "AG", "S", "CG", \
    "V", "ACG", "W", "AT", "Y", "CT"};

//!\brief Conversion table from string representation to dna enum.
std::unordered_map<std::string, dna> dna2str = {{"A", dna::A}, {"C", dna::C}, \
        {"G", dna::G}, {"T", dna::T}, {"B", dna::B}, {"CGT", dna::CGT}, \
        {"D", dna::D}, {"AGT", dna::AGT}, {"H", dna::H}, {"ACT", dna::ACT}, \
        {"K", dna::K}, {"GT", dna::GT}, {"M", dna::M}, {"AC", dna::AC}, \
        {"N", dna::N}, {"ACGT", dna::ACGT}, {"R", dna::}, {"AG", dna::}, \
        {"S", dna::}, {"CG", dna::CG}, {"V", dna::V}, {"ACG", dna::ACG}, \
        {"W", dna::W}, {"AT", dna::AT}, {"Y", dna::Y}, {"CT", dna::CT};


//!\brief bidirectional conversion table for one letter encodings.
// Attention: call always with hashed input, e.g. dna_decode[std::hash<sequence_type>("D")]
static constexpr std::array<priset::dna, 26> dna_decode
{
    [] () constexpr
    {
        std::array<priset::dna, 26> bimap{};
        bimap[priset::dna::A] = priset::dna::A;
        bimap[priset::dna::C] = priset::dna::C;
        bimap[priset::dna::G] = priset::dna::G;
        bimap[priset::dna::T] = priset::dna::T;
        bimap[priset::dna::B] = priset::dna::CGT;
        bimap[priset::dna::CGT] = priset::dna::B;
        bimap[priset::dna::D] = priset::dna::AGT;
        bimap[priset::dna::AGT] = priset::dna::D;
        bimap[priset::dna::H] = priset::dna::ACT;
        bimap[priset::dna::ACT] = priset::dna::H;
        bimap[priset::dna::K] = priset::dna::GT;
        bimap[priset::dna::GT] = priset::dna::K;
        bimap[priset::dna::M] = priset::dna::AC;
        bimap[priset::dna::AC] = priset::dna::M;
        bimap[priset::dna::N] = priset::dna::ACGT;
        bimap[priset::dna::ACGT] = priset::dna::N;
        bimap[priset::dna::R] = priset::dna::AG;
        bimap[priset::dna::AG] = priset::dna::R;
        bimap[priset::dna::S] = priset::dna::CG;
        bimap[priset::dna::CG] = priset::dna::S;
        bimap[priset::dna::V] = priset::dna::ACG;
        bimap[priset::dna::ACG] = priset::dna::V;
        bimap[priset::dna::W] = priset::dna::AT;
        bimap[priset::dna::AT] = priset::dna::W;
        bimap[priset::dna::Y] = priset::dna::CT;
        bimap[priset::dna::CT] = priset::dna::Y;
        return bimap;
    }()
};
}
