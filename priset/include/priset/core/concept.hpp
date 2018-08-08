// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================

/*!\file
 * \brief Core concepts.
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 */

#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <priset/taxtree.hpp>

namespace priset
{
//!\cond
template <typename PrimerConfig>
concept bool primer_config_concept = requires (PrimerConfig p)
{
    typename PrimerConfig::float_type;
    typename PrimerConfig::sequence_type;
    typename PrimerConfig::size_type;
    typename PrimerConfig::taxid_type;

    { p.check_primer_constraints(PrimerConfig::sequence_type) } -> bool;
    { p.check_primerpair_constraints(PrimerConfig::sequence_type, PrimerConfig::sequence_type) } -> bool;
    { p.set_root_taxid(PrimerConfig::taxid_type) } -> bool;
    { p.get_root_taxid() } -> PrimerConfig::taxid_type;
    { p.set_coverage_rate(PrimerConfig::float_type) } -> bool;
    { p.get_coverage_rate() } -> PrimerConfig::float_type;
    { p.set_mismatch_number(PrimerConfig::size_type) } -> bool;
    { p.get_mismatch_number() } -> PrimerConfig::size_type;
    { p.set_distance_range(std::pair<PrimerConfig::size_type, PrimerConfig::size_type>)} -> bool;
    { p.get_distance_range()} -> std::pair<PrimerConfig::size_type, PrimerConfig::size_type>;
};
//!\endcond

//!\cond
template <typename SequenceCollector>
concept bool sequence_collector_concept = requires (SequenceCollector s)
{
    typename SequenceCollector::taxmap_type;
    { s.run() } -> bool;
    { s.get_references() } -> ;
    { s.get_tax2acc_map() } -> SequenceCollector::taxmap_type;
    { s.get_taxtree() } -> TaxTree & ;
};
//!\endcond

//!\cond
template <typename HoReFinder>
concept bool homologous_region_finder_concept = requires (HoReFinder h)
{
    { h.run() } -> bool;
    typename HoReFinder::sequence_type;
    typename HoReFinder::block_type;
    { h.get_centroids() } -> std::unorderd_set<HoReFinder::sequence_type>;
};
//!\endcond

/*!\addtogroup combiner
 * \{
 */
//!\cond
template <typename Combiner>
concept bool combiner_concept = requires (Combiner c)
{
    typename Combiner::blockpair_type;
    { c.combine_blocks() } -> &std::unorderd_set<Combiner::blockpair_type>;
};
//!\endcond

//!\cond
template <typename Resolution>
concept bool resolution_concept = requires (Resolution r)
{
    typename Resolution::blockpair_type;
    typename Resolution::taxid_type;
    { r.get_primers() } -> &Resolution::blockpair_type;
    { r.get_lca_set() } -> &std::unorderd_set<Resolution::taxid_type>;
};
//!\endcond

//!\cond
template <typename Evaluator>
concept bool evaluator_concept = requires (Evaluator e)
{
    typename Evaluator::block_type;
    typename Evaluator::resolution_type;
    typename Evaluator::sequence_type;
    typename Evaluator::taxid_type;
    { e.eval() } -> std::unorderd_set<Evaluator::resolution_type>;
};
//!\endcond

//!\cond
template <typename Visualizer>
concept bool visualizer_concept = requires (Visualizer v)
{
    typename Visualizer::resolution_type;
    typename Evaluator::taxid_type;
    { v.display() } -> bool;
    { v.export() } -> bool;
    { v.close() } -> bool;
    { v.print_pdf() } -> bool;
};
//!\endcond

//!\}

}
