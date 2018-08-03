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
