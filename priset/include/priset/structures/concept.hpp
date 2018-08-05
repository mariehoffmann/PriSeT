// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================

/*!\file
 * \brief Core concepts.
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 */

#pragma once

#include <priset/taxdb_config.hpp>

namespace priset
{
/*!\addtogroup structures
 * \{
 */
//!\cond
template <typename taxid_type, template <typename TaxTree> >
concept bool taxtree_concept = requires (TaxTree t)
{
    typename TaxTree::taxid_type;
    typename TaxTree::node_type;
    typename TaxTree::rank_type;
    typename TaxTree::name_type;
    // require inner Node class
    { t.build_taxtree(TaxTree::taxid_type, &TaxDBConfig) } -> Bool;
    { t.get_root() } -> &TaxTree::node_type;
    { t.get_children(TaxTree::taxid_type) } -> std::vector<TaxTree::node_type>;
    { t.get_parent(TaxTree::taxid_type) } -> TaxTree::taxid_type;
    { t.get_rank(TaxTree::taxid_type) } -> TaxTree::rank_type;
    { t.get_scientific_name(TaxTree::taxid_type) } -> TaxTree::name_type;

};
//!\endcond

//!\}

}
