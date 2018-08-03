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

/*!\addtogroup combiner
 * \{
*/
//!\cond
template <typename DatabaseConnector>
concept bool database_connector_concept = requires (DatabaseConnector d)
{
    { d.query() } -> &std::String;
    { d.connection_open() } -> bool;
    { d.connection_close() } -> bool;
};
//!\endcond

//!\cond
template <typename AccessionCollector>
concept bool accession_collector_concept = requires (AccessionCollector a)
{
    typename AccessionCollector::taxmap_type;
    { v.get_tax2accs_map() } -> &AccessionCollector::taxmap_type;
    { v.run() } -> bool;
};
//!\endcond
//\}
