// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================

/*!\file
 * \brief Core concepts.
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 */

#pragma once

#include <vector>

namespace priset
{

/*!\addtogroup database
 * \{
*/

//!\cond
template <typename RefDBConfig>
concept bool refdb_config_concept = requires (RefDBConfig c)
{
    typename RefDBConfig::path_type;
    typename RefDBConfig::filename_type;
    { c.set_basepath(RefDBConfig::path_type) } -> bool;
    { c.get_basepath() } -> RefDBConfig::path_type;
    { c.set_src_files(&std::vector<RefDBConfig::filename_type>) } -> bool;
    { c.get_src_files() } -> &std::vector<RefDBConfig::filename_type>;
};
//!\endcond

//!\cond
template <typename RefDBConnector>
concept bool refdb_connector = requires (RefDBConnector c)
{
    typename RefDBConnector::accession_type;
    typename RefDBConnector::sequence_type;
    { c.get_reference(RefDBConnector::accession_type) } -> RefDBConnector::sequence_type;
};
//!\endcond

//!\cond
template <typename DatabaseConnector>
concept bool database_connector_concept = requires (DatabaseConnector d)
{
    { d.query() } -> std::String;
    { d.connection_open() } -> bool;
    { d.connection_close() } -> bool;
};
//!\endcond

//!\cond
template <typename AccessionCollector>
concept bool accession_collector_concept = requires (AccessionCollector a)
{
    typename AccessionCollector::taxmap_type;
    { a.get_tax2accs_map() } -> &AccessionCollector::taxmap_type;
    { a.run() } -> bool;
};
//!\endcond
//\}
