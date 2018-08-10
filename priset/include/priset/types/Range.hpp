// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================

/*!\file
* \brief Range type for storing upper and lower bounds.
* \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
*/

#pragma once

namespace priset
{
template<typename value_type>
struct Range
{
    value_type min;
    value_type max;
};
} // priset
