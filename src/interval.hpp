// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================

/*!\file
* \brief Interval type for storing upper and lower bounds.
* \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
*/

#pragma once

namespace priset
{
template<typename value_type>
struct interval
{
    bool in(value_type const value) const
    {
        assert(min <= max && "max set to a value less than min!");
        return min <= value && value <= max;
    }
    value_type min;
    value_type max;
};
} // priset
