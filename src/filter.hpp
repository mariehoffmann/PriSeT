#pragma once

namespace priset
{
/*
    using TLocations = std::map<seqan::Pair<TSeqNo, TSeqPos>,
             std::pair<std::vector<seqan::Pair<TSeqNo, TSeqPos> >,
                       std::vector<seqan::Pair<TSeqNo, TSeqPos> > > >;
                       */

// filter k-mers by frequency and chemical properties
template<typename io_config, typename primer_config, typename TLocations>
void filter(io_config & io_cfg, primer_config & primer_cfg, TLocations & locations)
{

}

}  // namespace priset
