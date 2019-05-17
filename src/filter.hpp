#pragma once

namespace priset
{
/*
    using TLocations = std::map<seqan::Pair<TSeqNo, TSeqPos>,
             std::pair<std::vector<seqan::Pair<TSeqNo, TSeqPos> >,
                       std::vector<seqan::Pair<TSeqNo, TSeqPos> > > >;
                       */

// TODO: globally or hierarchical?
// pre-filter candidates by their occurences independent of their chemical suitability
template<typename TLocations>
void pre_frequency_filter(TLocations & locations, float occurrence_freq)
{

}

// post-filter candidates fulfilling chemical constraints by their relative frequency
void post_frequency_filter(/*candidates with seq,*/ float occurrence_freq)
{

}

// filter k-mers by frequency and chemical properties
template<typename io_config, typename primer_config, typename TLocations>
void filter(io_config & io_cfg, primer_config & primer_cfg, TLocations & locations)
{
    pre_frequency_filter<TLocations>(locations, primer_cfg.get_occurence_freq());

    post_frequency_filter(/*locations, */ primer_cfg.get_occurence_freq());

}

}  // namespace priset
