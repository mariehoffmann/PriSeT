// ============================================================================
//             PriSet - Primer Search Tool for metagenomic Analyses
// ============================================================================

#pragma once

#include <string>

/*!\file
 * \brief A database connector for BLAST.
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 */

namespace priset
{
/*
 * output format 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    see here: http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
*/
struct blast
{

    struct result_builder
    {

    private:
        std::string _accession;
        float _pident;
        uint _length;
        uint _mismatch;
        uint _start;
        uint _end;
        double _evalue;

    public:
        // = sseqid?
        void set_accession(std::string accession) noexcept
        {
            _accession = accession;
        }
        // set percentage of identical matches
        void set_pident(float pident) noexcept
        {
            _pident = pident;
        }

        // set number of mismatches

        // set start and end position in sequence
        void set_start(uint start) noexcept
        {
            _start = start;
        }

        // set start and end position in sequence
        void set_end(uint end) noexcept
        {
            _end = end;
        }

        // set expect value see http://www.metagenomics.wiki/tools/blast/evalue
        void set_evalue(double evalue) noexcept
        {
            _evalue = evalue;
        }

    };

    result_builder result;

    struct query_builder
    {
    private:
        std::string query = "blastn -outfmt 6";

    public:
        // set path to file with accession
        void set_accession_file()
        {
            // throw exception if file does not exist
        }

        // launch query
        void launch()
        {


        }
    };

    query_builder query;
};

} // namespace priset
