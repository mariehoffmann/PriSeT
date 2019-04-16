#include <algorithm>
#include <unordered_map>
#include <vector>

/*
 * A taxonomy tree stores the hierarchy of taxonomic nodes given their integer
 * identifier.
 * Implementation details: It is stored as an unordered_map with the taxonomic
 * identifier of the parent as key and the list of child IDs as a sorted vector
 * to support binary search.
 */
struct taxonomy
{
    using key_type = unsigned int;
    using node_type = unordered_map<key_type, std::vector<key_type>>::value_type;
    using node_map_type = typename unordered_map<key_type, std::vector<key_type>>;
    using iterator_type = node_map_type::iterator;

    // tax file is comma separated tuple list [parent_id, child_id]
    taxonomy(std::string & tax_file)
    {
        std::ifstream ifs(tax_file, std::ifstream::in);
        if (!ifs.is_open())
            std::cout << "Error: could not open taxonomy file, namely '" << tax_file << "'\n", exit(0);
        std::cout << "ifs is open: " << ifs.is_open() << std::endl;

        // store flat, i.e., [pid1, cid1, pid2, cid2, ...]
        std::vector<key_type> raw_taxa;
        std::string cell = "";
        // Iterate through each line and split the content using delimeter
        while (getline(istr, cell, ','))
        {
            std::cout << cell << ", ";
            raw_taxa.push_back(cell);
        }
        std::cout << std::endl;
        std::cout << "raw_taxa.size() = " << raw_taxa.size() << std::endl;
        ifs.close();
        build_tree(raw_taxa);
        // destroy raw_taxa
        ~raw_taxa;
    };

private:

    void insert_node(key_type const parent_id, key_type const node_id)
    {
        iterator_type it = node_map.find(parent_id);
        if (it == node_map.end())
            node_map.insert(node_type{parent_id, std::vector<key_type>(node_id);
        else
            it->second.push_back(node_id);
    }

    void build_tree(std::vector<key_type> const & raw_taxa)
    {
        // for node pair in file apply insert_node

        // sort node lists
        for (iterator_type it = node_map.begin(); it != node_map.end(); ++it)
        {
            std::sort(it->second.begin(), it->second.end());
        }
    }

    key_type root_node;
    node_map_type node_map;

};
