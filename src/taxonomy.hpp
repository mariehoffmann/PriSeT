// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
//          Manual: https://github.com/mariehoffmann/PriSeT

#pragma once

#include <algorithm>
#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stack>
#include <unordered_map>
#include <vector>

namespace fs = std::experimental::filesystem;

namespace priset
{
/*
 * A taxonomy tree stores the hierarchy of taxonomic nodes given their integer
 * identifier. The taxonomy might be incomplete in terms of that there are multiple
 * root nodes. The taxonomy is retrieved from a text file where the first column
 * corresponds to the taxonomy identifier and the second column to its parent.
 * Implementation details: It is stored as an unordered_map with the taxonomic
 * identifier of the parent as key and an std::set of child IDs.
 */
struct taxonomy
{
    using key_type = unsigned int;
    using value_type = std::set<key_type>;
    using node_map_type = typename std::unordered_map<key_type, value_type>;
    using iterator_type = node_map_type::iterator;

    taxonomy(fs::path const & tax_file)
    {
        std::ifstream ifs{tax_file, std::ifstream::in};
        if (!ifs.is_open())
            std::cout << "Error: could not open taxonomy file, namely '" << tax_file << "'\n", exit(0);

        // Temporary set of all taxids assigned as child (for root list determination)
        std::set<key_type> someones_child{};
        // Iterate through tax_file and split content using delimeter
        std::string line, taxid_str;
        key_type taxid, p_taxid;
        while (std::getline(ifs, line))
        {
            std::istringstream iline(line);
            std::getline(iline, taxid_str, ',');
            taxid = char2key(taxid_str);
            std::getline(iline, taxid_str, ',');
            p_taxid = char2key(taxid_str);
            someones_child.insert(taxid);
            insert_node(taxid, p_taxid);
        }
        ifs.close();
        // determine root nodes, i.e. keys that are not in the children set
        std::vector<int> keys;
        std::transform(node_map.begin(), node_map.end(), std::back_inserter(keys),
            [](const node_map_type::value_type &keyval){return keyval.first;});
        for (key_type key : keys)
            if (someones_child.find(key) == someones_child.end())
                root_nodes.push_back(key);
    };

    // Display taxonomy by pre-order traversing tree/forest.
    void print_taxonomy()
    {
        std::stack<std::pair<key_type, key_type>> lifo;
        using tree_node_type = typename std::pair<key_type, key_type>;
        unsigned short tree_id = 0;
        unsigned short level = 0;
        for (key_type root : root_nodes)
        {
            std::cout << "################ Tree " << tree_id++ << " ################\n";
            level = 0;
            lifo.push(std::make_pair(root, static_cast<key_type>(0)));
            while (!lifo.empty())
            {
                    tree_node_type node = lifo.top();
                    lifo.pop();
                    level = node.second;
                    std::cout << "\n" << std::string(level, '\t') << "|__\t";
                    std::cout << node.first;
                    // push back children
                    for (key_type c : node_map[node.first])
                        lifo.push(std::make_pair(c, level + 1));
            }
            std::cout << std::endl;
        }
    }

private:

    key_type char2key(std::string taxid_str)
    {
        std::stringstream taxid_ss;
        taxid_ss << taxid_str;
        key_type taxid;
        taxid_ss >> taxid;
        return taxid;
    }

    void insert_node(key_type const node_id, key_type const parent_id)
    {
        iterator_type it = node_map.find(parent_id);
        if (it == node_map.end())
            node_map.insert({{parent_id, value_type{node_id}}});
        else
            it->second.insert(node_id);
    }

    // taxid of root node(s)
    std::vector<key_type> root_nodes{};
    // Store hierarchy as map over set {p_taxid: (c_taxid1, c_taxid2, ...)}
    node_map_type node_map;

};

}  // namespace priset
