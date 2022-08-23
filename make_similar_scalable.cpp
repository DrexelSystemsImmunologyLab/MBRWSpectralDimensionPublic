//
//  main.cpp
//  graph_randomizer
//
//  Created by Adam Craig on 6/7/19.
//  Copyright © 2019 Drexel Systems Immunology Laboratory. All rights reserved.
//
//  This code is optimized for a very specialized use case:
//  We want to randomize a simple undirected graph while keeping it connected and with minimum degree 2.
//  We assume we are starting with a graph meeting these criteria and disallow any change that results in a graph violating them.

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <boost/lexical_cast.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_01.hpp>

enum RandType { original, mean_degree, all_degrees, communities, intra_edges, inter_edges };

size_t get_num_coms(std::vector<size_t> com_assignments) {
    size_t num_coms = 0;
    for (size_t i = 0; i < com_assignments.size(); ++i) {
        num_coms = std::max( num_coms, com_assignments[i] );
    }
    num_coms++;
    return num_coms;
}

void rescale_com_assignments(std::vector<size_t> old_com_assignments, size_t target_num_nodes, size_t num_coms, std::vector<size_t> & new_com_assignments) {
    // Figure out the old community sizes.
    std::vector<size_t> old_com_sizes(num_coms,0);
    size_t old_num_nodes = old_com_assignments.size();
    for (size_t i = 0; i < old_com_assignments.size(); ++i) {
        old_com_sizes[ old_com_assignments[i] ]++;
    }
    // Make the new community sizes proportional to the old,
    // but round up to the nearest integer greater or equal
    // so that each community that was originally non-empty has at least one node.
    std::vector<size_t> new_com_sizes(num_coms);
    double num_nodes_ratio = ( (double) target_num_nodes )/( (double) old_num_nodes );
    for (size_t i = 0; i < old_com_sizes.size(); ++i) {
        new_com_sizes[i] = (size_t) std::ceil( (double) old_com_sizes[i] * num_nodes_ratio );
    }
    // Recalculate the total number of nodes, as it may not be exactly the target amount.
    size_t new_num_nodes = 0;
    for (size_t i = 0; i < new_com_sizes.size(); ++i) {
        new_num_nodes += new_com_sizes[i];
    }
    // Assign consecutively numbered blocks of nodes to the same community.
    new_com_assignments.resize(new_num_nodes);
    size_t n = 0;
    for (size_t c = 0; c < new_com_sizes.size(); ++c) {
        for (size_t i = 0; i < new_com_sizes[c]; ++i) {
            new_com_assignments[n] = c;
            n++;
        }
    }
}

// Make a vector such that edge_counts[i*num_coms+j] is the number of edges running between communities i and j,
// or twice the number if i==j.
std::vector<size_t> get_edge_counts(std::vector<size_t> sources, std::vector<size_t> targets, std::vector<size_t> com_assignments, size_t num_coms) {
    std::vector<size_t> edge_counts(num_coms*num_coms,0);
    size_t num_edges = sources.size();
    for (size_t edge = 0; edge < num_edges; ++edge) {
        size_t s = sources[edge];
        size_t t = targets[edge];
        size_t s_com = com_assignments[s];
        size_t t_com = com_assignments[t];
        edge_counts[s_com*num_coms+t_com]++;
        edge_counts[t_com*num_coms+s_com]++;
    }
    return edge_counts;
}

std::vector<size_t> rescale_edge_counts(std::vector<size_t> old_edge_counts, size_t old_num_edges, size_t target_num_edges, size_t num_coms) {
    size_t num_com_pairs = old_edge_counts.size();
    std::vector<size_t> new_edge_counts(num_com_pairs);
    double num_edges_ratio = ( (double) target_num_edges )/( (double) old_num_edges );
    for (size_t i = 0; i < num_com_pairs; ++i) {
        new_edge_counts[i] = (size_t) std::ceil( (double) old_edge_counts[i] * num_edges_ratio );
    }
    // Re-do the ones on the diagonal, since each is twice the number of intra-edges in a community.
    // We divide by 2 to get the number of edges in the old community.
    // We then rescale and ceil to get the desired number of edges in the new community.
    // We then double that.
    for (size_t c = 0; c < num_coms; ++c) {
        size_t diagonal_index = c*num_coms + c;
        new_edge_counts[diagonal_index] = (size_t) 2*std::ceil( (double) old_edge_counts[diagonal_index]/2 * num_edges_ratio );
    }
    return new_edge_counts;
}

size_t get_total_num_edges(std::vector<size_t> edge_counts) {
    size_t twice_num_edges = 0;
    for (size_t i = 0; i < edge_counts.size(); ++i) {
        twice_num_edges += edge_counts[i];
    }
    return twice_num_edges/2;
}

void rescale_edges( std::vector<size_t> old_sources, std::vector<size_t> old_targets, std::vector<size_t> old_com_assignments, size_t target_num_nodes, size_t target_num_edges, boost::random::mt19937 * rng, std::vector<size_t> & new_sources, std::vector<size_t> & new_targets, std::vector<size_t> & new_com_assignments ) {
    size_t num_coms = get_num_coms(old_com_assignments);
    std::cout << "num communities: " << num_coms << "\n";
    size_t old_num_edges = old_sources.size();
    size_t old_num_nodes = old_com_assignments.size();
    // Figure out the old community sizes.
    std::cout << "finding old community sizes...\n";
    std::vector<size_t> old_com_sizes( num_coms, 0);
    for (size_t i = 0; i < old_com_assignments.size(); ++i) {
        old_com_sizes[ old_com_assignments[i] ]++;
    }
    // Make the new community sizes proportional to the old.
    // Require each community to have at least 3 nodes so that we can make a loop through them.
    // This ensures that all nodes have degree at least 2.
    std::cout << "finding new community sizes...\n";
    std::vector<size_t> new_com_sizes(num_coms);
    double num_nodes_ratio = ( (double) target_num_nodes )/( (double) old_num_nodes );
    std::cout << "num nodes ratio: " << target_num_nodes << "/" << old_num_nodes << "=" << num_nodes_ratio << "\n";
    for (size_t i = 0; i < old_com_sizes.size(); ++i) {
        size_t num_desired = (size_t) std::round( (double) old_com_sizes[i] * num_nodes_ratio );
        new_com_sizes[i] = std::max( num_desired, (size_t) 3 );
        // std::cout << "com " << i << ": " << old_com_sizes[i] << "->" << new_com_sizes[i] << "\n";
    }
    // Recalculate the total number of nodes, as it may not be exactly the target amount.
    size_t new_num_nodes = 0;
    for (size_t i = 0; i < new_com_sizes.size(); ++i) {
        new_num_nodes += new_com_sizes[i];
    }
    std::cout << "num nodes: " << new_num_nodes << "\n";
    // Assign consecutively numbered blocks of nodes to the same community.
    std::cout << "setting new community assignments...\n";
    new_com_assignments.resize(new_num_nodes);
    std::vector<size_t> com_start_node(num_coms);
    size_t n = 0;
    for (size_t c = 0; c < new_com_sizes.size(); ++c) {
        com_start_node[c] = n;
        for (size_t i = 0; i < new_com_sizes[c]; ++i) {
            new_com_assignments[n] = c;
            n++;
        }
    }
    // Make the new number of edges between two communities proportional to the old.
    // Use ceil() so that any communities that were connected in the old graph still have at least 1 edge between them.
    // If we end up with a number greater than the number of distinct edges we can make between the communities,
    // use that maximum number possible instead.
    std::cout << "finding old inter-community intra-community edge counts...\n";
    std::vector<size_t> old_edge_counts = get_edge_counts(old_sources, old_targets, old_com_assignments, num_coms);
    std::cout << "finding new inter-community edge counts...\n";
    size_t num_com_pairs = old_edge_counts.size();
    std::vector<size_t> new_edge_counts(num_com_pairs);
    double num_edges_ratio = ( (double) target_num_edges )/( (double) old_num_edges );
    std::cout << "num edges ratio: " << target_num_edges << "/" << old_num_edges << "=" << num_edges_ratio << "\n";
    for (size_t s_com = 0; s_com < num_coms; ++s_com) {
        for (size_t t_com = 0; t_com < num_coms; ++t_com) {
            size_t com_pair_index = s_com*num_coms + t_com;
            size_t max_desired = (size_t) std::ceil( (double) old_edge_counts[com_pair_index] * num_edges_ratio );
            size_t max_possible = new_com_sizes[s_com]*new_com_sizes[t_com];
            new_edge_counts[com_pair_index] = std::min(max_desired,max_possible);
            // if (s_com != t_com) {
                // std::cout << "com " << s_com << " <-> com " << t_com << ": " << old_edge_counts[com_pair_index] << "->" << new_edge_counts[com_pair_index] << ", desired: " << max_desired << ", max possible: " << max_possible << "\n";
                // std::cout << new_edge_counts[com_pair_index] << "+";
            // }
        }
    }
    // Re-do the ones on the diagonal, since each is twice the number of intra-edges in a community.
    // We divide by 2 to get the number of edges in the old community.
    // We then rescale and ceil to get the desired number of edges in the new community.
    // We then double that.
    // Also, the maximum intra-edges possible is com_size*(com_size-1)/2 instead of com_size*com_size.
    // We also want to require at least com_size edges so that we can make a loop through all of them.
    std::cout << "finding new intra-edge counts...\n";
    for (size_t c = 0; c < num_coms; ++c) {
        size_t diagonal_index = c*num_coms + c;
        size_t max_desired = (size_t) 2*std::ceil( (double) old_edge_counts[diagonal_index]/2 * num_edges_ratio );
        size_t max_possible = new_com_sizes[c]*( new_com_sizes[c] - 1 );
        size_t min_required = 2*new_com_sizes[c];
        new_edge_counts[diagonal_index] = std::max( std::min(max_desired,max_possible), min_required );
        // std::cout << "in com " << c << ": "  <<  old_edge_counts[diagonal_index] << "->" << new_edge_counts[diagonal_index] << ", desired: " << max_desired << ", possible: " << max_possible << ", required: " << min_required << "\n";
        // std::cout << new_edge_counts[diagonal_index] << "+";
    }

    // Count the total number of edges we plan to make.
    size_t twice_num_edges = 0;
    for (size_t i = 0; i < new_edge_counts.size(); ++i) {
        twice_num_edges += new_edge_counts[i];
    }
    size_t new_num_edges = twice_num_edges/2;
    std::cout << "new num edges: " << new_num_edges << "\n";

    new_sources.resize(new_num_edges);
    new_targets.resize(new_num_edges);
    size_t edge = 0;

    // First, we deterministically create a skeleton for the network:
    // Create a ring through the nodes of each community
    // and connect the first nodes of any pair of communities that should be connected.
    // This ensures that each node has degree at least 2 and that the network is connected.

    // Add a ring through the nodes of each community.
    std::cout << "Adding circuits through communities...\n";
    for (int c = 0; c < num_coms; ++c) {
        size_t first_node = com_start_node[c];
        size_t last_node = com_start_node[c] + new_com_sizes[c] - 1;
        for (int n = first_node; n < last_node; ++n) {
            new_sources[edge] = n;
            new_targets[edge] = n+1;
            edge++;
        }
        new_sources[edge] = first_node;
        new_targets[edge] = last_node;
        edge++;
    }

    // Connect the first nodes of any two communities that should be connected.
    std::cout << "Connecting connected communities...\n";
    std::vector<bool> is_connected(num_com_pairs);
    for (size_t s_com = 0; s_com < num_coms-1; ++s_com) {
        size_t first_s = com_start_node[s_com];
        for (size_t t_com = s_com+1; t_com < num_coms; ++t_com) {
            size_t first_t = com_start_node[t_com];
            size_t com_pair_index = s_com*num_coms + t_com;
            bool ic = new_edge_counts[s_com*num_coms+t_com] > 0;
            is_connected[com_pair_index] = ic;
            if (ic) {
                new_sources[edge] = first_s;
                new_targets[edge] = first_t;
                edge++;
            }
        }
    }

    // For each pair of communities, find the probability that two nodes will have an edge between them,
    // given their community assignments.
    // Choose this value so that the expected number of edges is the desired number,
    // taking into account the edges we already added.
    std::cout << "Computing probabilities of additional inter-edges...\n";
    std::vector<double> prob_connect(num_com_pairs);
    for (int s_com = 0; s_com < num_coms; ++s_com) {
        double s_num_nodes = (double) new_com_sizes[s_com];
        for (int t_com = 0; t_com < num_coms; ++t_com) {
            size_t com_pair_index = s_com*num_coms+t_com;
            double t_num_nodes = (double) new_com_sizes[t_com];
            double ic_double = (double) is_connected[com_pair_index];
            // Note that we set the community sizes to be at least 3, so we do not worry about singleton communities.
            prob_connect[s_com*num_coms+t_com] = ( (double) new_edge_counts[com_pair_index] - ic_double )/( s_num_nodes * t_num_nodes - ic_double );
        }
    }
    // Handle the diagonals of the matrix separately.
    std::cout << "Computing probabilities of additional intra-edges...\n";
    for (int c = 0; c < num_coms; ++c) {
        size_t com_pair_index = c*num_coms + c;
        double c_num_nodes = (double) new_com_sizes[c];
        prob_connect[com_pair_index] = ( (double) new_edge_counts[com_pair_index]/2 - c_num_nodes )/( c_num_nodes*(c_num_nodes-1)/2.0 - c_num_nodes );
    }

    boost::random::uniform_01<double> yes_no_dist;
    // If we end up with more edges than expected, increase the vector sizes to accommodate the additional edges.
    // It is generally more efficient to increase the sizes by a few more at once instead of 1 each time.
    const size_t size_increment = 100;

    // For each community, randomly decide whether to add each of the intra-edges we have not already added.
    std::cout << "Randomly adding intra-edges...\n";
    for (size_t c = 0; c < num_coms; ++c) {
        size_t first_node = com_start_node[c];
        size_t last_node = com_start_node[c] + new_com_sizes[c] - 1;
        double pc_c = prob_connect[c*num_coms+c];
        // Handle edges connecting to the first node separately.
        // Instead of connecting to the one before it and the one after it,
        // it connects to the second and last nodes.
        for (size_t t = first_node+2; t <= last_node-1; ++t) {
            if ( yes_no_dist(*rng) <= pc_c ) {
                if ( edge == new_sources.size() ) {
                    new_sources.resize( new_sources.size() + size_increment );
                    new_targets.resize( new_targets.size() + size_increment );
                }
                new_sources[edge] = first_node;
                new_targets[edge] = t;
                // std::cout << edge << ": " << first_node << "<->" << t << "\n";
                edge++;
            }
        }
        for (size_t s = first_node+1; s <= last_node-2; ++s) {
            for (size_t t = s+2; t <= last_node; ++t) {
                if ( yes_no_dist(*rng) <= pc_c ) {
                    if ( edge == new_sources.size() ) {
                        new_sources.resize( new_sources.size() + size_increment );
                        new_targets.resize( new_targets.size() + size_increment );
                    }
                    new_sources[edge] = s;
                    new_targets[edge] = t;
                    // std::cout << edge << ": " << s << "<->" << t << "\n";
                    edge++;
                }
            }
        }
    }
    // Randomly decide whether to include each of the inter-edges we have not already added.
    std::cout << "Randomly adding inter-edges...\n";
    for (size_t s_com = 0; s_com < num_coms-1; ++s_com) {
        size_t first_s = com_start_node[s_com];
        size_t last_s = first_s + new_com_sizes[s_com] - 1;
        for (size_t t_com = s_com+1; t_com < num_coms; ++t_com) {
            size_t first_t = com_start_node[t_com];
            size_t last_t = first_t + new_com_sizes[t_com] - 1;
            double pc_c = prob_connect[s_com*num_coms+t_com];
            // std::cout << "com " << s_com << " <-> " << " com " << t_com << ", edge indices: " << edge << " through " << stop_at_edge-1 << ", " << new_edge_counts[s_com*num_coms+t_com] << " edges" << "\n";
            // If the communities are connected, we have already connected the first nodes.
            for (size_t s = first_s+1; s <= last_s; ++s) {
                for (size_t t = first_t+1; t <= last_t; ++t) {
                    if ( yes_no_dist(*rng) <= pc_c ) {
                        if ( edge == new_sources.size() ) {
                            new_sources.resize( new_sources.size() + size_increment );
                            new_targets.resize( new_targets.size() + size_increment );
                        }
                        new_sources[edge] = s;
                        new_targets[edge] = t;
                        // std::cout << edge << ": " << s << "<->" << t << "\n";
                        edge++;
                    }
                }
            }
        }
    }

    size_t final_num_edges = edge;
    std::cout << "Resizing to final edge count: " << final_num_edges << "\n";
    new_sources.resize(final_num_edges);
    new_targets.resize(final_num_edges);
    std::cout << "done creating rescaled network\n";
    // If the original graph was connected, we should now have a new connected graph.
    // Each community is a connected subgraph,
    // and any two communities that were connected in the original are still connected.
}

bool read_string_param(int argc, const char * argv[], const char * name, std::string & value) {
    bool found_value = false;
    for (size_t arg = 0; (arg < argc-1) && !found_value; arg++) {
        if ( strcmp(argv[arg], name) == 0 ) {
            value = std::string( argv[arg+1] );
            found_value = true;
        }
    }
    return found_value;
}

bool read_int_param(int argc, const char * argv[], const char * name, int & value) {
    std::string string_value;
    bool found_value = read_string_param(argc, argv, name, string_value);
    value = boost::lexical_cast<int>(string_value);
    return found_value;
}

bool read_size_param(int argc, const char * argv[], const char * name, size_t & value) {
    std::string string_value;
    bool found_value = read_string_param(argc, argv, name, string_value);
    value = boost::lexical_cast<size_t>(string_value);
    return found_value;
}

bool read_rand_type_param(int argc, const char * argv[], const char * name, RandType & value) {
    std::string string_value;
    bool found_value = read_string_param(argc, argv, name, string_value);
    if (string_value == "mean_degree") {
        value = mean_degree;
    }
    else if (string_value == "all_degrees") {
        value = all_degrees;
    }
    else if (string_value == "communities") {
        value = communities;
    }
    else if (string_value == "intra_edges") {
        value = intra_edges;
    }
    else if (string_value == "inter_edges") {
        value = inter_edges;
    }
    else {
        std::cerr << "Unrecognized randomization type " << string_value << ".\n";
    }
    return found_value;
}

// We assume that the edges are undirected so each connected pair is only listed once.
bool read_edges(std::string file_name, std::vector<size_t> & sources, std::vector<size_t> & targets) {
    bool is_valid;
    const std::string delimiter = "\t";
    size_t num_edges = 0;
    std::ifstream edge_file (file_name);
    std::string line;
    // std::cout << "reading edge list from " << file_name << "\n";
    if ( edge_file.is_open() ) {
        // First, count how many lines/edges we have.
        // Suppose each line is valid until
        // we find one without the delimiter.
        // std::cout << "counting edges\n";
        is_valid = true;
        while ( std::getline(edge_file, line) ) {
            // std::cout << line << "\n";
            // Check that each contains the delimiter.
            if ( line.find(delimiter) == std::string::npos ) {
                std::cerr << "Edge file " << file_name << " contained misformatted line " << line << ".\n";
                is_valid = false;
                break;
            }
            num_edges += 1;
        }
        // std::cout << "found " << num_edges << " edges\n";
        if (is_valid) {
            // std::cout << "reading to edge list\n";
            // Pre-allocate a vector of appropriate size.
            sources.resize(num_edges);
            targets.resize(num_edges);
            // Reset to the beginning of the file.
            edge_file.clear();
            edge_file.seekg(0);
            // Read in the results,
            // this time saving the results to the edge list.
            size_t edge = 0;
            while ( std::getline(edge_file, line) ) {
                // std::cout << line << "\n";
                // Find the first occurrence of the tab delimiter.
                size_t tab_pos = line.find(delimiter);
                // std::cout << line.substr(0, tab_pos) << "\t" << line.substr(tab_pos+1) << "\n";
                size_t node1 = boost::lexical_cast<size_t>( line.substr(0, tab_pos) ) - 1;
                size_t node2 = boost::lexical_cast<size_t>( line.substr(tab_pos+1) ) - 1;
                if (node1 < node2) {
                    sources[edge] = node1;
                    targets[edge] = node2;
                }
                else {
                    sources[edge] = node2;
                    targets[edge] = node1;
                }
                edge++;
            }
        }
        edge_file.close();
    }
    else {
        std::cerr << "Failed to open file " << file_name << ".\n";
        is_valid = false;
    }
    return is_valid;
}

// Assume one assignment per node.
bool read_communities(std::string file_name, std::vector<size_t> & com_assignments) {
    bool is_valid;
    const std::string delimiter = " ";
    const std::string commenter = "#";
    size_t num_nodes = 0;
    std::ifstream com_file (file_name);
    std::string line;
    // std::cout << "reading community assignments from " << file_name << "\n";
    if ( com_file.is_open() ) {
        // First, count how many lines/assignments we have.
        // Suppose each line is valid until
        // we find one without the delimiter.
        // std::cout << "counting assignments\n";
        is_valid = true;
        while ( std::getline(com_file, line) ) {
            // std::cout << line << "\n";
            // Check that each contains the delimiter.
            if (  ( line.find(commenter) != 0 ) && ( line.find(delimiter) == std::string::npos )  ) {
                std::cerr << "Community file " << file_name << " contained misformatted line " << line << ".\n";
                is_valid = false;
                break;
            }
            // Only count non-comment lines.
            num_nodes += ( line.find(commenter) != 0 );
        }
        // std::cout << "found " << num_nodes << " nodes\n";
        if (is_valid) {
            // std::cout << "reading to edge list\n";
            // Pre-allocate a vector of appropriate size.
            com_assignments.resize(num_nodes);
            // Reset to the beginning of the file.
            com_file.clear();
            com_file.seekg(0);
            // Read in the results,
            // this time saving the results to the edge list.
            while ( std::getline(com_file, line) ) {
                // Skip comments, which start with #.
                if ( line.find(commenter) != 0 ) {
                    // std::cout << line << "\n";
                    // Find the first occurrence of the tab delimiter.
                    size_t tab_pos = line.find(delimiter);
                    // std::cout << line.substr(0, tab_pos) << "\t" << line.substr(tab_pos+1) << "\n";
                    size_t node = boost::lexical_cast<size_t>( line.substr(0, tab_pos) ) - 1;
                    size_t next_tab_pos = line.find(delimiter, tab_pos+1);
                    size_t community = boost::lexical_cast<size_t>( line.substr(tab_pos+1, next_tab_pos - tab_pos - 1) ) - 1;
                    com_assignments[node] = community;
                }
            }
        }
        com_file.close();
    }
    else {
        std::cerr << "Failed to open file " << file_name << ".\n";
        is_valid = false;
    }
    return is_valid;
}

bool write_edges(std::string file_name, std::vector<size_t> sources, std::vector<size_t> targets) {
    bool is_valid;
    if ( sources.size() != targets.size() ) {
        std::cerr << "write_edges received list of sources of length " << sources.size() << " and list of targets of length " << targets.size() << "\n";
        return false;
    }
    std::ofstream edge_file(file_name);
    if ( edge_file.is_open() ) {
        for (size_t edge = 0; edge < sources.size(); edge++) {
            edge_file << sources[edge] + 1 << "\t" << targets[edge] + 1 << "\n";
        }
        if ( edge_file.fail() ) {
            std::cerr << strerror(errno);
            is_valid = false;
        }
        else {
            is_valid = true;
        }
        edge_file.close();
    }
    else {
        std::cerr << "Failed to open file " << file_name << "\n";
        is_valid = false;
    }
    return is_valid;
}

bool write_communities(std::string file_name, std::vector<size_t> com_assignments) {
    bool is_valid;
    std::ofstream edge_file(file_name);
    if ( edge_file.is_open() ) {
        for (size_t node = 0; node < com_assignments.size(); node++) {
            edge_file << node + 1 << "\t" << com_assignments[node] + 1 << "\n";
        }
        if ( edge_file.fail() ) {
            std::cerr << strerror(errno);
            is_valid = false;
        }
        else {
            is_valid = true;
        }
        edge_file.close();
    }
    else {
        std::cerr << "Failed to open file " << file_name << "\n";
        is_valid = false;
    }
    return is_valid;
}

int main(int argc, const char * argv[]) {
    std::cout << "randomize graph:\n";
    int rand_seed;
    bool found_rand_seed = read_int_param(argc, argv, "-r", rand_seed);
    if (!found_rand_seed) {
        std::cerr << "Argument -r (random number generator seed) is required.\n";
        return 1;
    }
    std::string original_edge_file_name;
    bool found_original_file = read_string_param(argc, argv, "-g", original_edge_file_name);
    if (!found_original_file) {
        std::cerr << "Argument -g (original graph file) is required.\n";
        return 1;
    }
    std::string community_file_name;
    bool found_community_file = read_string_param(argc, argv, "-c", community_file_name);
    if (!found_community_file) {
        std::cerr << "Argument -c (community assignment file) is required.\n";
        return 1;
    }
    size_t target_num_nodes;
    bool found_target_num_nodes = read_size_param(argc, argv, "-n", target_num_nodes);
    if (!found_target_num_nodes) {
        std::cerr << "Argument -n (target number of nodes in the new graph) is required.\n";
        return 1;
    }
    std::string rand_edge_file_name;
    bool found_rand_file = read_string_param(argc, argv, "-o", rand_edge_file_name);
    if (!found_rand_file) {
        std::cerr << "Argument -o (rescaled graph file) is required.\n";
        return 1;
    }
    std::string rand_com_file_name;
    bool found_rand_com_file = read_string_param(argc, argv, "-k", rand_com_file_name);
    if (!found_rand_com_file) {
        std::cerr << "Argument -k (rescaled community assignments file) is required.\n";
        return 1;
    }
    std::cout << "found params\n";
    std::cout << "original graph file: " << original_edge_file_name << "\n";
    std::cout << "original community assignment file: " << community_file_name << "\n";
    std::cout << "rescaled graph file: " << rand_edge_file_name << "\n";
    std::cout << "rescaled community assignment file: " << rand_com_file_name << "\n";
    std::cout << "target num nodes: " << target_num_nodes << "\n";
    std::cout << "rand seed: " << rand_seed << "\n";
    boost::random::mt19937 rng (rand_seed);
    std::vector<size_t> old_sources;
    std::vector<size_t> old_targets;
    bool edge_read_succeeded = read_edges(original_edge_file_name, old_sources, old_targets);
    if (!edge_read_succeeded) {
        return 1;
    }
    std::cout << "read in edges\n";
    std::vector<size_t> old_com_assignments;
    bool com_read_succeeded = read_communities(community_file_name, old_com_assignments);
    if (!com_read_succeeded) {
        return 1;
    }
    std::cout << "read in community assignments\n";
    size_t old_num_edges = old_sources.size();
    std::cout << "old num edges: " << old_num_edges << "\n";
    // Find the number of nodes.
    // This assumes the nodes are consecutively numbered and connected, so we only need to find the maximum node index.
    // This assumes that the graph is never empty and always connected, so an empty edge list maps to a singleton node.
    size_t old_num_nodes = 0;
    for (size_t edge = 0; edge < old_num_edges; ++edge) {
        if (old_sources[edge] > old_num_nodes) {
            old_num_nodes = old_sources[edge];
        }
        if (old_targets[edge] > old_num_nodes) {
            old_num_nodes = old_targets[edge];
        }
    }
    old_num_nodes++;
    std::cout << "old num nodes: " << old_num_nodes << "\n";
    size_t target_num_edges = (size_t) std::round(  ( (double) old_num_edges )*( (double) target_num_nodes )/( (double) old_num_nodes )  );
    std::cout << "target num edges: " << target_num_edges << "\n";
    std::vector<size_t> new_sources;
    std::vector<size_t> new_targets;
    std::vector<size_t> new_com_assignments;
    rescale_edges( old_sources, old_targets, old_com_assignments, target_num_nodes, target_num_edges, &rng, new_sources, new_targets, new_com_assignments );
    bool write_succeeded;
    write_succeeded = write_edges(rand_edge_file_name, new_sources, new_targets);
    if (!write_succeeded) {
        return 1;
    }
    std::cout << "wrote rescaled graph\n";
    write_succeeded = write_communities(rand_com_file_name, new_com_assignments);
    if (!write_succeeded) {
        return 1;
    }
    std::cout << "wrote rescaled community assignments\n";
    return 0;
}
