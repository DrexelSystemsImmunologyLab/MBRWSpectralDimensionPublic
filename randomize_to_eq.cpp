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

enum RandType { original, mean_degree, all_degrees, communities, intra_edges, inter_edges };

class Graph {
    
    const size_t ABSENT = SIZE_MAX;
    size_t num_nodes;
    size_t num_edges;
    std::vector<size_t> sources;
    std::vector<size_t> targets;
    std::vector<size_t> degrees;
    std::vector<size_t> adjacency_matrix;
    
    size_t get_loc(size_t source, size_t target);
    // These do not perform any checks.
    void remove_edge_at(size_t edge_index);
    void add_edge_at(size_t edge_index, size_t source, size_t target);
    void replace_edge_at(size_t edge_index, size_t source, size_t target);
    
public:
    Graph();
    Graph(std::vector<size_t> new_sources, std::vector<size_t> new_targets);
    bool is_present(size_t source, size_t target);
    bool is_absent(size_t source, size_t target);
    size_t get_source(size_t edge_index);
    size_t get_target(size_t edge_index);
    size_t get_num_nodes();
    size_t get_num_edges();
    bool is_connected();
    // Returns false if the replacement could not happen for one of the following reasons:
    // The edge to be added is already in the graph.
    // The edge to be added is a self-loop (new_source == new_target).
    // The replacement would drop the degree of one of the nodes below 2.
    // The replacement would make the graph not connected.
    bool replace_edge(size_t edge_index, size_t new_source, size_t new_target);
    // Returns false if the swap could not happen for one of the following reasons:
    // Edge source_a<->target_a or edge source_b<->target_b is not in the graph.
    // Edge source_a<->target_b or edge source_b<->target_a is already in the graph.
    // source_a<->target_b is a self-loop (source_a==target_b).
    // source_b<->target_a is a self-loop (source_b==target_a).
    // The replacement would make the graph unconnected.
    bool swap_endpoints(size_t edge_index_a, size_t edge_index_b);
    // Make a deep copy of the Graph.
    Graph get_copy();
};

size_t Graph::get_loc(size_t source, size_t target) {
    return num_nodes*source + target;
}

void Graph::remove_edge_at(size_t edge_index) {
    size_t source = sources[edge_index];
    size_t target = targets[edge_index];
    adjacency_matrix[ this->get_loc(source, target) ] = ABSENT;
    degrees[source]--;
    degrees[target]--;
}

void Graph::add_edge_at(size_t edge_index, size_t source, size_t target) {
    sources[edge_index] = source;
    targets[edge_index] = target;
    adjacency_matrix[ this->get_loc(source, target) ] = edge_index;
    degrees[source]++;
    degrees[target]++;
}

void Graph::replace_edge_at(size_t edge_index, size_t new_source, size_t new_target) {
    this->remove_edge_at(edge_index);
    this->add_edge_at(edge_index, new_source, new_target);
}

Graph::Graph() {}

Graph::Graph(std::vector<size_t> new_sources, std::vector<size_t> new_targets) {
    num_edges = new_sources.size();
    sources.resize(num_edges);
    targets.resize(num_edges);
    // std::cout << "allocated space for " << num_edges << " edges:\n";
    size_t max_node = 0;
    for (size_t edge = 0; edge < num_edges; edge++) {
        size_t source = new_sources[edge];
        size_t target = new_targets[edge];
        // std::cout << "edge[" << edge << "] = (" << source << "," << target << ")\n";
        if (source > max_node) {
            max_node = source;
        }
        if (target > max_node) {
            max_node = target;
        }
    }
    num_nodes = max_node+1;
    // std::cout << "found max node index " << max_node << "\n";
    degrees.resize(num_nodes);
    adjacency_matrix.resize(num_nodes*num_nodes);
    // std::cout << "allocated space for " << num_nodes << " nodes\n";
    for (auto it = adjacency_matrix.begin(); it != adjacency_matrix.end(); ++it) {
        *it = ABSENT;
    }
    // std::cout << "initialized adjacency matrix as empty\n";
    // Filter out any self-loops or multi-edges.
    size_t internal_edge = 0;
    for (size_t edge= 0; edge < num_edges; edge++) {
        size_t source = new_sources[edge];
        size_t target = new_targets[edge];
        if (source > target) {
            size_t temp = target;
            target = source;
            source = temp;
        }
        bool is_not_self_loop = source != target;
        bool is_not_multi_edge = !( this->is_present(source, target) );
        if (is_not_self_loop && is_not_multi_edge) {
            this->add_edge_at(internal_edge, source, target);
            internal_edge++;
        }
    }
    num_edges = internal_edge;
    // std::cout << "found " << num_edges << " valid edges\n";
}

bool Graph::is_present(size_t source, size_t target) {
    return adjacency_matrix[ this->get_loc(source, target) ] != ABSENT;
}

bool Graph::is_absent(size_t source, size_t target) {
    return adjacency_matrix[ this->get_loc(source, target) ] == ABSENT;
}

size_t Graph::get_source(size_t edge) {
    return sources[edge];
}

size_t Graph::get_target(size_t edge) {
    return targets[edge];
}

size_t Graph::get_num_nodes() {
    return num_nodes;
}

size_t Graph::get_num_edges() {
    return num_edges;
}

bool Graph::is_connected() {
    // Create a representation of the graph as a set of neighbor lists.
    // This is easier to use for the breadth-first-search process we use to check whether the graph is connected.
    std::vector<size_t> neighbor_start_indices( num_nodes );
    neighbor_start_indices[0] = 0;
    for (size_t node = 0; node < num_nodes-1; node++) {
        neighbor_start_indices[node+1] = neighbor_start_indices[node] + degrees[node];
    }
    std::vector<size_t> neighbor_end_indices( neighbor_start_indices );
    std::vector<size_t> neighbors( 2*num_edges );
    for (size_t edge = 0; edge < num_edges; edge++) {
        size_t source = this->get_source(edge);
        size_t target = this->get_target(edge);
        neighbors[ neighbor_end_indices[source] ] = target;
        neighbor_end_indices[source]++;
        neighbors[ neighbor_end_indices[target] ] = source;
        neighbor_end_indices[target]++;
    }
    // Keep track of which nodes we have reached and how many.
    std::vector<bool> is_unreached (num_nodes, true);
    size_t num_reached = 0;
    // In an undirected graph, if we can reach all points from one node, we can reach all points from any node, so the starting point does not matter.
    // For simplicity, start at node 0.
    const size_t start_node = 0;
    is_unreached[start_node] = false;
    num_reached++;
    // Keep a queue of ones we have reached but the neighbors of which we have not checked.
    std::queue<size_t> nodes_to_explore;
    nodes_to_explore.push(start_node);
    while ( !nodes_to_explore.empty() ) {
        size_t current_node = nodes_to_explore.front();
        nodes_to_explore.pop();
        // If we have reached the current node, then we can reach all of its neighbors.
        for (size_t neighbor_index = neighbor_start_indices[current_node]; neighbor_index < neighbor_end_indices[current_node]; neighbor_index++) {
            size_t neighbor = neighbors[neighbor_index];
            if (is_unreached[neighbor]) {
                is_unreached[neighbor] = false;
                nodes_to_explore.push(neighbor);
                num_reached++;
            }
        }
    }
    // std::cout << "BFS complete: reached " << num_reached << " of " << num_nodes << " nodes\n";
    return num_reached == num_nodes;
}

bool Graph::replace_edge(size_t edge, size_t new_source, size_t new_target) {
    size_t old_source = this->get_source(edge);
    size_t old_target = this->get_target(edge);
    if (new_source > new_target) {
        size_t temp = new_target;
        new_target = new_source;
        new_source = temp;
    }
    bool new_is_absent = this->is_absent(new_source, new_target);
    // if (!new_is_absent) {
    //     std::cout << "replacement of " << old_source << " <-> " << old_target << " with " << new_source << "<->" << new_target << " creates a multi-edge\n";
    // }
    bool is_not_self_loop = new_source != new_target;
    // if (!is_not_self_loop) {
    //     std::cout << "replacement of " << old_source << " <-> " << old_target << " with " << new_source << "<->" << new_target << " creates a self-loop\n";
    // }
    bool source_deg_above_2 = (degrees[old_source] > 2) || (old_source == new_source);
    // if (!source_deg_above_2) {
    //    std::cout << "replacement of " << old_source << " <-> " << old_target << " with " << new_source << "<->" << new_target << " puts degree of " << old_source << " below 2\n";
    // }
    bool target_deg_above_2 = (degrees[old_target] > 2) || (old_target == new_target);
    // if (!target_deg_above_2) {
    //     std::cout << "replacement of " << old_source << " <-> " << old_target << " with " << new_source << "<->" << new_target << " puts degree of " << old_target << " below 2\n";
    // }
    bool result = new_is_absent && is_not_self_loop && source_deg_above_2 && target_deg_above_2;
    if (result) {
        this->replace_edge_at(edge, new_source, new_target);
        // If this replacement makes the graph unconnected, reverse it.
        if (  !( this->is_connected() )  ) {
            result = false;
            this->replace_edge_at(edge, old_source, old_target);
        }
        // else {
        //     std::cout << "replacement of " << old_source << " <-> " << old_target << " with " << new_source << "<->" << new_target << " is good\n";
        // }
    }
    return result;
}

bool Graph::swap_endpoints(size_t edge_a, size_t edge_b) {
    size_t source_a = this->get_source(edge_a);
    size_t target_a = this->get_target(edge_a);
    size_t source_b = this->get_source(edge_b);
    size_t target_b = this->get_target(edge_b);
    size_t source_c = source_a;
    size_t target_c = target_b;
    if (source_c > target_c) {
        size_t temp = target_c;
        target_c = source_c;
        source_c = temp;
    }
    size_t source_d = source_b;
    size_t target_d = target_a;
    if (source_d > target_d) {
        size_t temp = target_d;
        target_d = source_d;
        source_d = temp;
    }
    bool c_is_absent = this->is_absent(source_c, target_c);
    bool d_is_absent = this->is_absent(source_d, target_d);
    // if ( !(c_is_absent && d_is_absent) ) {
    //     std::cout << "swap of " << source_a << " <-> " << target_a << " with " << source_b << "<->" << target_b << " is creates a multi-edge\n";
    // }
    bool c_is_not_self_loop = source_c != target_c;
    bool d_is_not_self_loop = source_d != target_d;
    // if ( !(c_is_not_self_loop && d_is_not_self_loop) ) {
    //     std::cout << "swap of " << source_a << " <-> " << target_a << " with " << source_b << "<->" << target_b << " is creates a self-loop\n";
    // }
    bool result = c_is_absent && d_is_absent && c_is_not_self_loop && d_is_not_self_loop;
    if (result) {
        size_t edge_c = edge_a;
        size_t edge_d = edge_b;
        this->replace_edge_at(edge_c, source_c, target_c);
        this->replace_edge_at(edge_d, source_d, target_d);
        if (  !( this->is_connected() )  ) {
            // If this swap makes the graph unconnected, reverse it.
            result = false;
            this->replace_edge_at(edge_a, source_a, target_a);
            this->replace_edge_at(edge_b, source_b, target_b);
            // std::cout << "swap of " << source_a << " <-> " << target_a << " with " << source_b << "<->" << target_b << " makes graph unconnected\n";
        }
        // else {
        //     std::cout << "swap of " << source_a << " <-> " << target_a << " with " << source_b << "<->" << target_b << " is good\n";
        // }
    }
    return result;
}

Graph Graph::get_copy() {
    return Graph( sources, targets );
}

// Count the number of edges present in both graphs.
// Do this by iterating over the edges in a and checking which ones are also in b.
size_t get_num_shared(Graph *a, Graph *b) {
    size_t num_shared = 0;
    for (size_t edge = 0; edge < a->get_num_edges(); edge++) {
        size_t source = a->get_source(edge);
        size_t target = a->get_target(edge);
        num_shared += b->is_present(source, target);
    }
    return num_shared;
}

// For each edge, attempt num_passes swaps and num_passes replacements.
// Nodes with degree exactly 2 can only participate in swaps,
// so this approach ensures that this randomization can randomize their edges, too.
// Going through the edges in order multiple times lets us come back to ones we failed to randomize before.
// Other randomizations may have made it easier to replace an edge that was previously critical to keeping min degree >= 2 or the graph connected.
void replace_all_with_any(Graph * graph, boost::random::mt19937 * rng, size_t num_passes) {
    std::cout << "replacing edges while preserving total number of nodes and edges...\n";
    size_t num_edges = graph->get_num_edges();
    size_t num_nodes = graph->get_num_nodes();
    boost::random::uniform_int_distribution<size_t> all_edges_dist(0, num_edges-1);
    boost::random::uniform_int_distribution<size_t> all_nodes_dist(0,num_nodes-1);
    Graph copy = graph->get_copy();
    size_t min_num_shared = num_edges;
    size_t passes_since_min = 0;
    size_t pass_num = 0;
    do {
        size_t num_swaps = 0;
        size_t num_replaces = 0;
        for (size_t edge = 0; edge < num_edges; edge++) {
            size_t swap_partner = all_edges_dist(*rng);
            num_swaps += graph->swap_endpoints(edge, swap_partner);
            size_t new_source = all_nodes_dist(*rng);
            size_t new_target = all_nodes_dist(*rng);
            num_replaces += graph->replace_edge(edge, new_source, new_target);
        }
        pass_num++;
        // Stop when we have gone num_passes passes without finding a randomization further away from the original.
        // The longer we have gone since seeing a new minimum, the more likely it is that we have reached a point of diminishing returns.
        size_t num_shared = get_num_shared( graph, &copy );
        if (num_shared < min_num_shared) {
            min_num_shared = num_shared;
            passes_since_min = 0;
        } else {
            passes_since_min++;
        }
        std::cout << "pass " << pass_num << " complete: " << num_swaps << " swaps, " << num_replaces << " replacements, " << num_shared << " edges also in original\n";
    } while (passes_since_min < num_passes);
    std::cout << "done\n";
}

// Attempt num_passes swaps for each edge.
void swap_all_with_all(Graph * graph, boost::random::mt19937 * rng, size_t num_passes) {
    std::cout << "swapping endpoints while preserving degrees...\n";
    size_t num_edges = graph->get_num_edges();
    boost::random::uniform_int_distribution<size_t> all_edges_dist(0, num_edges-1);
    Graph copy = graph->get_copy();
    size_t min_num_shared = num_edges;
    size_t passes_since_min = 0;
    size_t pass_num = 0;
    do {
        size_t num_swaps = 0;
        for (size_t edge = 0; edge < num_edges; edge++) {
            size_t swap_partner = all_edges_dist(*rng);
            num_swaps += graph->swap_endpoints(edge, swap_partner);
        }
        pass_num++;
        // Stop when we have gone num_passes passes without finding a randomization further away from the original.
        // The longer we have gone since seeing a new minimum, the more likely it is that we have reached a point of diminishing returns.
        size_t num_shared = get_num_shared( graph, &copy );
        if (num_shared < min_num_shared) {
            min_num_shared = num_shared;
            passes_since_min = 0;
        } else {
            passes_since_min++;
        }
        std::cout << "pass " << pass_num << " complete: " << num_swaps << " swaps, " << num_shared << " edges also in original\n";
    } while (passes_since_min < num_passes);
    std::cout << "done\n";
}

// Attempt num_passes swaps and num_passes replacements for each edge with endpoints in different communities (inter-edge).
// Only swap it if the swap creates new inter-edges.
// Only replace it with another inter-edge.
void replace_inter_with_inter(Graph * graph, boost::random::mt19937 * rng, size_t num_passes, std::vector<size_t> com_assignments) {
    std::cout << "replacing inter-edges only...\n";
    size_t num_edges = graph->get_num_edges();
    size_t num_nodes = graph->get_num_nodes();
    std::vector<size_t> inter_edges(num_edges);
    size_t num_inter_edges;
    boost::random::uniform_int_distribution<size_t> all_nodes_dist(0,num_nodes-1);
    // To save time, compile a list of inter-edges once, and select edges to randomize from it.
    num_inter_edges = 0;
    for (size_t edge_index = 0; edge_index < num_edges; edge_index++) {
        inter_edges[num_inter_edges] = edge_index;
        size_t source = graph->get_source(edge_index);
        size_t target = graph->get_target(edge_index);
        num_inter_edges += (com_assignments[source] != com_assignments[target]);
    }
    boost::random::uniform_int_distribution<size_t> inter_edges_dist(0,num_inter_edges-1);
    Graph copy = graph->get_copy();
    size_t min_num_shared = num_edges;
    size_t passes_since_min = 0;
    size_t pass_num = 0;
    do {
        size_t num_swaps = 0;
        size_t num_replaces = 0;
        for (size_t edge_index_index = 0; edge_index_index < num_inter_edges; edge_index_index++) {
            size_t edge_index = inter_edges[edge_index_index];
            size_t swap_partner_index = inter_edges_dist(*rng);
            size_t swap_partner = inter_edges[swap_partner_index];
            size_t source_a = graph->get_source(edge_index);
            size_t target_a = graph->get_target(edge_index);
            size_t source_b = graph->get_source(swap_partner);
            size_t target_b = graph->get_target(swap_partner);
            if ( (com_assignments[source_a] != com_assignments[target_b]) && (com_assignments[source_b] != com_assignments[target_a]) ) {
                num_swaps += graph->swap_endpoints(edge_index, swap_partner);
            }
            size_t new_source = all_nodes_dist(*rng);
            size_t new_target = all_nodes_dist(*rng);
            if (com_assignments[new_source] != com_assignments[new_target]) {
                num_replaces += graph->replace_edge(edge_index, new_source, new_target);
            }
        }
        pass_num++;
        // Stop when we have gone num_passes passes without finding a randomization further away from the original.
        // The longer we have gone since seeing a new minimum, the more likely it is that we have reached a point of diminishing returns.
        size_t num_shared = get_num_shared( graph, &copy );
        if (num_shared < min_num_shared) {
            min_num_shared = num_shared;
            passes_since_min = 0;
        } else {
            passes_since_min++;
        }
        std::cout << "pass " << pass_num << " complete: " << num_swaps << " swaps, " << num_replaces << " replacements, " << num_shared << " edges also in original\n";
    } while (passes_since_min < num_passes);
    std::cout << "done\n";
}

// Attempt num_passes swaps and num_passes replacements for each edge with endpoints in the same community (intra-edge).
// Only swap it if the swap creates new intra-edges.
// Only replace it with another intra-edge in the same community.
void replace_intra_with_intra(Graph * graph, boost::random::mt19937 * rng, size_t num_passes, std::vector<size_t> com_assignments) {
    std::cout << "replacing intra-edges only...\n";
    const size_t num_edges = graph->get_num_edges();
    const size_t num_nodes = graph->get_num_nodes();
    size_t max_community_index = 0;
    for (size_t node = 0; node < num_nodes; node++) {
        if (max_community_index < com_assignments[node]) {
            max_community_index = com_assignments[node];
        }
    }
    size_t num_communities = max_community_index+1;
    // Compile a list of all the nodes in this community so that we can select replacement endpoints from it.
    std::vector<size_t> nodes_in_community(num_communities*num_nodes);
    std::vector<size_t> num_nodes_in_community(num_communities, 0);
    // Compile a list of all intra-edges in this community so that we can select swap partners from it.
    std::vector<size_t> intra_edges(num_communities*num_edges);
    std::vector<size_t> num_intra_edges (num_communities, 0);
    // Doing this once before all the passes saves time compared to selecting from all the nodes/edges and testing
    // or re-compiling the list on each pass.
    for (size_t community = 0; community < num_communities; community++) {
        size_t node_community_offset = num_nodes*community;
        for (size_t node = 0; node < num_nodes; node++) {
            nodes_in_community[ node_community_offset + num_nodes_in_community[community] ] = node;
            num_nodes_in_community[community] += (com_assignments[node] == community);
        }
        size_t edge_community_offset = num_edges*community;
        for (size_t edge = 0; edge < num_edges; edge++) {
            intra_edges[ edge_community_offset + num_intra_edges[community] ] = edge;
            size_t source = graph->get_source(edge);
            size_t target = graph->get_target(edge);
            num_intra_edges[community] += ( com_assignments[source] == community ) && ( com_assignments[target] == community );
        }
    }
    Graph copy = graph->get_copy();
    size_t min_num_shared = num_edges;
    size_t passes_since_min = 0;
    size_t pass_num = 0;
    do {
        size_t num_swaps = 0;
        size_t num_replaces = 0;
        for (size_t community = 0; community < num_communities; community++) {
            boost::random::uniform_int_distribution<size_t> nodes_in_community_dist(0,num_nodes_in_community[community]-1);
            boost::random::uniform_int_distribution<size_t> intra_edges_dist(0,num_intra_edges[community]-1);
            size_t node_community_offset = num_nodes*community;
            size_t edge_community_offset = num_edges*community;
            for (size_t edge_index = 0; edge_index < num_intra_edges[community]; edge_index++) {
                size_t edge = intra_edges[edge_community_offset + edge_index];
                size_t swap_partner_index = intra_edges_dist(*rng);
                size_t swap_partner = intra_edges[edge_community_offset + swap_partner_index];
                num_swaps += graph->swap_endpoints(edge, swap_partner);
                size_t new_source_index = nodes_in_community_dist(*rng);
                size_t new_source = nodes_in_community[node_community_offset + new_source_index];
                size_t new_target_index = nodes_in_community_dist(*rng);
                size_t new_target = nodes_in_community[node_community_offset + new_target_index];
                num_replaces += graph->replace_edge(edge, new_source, new_target);
            }
        }
        pass_num++;
        // Stop when we have gone num_passes passes without finding a randomization further away from the original.
        // The longer we have gone since seeing a new minimum, the more likely it is that we have reached a point of diminishing returns.
        size_t num_shared = get_num_shared( graph, &copy );
        if (num_shared < min_num_shared) {
            min_num_shared = num_shared;
            passes_since_min = 0;
        } else {
            passes_since_min++;
        }
        std::cout << "pass " << pass_num << " complete: " << num_swaps << " swaps, " << num_replaces << " replacements, " << num_shared << " edges also in original\n";
    } while (passes_since_min < num_passes);
    std::cout << "done\n";
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

bool write_edges(std::string file_name, Graph * graph) {
    bool is_valid;
    std::ofstream edge_file(file_name);
    if ( edge_file.is_open() ) {
        for (size_t edge = 0; edge < graph->get_num_edges(); edge++) {
            // std::cout << edges[edge_index].get_source() + 1 << "\t" << edges[edge_index].get_target() + 1 << "\n";
            edge_file << graph->get_source(edge) + 1 << "\t" << graph->get_target(edge) + 1 << "\n";
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
    RandType rand_type;
    bool found_rand_type = read_rand_type_param(argc, argv, "-t", rand_type);
    if (!found_rand_type) {
        std::cerr << "Argument -t (randomization type) is required.\n";
        return 1;
    }
    size_t num_passes;
    bool found_num_passes = read_size_param(argc, argv, "-n", num_passes);
    if (!found_num_passes) {
        std::cerr << "Argument -n (number of randomization passes to go without finding a new most different graph) is required.\n";
        return 1;
    }
    std::string rand_edge_file_name;
    bool found_rand_file = read_string_param(argc, argv, "-o", rand_edge_file_name);
    if (!found_rand_file) {
        std::cerr << "Argument -o (randomized graph file) is required.\n";
        return 1;
    }
    std::cout << "found params\n";
    std::cout << "original graph file: " << original_edge_file_name << "\n";
    std::cout << "randomized graph file: " << rand_edge_file_name << "\n";
    std::cout << "rand seed: " << rand_seed << "\n";
    std::cout << "num passes: " << num_passes << "\n";
    boost::random::mt19937 rng (rand_seed);
    std::vector<size_t> sources;
    std::vector<size_t> targets;
    bool edge_read_succeeded = read_edges(original_edge_file_name, sources, targets);
    if (!edge_read_succeeded) {
        return 1;
    }
    std::cout << "read in edges\n";
    Graph graph (sources, targets);
    std::cout << "initialized graph with " << graph.get_num_nodes() << " nodes and " << graph.get_num_edges() << " edges\n";
    // Only read in a community file for the randomizations that preserve communities.
    std::vector<size_t> com_assignments;
    if ( (rand_type == communities) || (rand_type == inter_edges) || (rand_type == intra_edges) ) {
        std::string community_file_name;
        bool found_community_file = read_string_param(argc, argv, "-c", community_file_name);
        if (!found_community_file) {
            std::cerr << "Argument -c (community assignment file) is required.\n";
            return 1;
        }
        std::cout << "community assignment file: " << community_file_name << "\n";
        bool com_read_succeeded = read_communities(community_file_name, com_assignments);
        std::cout << "read in community assignments\n";
        if (!com_read_succeeded) {
            return 1;
        }
    }
    switch (rand_type) {
        case mean_degree:
            replace_all_with_any(&graph, &rng, num_passes);
            break;
        case all_degrees:
            swap_all_with_all(&graph, &rng, num_passes);
            break;
        case communities:
            replace_inter_with_inter(&graph, &rng, num_passes, com_assignments);
            replace_intra_with_intra(&graph, &rng, num_passes, com_assignments);
            break;
        case inter_edges:
            // Preserve inter-edges, randomize intra-edges.
            replace_intra_with_intra(&graph, &rng, num_passes, com_assignments);
            break;
        case intra_edges:
            // Preserve intra-edges, randomize inter-edges.
            replace_inter_with_inter(&graph, &rng, num_passes, com_assignments);
            break;
        default:
            std::cerr << "There is a bug in the code. The switch statement in main should cover all randomization types.\n";
            break;
    }
    bool write_succeeded = write_edges(rand_edge_file_name, &graph);
    if (!write_succeeded) {
        return 1;
    }
    std::cout << "wrote randomized graph\n";
    return 0;
}
