//
//  main.cpp
//  mbrw4cluster
//
//  Created by Adam Craig on 2/23/19.
//  Copyright © 2019 Drexel Systems Immunology Laboratory. All rights reserved.
//

#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <random>
#include <cmath>
#include <chrono>
#include <ctime>
// #include <boost/filesystem.hpp> -does not compile
#include <boost/lexical_cast.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

// Set missing values to this value.
// Since size_t is usually an unsigned value,
// this should wrap around to SIZE_T_MAX.
// Whether signed or unsigned, this value should not occur
// unless we explicitly set it.
// The two cases where we use it:
// 1. step_last_visited[i] when i has not been visited.
// We check for this value explicitly when
// checking whether an edge is remembered.
// If we actually end up running for SIZE_T_MAX or more steps,
// the roll-over in step count will affect our ability
// to remember which vertices we have visited anyway.
// 2. edge_states[i].back_edge when we do not find a back-edge.
// We should have fewer than SIZE_T_MAX edges in the graph
// so that this is not the index of an actual edge.
// When the forbidden back-edge is MBRW_NONE,
// none of the actual edges will be seen as forbidden.
#ifndef MBRW_NONE
#define MBRW_NONE 0-1
#endif

typedef struct {
    bool is_valid;
    const char *input_file;
    const char *output_file;
    size_t bias;
    size_t memory_length;
    size_t start_vertex;
    size_t max_doublings;
    int rand_seed;
} run_params_t;

run_params_t parse_args(const int argc, const char *argv[]) {
    run_params_t params;
    bool expect_input_file = false;
    bool expect_output_file = false;
    bool expect_bias = false;
    bool expect_memory_length = false;
    bool expect_max_doublings = false;
    bool expect_start_vertex = false;
    bool expect_rand_seed = false;
    bool found_input_file = false;
    bool found_output_file = false;
    bool found_bias = false;
    bool found_memory_length = false;
    bool found_start_vertex = false;
    bool found_max_doublings = false;
    bool found_rand_seed = false;
    // Expect one agrument after each flag.
    // e.g. ./mbrw4cluster -i ~/real_data/c_elegans_ppi/ -b 1000 -m 5 -s 10
    for (size_t arg_index = 0; arg_index < argc; arg_index++) {
        const char * arg = argv[arg_index];
        // std::cout << "arg " << arg_index << ": " << arg << "\n";
        if (expect_input_file) {
            params.input_file = arg;
            found_input_file = true;
            // Does not compile:
            // boost::filesystem::path my_file( params.input_file );
            // found_input_file = boost::filesystem::is_regular_file( my_file );
            // if (!found_input_file) {
            //     std::cout << "Error: File: " << params.input_file << " does not exist.\n";
            // }
            // else {
            // std::cout << "found input file: " << params.input_file << "\n";
            // }
            //
            expect_input_file = false;
        }
        if (expect_output_file) {
            params.output_file = arg;
            found_output_file = true;
            expect_output_file = false;
        }
        else if (expect_bias) {
            params.bias = boost::lexical_cast<size_t>(arg);
            // std::cout << "found bias: " << params.bias << "\n";
            found_bias = true;
            expect_bias = false;
        }
        else if (expect_memory_length) {
            params.memory_length = boost::lexical_cast<size_t>(arg);
            // std::cout << "found memory length: " << params.memory_length << "\n";
            found_memory_length = true;
            expect_memory_length = false;
        }
        else if (expect_start_vertex) {
            params.start_vertex = boost::lexical_cast<size_t>(arg);
            // std::cout << "found start vertex: " << params.start_vertex << "\n";
            found_start_vertex = true;
            expect_start_vertex = false;
        }
        else if (expect_max_doublings) {
            params.max_doublings = boost::lexical_cast<size_t>(arg);
            // std::cout << "found num circulations: " << params.num_circulations << "\n";
            found_max_doublings = true;
            expect_max_doublings = false;
        }
        else if (expect_rand_seed) {
            params.rand_seed = boost::lexical_cast<int>(arg);
            // std::cout << "found rand seed: " << params.rand_seed << "\n";
            found_rand_seed = true;
            expect_rand_seed = false;
        }
        else if ( !strcmp(arg, "-i") ) {
            expect_input_file = true;
        }
        else if ( !strcmp(arg, "-o") ) {
            expect_output_file = true;
        }
        else if ( !strcmp(arg, "-b") ) {
            expect_bias = true;
        }
        else if ( !strcmp(arg, "-m") ) {
            expect_memory_length = true;
        }
        else if ( !strcmp(arg, "-s") ) {
            expect_start_vertex = true;
        }
        else if ( !strcmp(arg, "-c") ) {
            expect_max_doublings = true;
        }
        else if ( !strcmp(arg, "-r") ) {
            expect_rand_seed = true;
        }
        // else {
        //     std::cout << "skipping arg " << arg_index << ", " << arg << "\n";
        // }
    }
    // We need to have a directory, bias, and memory length.
    params.is_valid = found_bias && found_memory_length;
    // If no input directory is provided, assume the current directory.
    if (!found_input_file) {
        params.input_file = "edges.txt";
    }
    // If no output directory is provided, assume the current directory.
    if (!found_output_file) {
        params.output_file = "entropies.csv";
    }
    // If no starting vertex is provided, use MBRW_NONE to tell the constructor to randomly select one.
    if (!found_start_vertex) {
        params.start_vertex = MBRW_NONE;
    }
    // If no max fractional change is provided, use 31 => stop after 2^31 = around 2 billion steps,
    // takes around 2 hours and seems to be big enough for the networks we are using when bias = memory = 100k.
    if (!found_max_doublings) {
        params.max_doublings = 31;
    }
    // If no rand seed is provided, use a number from the bad RNG to initialize the good one.
    if (!found_rand_seed) {
        params.rand_seed = std::rand();
    }
    return params;
}

// some edge in the graph A->B
class Edge {
    
    // ID of A
    std::string source_id;
    // ID of B
    std::string target_id;
    
public:
    
    Edge();
    Edge(std::string source, std::string target);
    void set_source(std::string source);
    void set_target(std::string target);
    void set_nodes(std::string source, std::string target);
    std::string get_source();
    std::string get_target();
    
};

Edge::Edge() {
}

Edge::Edge(std::string source, std::string target) {
    source_id = source;
    target_id = target;
}

void Edge::set_source(std::string source) {
    source_id = source;
}

void Edge::set_target(std::string target) {
    target_id = target;
}

void Edge::set_nodes(std::string source, std::string target) {
    source_id = source;
    target_id = target;
}

std::string Edge::get_source() {
    return source_id;
}

std::string Edge::get_target() {
    return target_id;
}

// We assume that the edges are undirected so each connected pair is only listed once.
bool read_edges(std::string file_name, std::vector<Edge> & edges) {
    bool is_valid;
    const std::string delimiter = "\t";
    size_t num_edges = 0;
    std::ifstream edge_file;
    std::string line;
    edge_file.open(file_name);
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
            if ( line.find(delimiter) == line.npos ) {
                std::cerr << "Error: edge file " << file_name << " contained misformatted line " << line << "\n";
                is_valid = false;
                break;
            }
            num_edges += 2;
        }
        // std::cout << "found " << num_edges << " edges\n";
        if (is_valid) {
            // std::cout << "reading to edge list\n";
            // Pre-allocate a vector of appropriate size.
            edges.resize(num_edges);
            // Reset to the beginning of the file.
            edge_file.clear();
            edge_file.seekg(0);
            // Read in the results,
            // this time saving the results to the edge list.
            size_t edge_index = 0;
            while ( std::getline(edge_file, line) ) {
                // std::cout << line << "\n";
                // Find the first occurrence of the tab delimiter.
                size_t tab_pos = line.find("\t");
                // std::cout << line.substr(0, tab_pos) << "\t" << line.substr(tab_pos+1) << "\n";
                // Add edges in both directions.
                std::string node1 = line.substr(0, tab_pos);
                std::string node2 = line.substr(tab_pos+1);
                if ( node2[ node2.length()-1 ] == '\r' ) {
                    // std::cout << "found carriage return on line " << edge_index << "\n";
                    node2 = node2.substr( 0, node2.length()-1 );
                }
                edges[edge_index].set_nodes( node1, node2 );
                edge_index++;
                edges[edge_index].set_nodes( node2 , node1 );
                edge_index++;
            }
        }
        edge_file.close();
    }
    else {
        std::cerr << "Error: failed to open file " << file_name << "\n";
        is_valid = false;
    }
    return is_valid;
}

class MBRWAgent {
    struct node_state_t {
        std::string id;
        // All edges leading away from this node
        // represented by their indices in edge_states.
        std::vector<size_t> out_edges;
    };
    struct edge_state_t {
        // indices of the nodes in node_states
        size_t source;
        size_t target;
        // index of the reverse edge B->A
        size_t back_edge;
        // step at which the agent last visited this edge
        size_t step_last_taken;
    };
    // states of all the edges in the graph
    std::vector<struct edge_state_t> edge_states;
    // states of all nodes in the graph
    std::vector<struct node_state_t> node_states;
    // number of steps back remembered
    size_t memory_length;
    // weight given to remembered edges - 1
    size_t bias;
    // Keep track of the step number.
    size_t steps_taken;
    // Keep track of the last edge taken.
    size_t last_edge;
    // The random number generator.
    boost::random::mt19937 generator;
    // Size these vectors so that we can re-use them at each step.
    // All edges leading away from the current node
    // that are allowed but not remembered.
    std::vector<size_t> allowed_unremembered_edges;
    // All edges leading away from the current node
    // that are allowed and remembered.
    std::vector<size_t> allowed_remembered_edges;
    // Keep track of whether the last edge taken was a remembered edge.
    bool took_remembered;
    // Keep track of when the last edge taken was taken before that.
    size_t last_step_last_taken;
    // Whether to let the agent take the reverse of the previous transition.
    bool allow_backtracking;
public:
    // Create the agent state from the
    MBRWAgent(std::vector<Edge> edges, size_t new_bias, size_t new_memory_length, size_t start_vertex, int rand_seed);
    void take_step();
    size_t get_last_edge_index();
    size_t get_last_step_last_taken();
    size_t get_num_nodes();
    size_t get_previous_node();
    size_t get_current_node();
    bool get_took_remembered();
    size_t get_steps_taken();
};

MBRWAgent::MBRWAgent (std::vector<Edge> edges, size_t new_bias, size_t new_memory_length, size_t start_vertex, int rand_seed) {
    // Set the bias and memory length.
    bias = new_bias;
    memory_length = new_memory_length;
    // std::cout << "initializing MBRW agent\n";
    // Seed the generator.
    generator = boost::random::mt19937 (rand_seed);
    // std::cout << "initialized RNG\n";
    // Create a state for each edge.
    edge_states.resize( edges.size() );
    // std::cout << "resized edge states array to " << edge_states.size() << "\n";
    // Create a map from the node ID
    // to what will be the index of the node.
    std::unordered_map<std::string, size_t> id_to_index;
    size_t num_nodes = 0;
    for (size_t edge_i = 0; edge_i < edges.size(); edge_i++) {
        // Get the index for each node.
        // If it does not yet have an index,
        // assign it the next index in order.
        std::string source = edges[edge_i].get_source();
        size_t source_i;
        if ( id_to_index.count(source) ) {
            source_i = id_to_index[source];
            // std::cout << "already have node " << source << "\n";
        }
        else {
            source_i = num_nodes;
            id_to_index[source] = source_i;
            num_nodes++;
            // std::cout << "added node " << source << "\n";
        }
        std::string target = edges[edge_i].get_target();
        size_t target_i;
        if ( id_to_index.count(target) ) {
            target_i = id_to_index[target];
            // std::cout << "already have node " << target << "\n";
        } else {
            target_i = num_nodes;
            id_to_index[target] = target_i;
            num_nodes++;
            // std::cout << "added node " << target << "\n";
        }
        edge_states[edge_i].source = source_i;
        edge_states[edge_i].target = target_i;
        // The edge has not been visited.
        // As a convention, we have an artificial step 0 where every edge was taken all at once.
        // We then consider the first step and the start of the first circulation step 1.
        // This lets us count first visits properly during the first circulation
        // and include the time interval up to the first visit to each edge in the entropy calculation.
        edge_states[edge_i].step_last_taken = 0;
        // The back-edge is unset.
        edge_states[edge_i].back_edge = MBRW_NONE;
    }
    // std::cout << "set up edge states\n";
    // Create a state for each node.
    node_states.resize(num_nodes);
    // std::cout << "resized node states to " << node_states.size() << "\n";
    for (auto it = id_to_index.begin(); it != id_to_index.end(); ++it) {
        // Record the ID of each node.
        node_states[it->second].id = it->first;
    }
    // std::cout << "assigned IDs to nodes\n";
    // Add the out-edges to each node state.
    for (size_t edge_i = 0; edge_i < edge_states.size(); edge_i++) {
        size_t source_i = edge_states[edge_i].source;
        node_states[source_i].out_edges.push_back(edge_i);
    }
    // std::cout << "added out-edges to node states\n";
    // Find the maximum degree.
    size_t max_out_degree = 0;
    for (size_t node_i = 0; node_i < node_states.size(); node_i++) {
        size_t out_degree = node_states[node_i].out_edges.size();
        max_out_degree = std::max(out_degree, max_out_degree);
        if ( node_states[node_i].out_edges.size() < 2 ) {
            std::cout << "Node " << node_i << " has degreee < 2. This will end badly.\n";
        }
    }
    // std::cout << "found max out degree " << max_out_degree << "\n";
    // Resize the out-edge working spaces to accommodate this many.
    allowed_unremembered_edges.resize(max_out_degree);
    allowed_remembered_edges.resize(max_out_degree);
    // std::cout << "resized allowed edge workspaces\n";
    // If back-tracking is forbidden, find the back-edge for each edge.
    // std::cout << "took first step\n";
    allow_backtracking = memory_length == 0;
    if (!allow_backtracking) {
        // Find the back-edge for each edge.
        // For each edge A->B,
        for (size_t edge_i = 0; edge_i < edge_states.size(); edge_i++) {
            size_t source_i = edge_states[edge_i].source;
            size_t target_i = edge_states[edge_i].target;
            // search all edges B->C
            for (size_t back_i_i = 0;
                 back_i_i < node_states[target_i].out_edges.size();
                 back_i_i++) {
                size_t back_i = node_states[target_i].out_edges[back_i_i];
                // The back-edge is the edge such that C == A.
                if (edge_states[back_i].target == source_i) {
                    edge_states[edge_i].back_edge = back_i;
                    // Each edge must be unique,
                    // so there can only be one back-edge.
                    break;
                }
            }
        }
    }
    // std::cout << "assigned back-edges\n";
    // If the start_vertex is MBRW_NONE, randomly select one.
    if (start_vertex == MBRW_NONE) {
        // Randomly select a starting node.
        boost::random::uniform_int_distribution<size_t> dist_all( 0, node_states.size() - 1 );
        start_vertex = dist_all(generator);
    }
    // Randomly select one of its edges to take.
    boost::random::uniform_int_distribution<size_t> dist_targets( 0, node_states[start_vertex].out_edges.size() - 1 );
    last_edge = node_states[start_vertex].out_edges[ dist_targets(generator) ];
    // The edge was not taken before.
    last_step_last_taken = 0;
    // Count this as the first step.
    steps_taken = 1;
    // Record that we took it on the first step.
    edge_states[last_edge].step_last_taken = steps_taken;
}

void MBRWAgent::take_step() {
    // The source of the next edge is the target of the previous one.
    size_t source = edge_states[last_edge].target;
    // std::cout << "starting at " << source << "\n";
    // The forbidden edge is the one that leads back to the source of the last edge,
    // unless we allow backtracking.
    size_t forbidden_edge = edge_states[last_edge].back_edge;
    // Accumulate allowed edges and remembered allowed edges.
    // std::cout << "choices:\n";
    size_t num_unremembered = 0;
    size_t num_remembered = 0;
    for (size_t edge_i_i = 0; edge_i_i < node_states[source].out_edges.size(); edge_i_i++) {
        size_t edge_i = node_states[source].out_edges[edge_i_i];
        // Do it this way to avoid using if statements/branch prediction.
        bool is_allowed = edge_i != forbidden_edge;
        size_t step_last_taken = edge_states[edge_i].step_last_taken;
        bool is_remembered = (step_last_taken != 0) && (steps_taken - step_last_taken < memory_length);
        allowed_remembered_edges[num_remembered] = edge_i;
        num_remembered += (is_allowed && is_remembered);
        allowed_unremembered_edges[num_unremembered] = edge_i;
        num_unremembered += (is_allowed && !is_remembered);
        // std::cout << "\t" << edge_i << "\t";
        // if (edge_i != forbidden_edge) {
        // Examples:
        // 3 - 1 > 3 -> true
        // 4 < 3 + 1 -> false
        // if ( (step == 0) || (steps_taken - step >= memory_length) ) {
        //     allowed_unremembered_edges[num_unremembered] = edge_i;
        //     num_unremembered++;
        //     // std::cout << "unremembered\n";
        // }
        // else {
        //     allowed_remembered_edges[num_remembered] = edge_i;
        //     num_remembered++;
        //     // std::cout << "remembered\n";
        // }
        // }
        // else {
        //     std::cout << "forbidden\n";
        // }
    }
    // Set the range of values for the random number.
    // Imagine the space of numbers from 1 to num_unremembered+bias*num_remembered
    // as divided into bins.
    // The first num_unremembered bins are of size 1 and belong to the unremembered edges.
    // The subsequent bins are of size bias, and each belongs to a remembered edge.
    // Note: The range for uniform_int_distribution is inclusive at both endpoints.
    // Example: bias=1k, 5 unremembered and 2 remembered -> range [0, 2004]
    // Generate 0 < num_unremembered, so we take unremembered_edges[0].
    // Generate 4 < num_unremembered, so we take unremembered_edges[4].
    // Generate 5 == num_unremembered and (5 - 5)/1000 = 0, so remembered_edges[0].
    // Generate 1004 > num_unremembered and (1004 - 5)/1000 = 0, so remembered_edges[0].
    // Generate 1005 > num_unremembered and (1005 - 5)/1000 = 1, so remembered_edges[1].
    // Generate 2004 > num_unremembered and (2004 - 5)/1000 = 1, so remembered_edges[1].
    boost::random::uniform_int_distribution<size_t> dist(0, num_unremembered + bias*num_remembered - 1);
    size_t rand_num = dist(generator);
    // std::cout << "random number: " << rand_num << "\n";
    took_remembered = rand_num >= num_unremembered;
    if (took_remembered) {
        last_edge = allowed_remembered_edges[(rand_num - num_unremembered)/bias];
        // std::cout << "chose remembered: ";
    } else {
        last_edge = allowed_unremembered_edges[rand_num];
        // std::cout << "chose unremembered: ";
    }
    // std::cout << last_edge << ":"
    // << edge_states[last_edge].source << "->" << edge_states[last_edge].target << "\n";
    // Record when we last visited this edge.
    last_step_last_taken = edge_states[last_edge].step_last_taken;
    // Update the count of steps and the step on which we took this edge.
    steps_taken++;
    edge_states[last_edge].step_last_taken = steps_taken;
}

size_t MBRWAgent::get_last_edge_index() {
    return last_edge;
}

size_t MBRWAgent::get_last_step_last_taken() {
    return last_step_last_taken;
}

size_t MBRWAgent::get_num_nodes() {
    return node_states.size();
}

size_t MBRWAgent::get_previous_node() {
    return edge_states[last_edge].source;
}

size_t MBRWAgent::get_current_node() {
    return edge_states[last_edge].target;
}

bool MBRWAgent::get_took_remembered() {
    return took_remembered;
}

size_t MBRWAgent::get_steps_taken() {
    return steps_taken;
}

typedef struct {
    std::string source;
    std::string target;
    double entropy;
    size_t num_intervals;
    double normed_entropy;
    size_t num_coms;
    size_t num_singletons;
    size_t max_size;
    double modularity;
} EdgeEntropy;

std::vector<EdgeEntropy> get_sorted_entropies(std::vector<Edge> edges,
                                              std::vector<double> & interval_sums, std::vector<double> & entropy_sums,
                                              std::vector<size_t> & num_intervals, size_t num_steps) {
    std::vector<EdgeEntropy> edge_entropies ( edges.size() );
    double log_num_steps = log2( (double) num_steps );
    for (size_t edge_index = 0; edge_index < interval_sums.size(); edge_index++) {
        // Calculate the entropy of time intervals between visits to some node.
        // Let the intervals be x1, x2, ..., xn and the total number of steps be N.
        // entropy = -sum( (xi/N)*log(xi/N) )
        //         = -sum(  (xi/N)*( log(xi)-log(N) )  )
        //         = -sum( (xi/N)*log(xi)-(xi/N)*log(N) )
        //         = -sum( (xi/N)*log(xi) ) + -sum( -(xi/N)*log(N) )
        //         = -sum( (xi/N)*log(xi) ) + sum( (xi/N)*log(N) )
        //         = -sum( xi*log(xi) )/N + ( sum(xi)/N )*log(N)
        //         = -sum( xi*log(xi) )/N + sum(xi)*log(N)/N
        //         = sum(xi)*log(N)/N - sum( xi*log(xi) )/N
        //         = (  sum(xi)*log(N) - sum( xi*log(xi) )  )/N
        // Store sum(xi) in interval_sums, xi*log(xi) in entropy_sums, N in num_steps.
        edge_entropies[edge_index].source = edges[edge_index].get_source();
        edge_entropies[edge_index].target = edges[edge_index].get_target();
        edge_entropies[edge_index].entropy = ( interval_sums[edge_index]*log_num_steps - entropy_sums[edge_index] )/( (double) num_steps );
        // Normalize the entropies so that we can compare them.
        // The entropy is at its maximum when all intervals are of equal size,
        // that is, when visits are the least clumpy,
        // or, when a randomly selected step is equally likely to be in any interval.
        // In this case, if there are n intervals and N steps, every xi = N/n, so xi/N = 1/n.
        // entropy = -sum( (xi/N)*log(xi/N) ) = -sum( (1/n)*log(1/n) ) = -n*(1/n)*log(1/n) = -log(1/n) = log(n)
        edge_entropies[edge_index].num_intervals = num_intervals[edge_index];
        edge_entropies[edge_index].normed_entropy = edge_entropies[edge_index].entropy / log2( (double) edge_entropies[edge_index].num_intervals );
    }
    // Use insertion sort to put the transitions in order of increasing entropy => greater betweenness/lower embeddedness in communities.
    for (size_t edge_index = 1; edge_index < interval_sums.size(); edge_index++) {
        size_t new_index = edge_index;
        while ( (new_index > 0) && (edge_entropies[new_index].normed_entropy < edge_entropies[new_index-1].normed_entropy) ) {
            EdgeEntropy temp = edge_entropies[new_index];
            edge_entropies[new_index] = edge_entropies[new_index-1];
            edge_entropies[new_index-1] = temp;
            new_index--;
        }
    }
    return edge_entropies;
}

typedef struct {
    size_t degree;
    size_t com_index;
} NodeComInfo;

typedef std::unordered_map<std::string, NodeComInfo> CommunityAssignments;

// Use the formula from
// https://en.wikipedia.org/wiki/Louvain_Modularity
// with all edge weights set to 1 if the edge is in the graph or 0 if it is not.
double get_modularity( std::vector<EdgeEntropy> edge_entropies, CommunityAssignments community_assignments ) {
    // This is already 2m, since we count both directions of each undirected edge as separate edges here.
    size_t num_edges = edge_entropies.size();
    // Count the number of intra-community edges.
    // This may be more than the number of edges that have entropies <= the threshold,
    // which are included as intra-edges by our definition of communities.
    size_t num_intra_edges = 0;
    for (size_t edge_index = 0; edge_index < edge_entropies.size(); edge_index++) {
        size_t source_com_index = community_assignments[ edge_entropies[edge_index].source ].com_index;
        size_t target_com_index = community_assignments[ edge_entropies[edge_index].target ].com_index;
        num_intra_edges += (source_com_index == target_com_index);
    }
    size_t degree_product_sum = 0;
    for (auto it1 = community_assignments.begin(); it1 != community_assignments.end(); ++it1) {
        size_t com_index_1 = (it1->second).com_index;
        size_t degree_1 = (it1->second).degree;
        for (auto it2 = community_assignments.begin(); it2 != community_assignments.end(); ++it2) {
            degree_product_sum += ( com_index_1 == (it2->second).com_index )
            * degree_1 * (it2->second).degree;
        }
    }
    return ( (double) num_intra_edges )/( (double) num_edges )
    - ( (double) degree_product_sum )/( (double) num_edges*num_edges );
}

typedef struct {
    size_t community;
    size_t size;
} ComSizePair;

void standardize_community_nums(CommunityAssignments & community_assignments, std::vector<size_t> community_sizes) {
    // First, sort the communities by size from largest to smallest.
    // Use insertion sort.
    size_t num_coms = community_sizes.size();
    std::vector<ComSizePair> com_size_pairs ( num_coms );
    com_size_pairs[0].community = 0;
    com_size_pairs[0].size = community_sizes[0];
    for (size_t com_index = 1; com_index < num_coms; com_index++) {
        ComSizePair new_pair;
        new_pair.community = com_index;
        new_pair.size = community_sizes[com_index];
        size_t pair_index = com_index;
        while ( (pair_index > 0) && (new_pair.size > com_size_pairs[pair_index-1].size) ) {
            com_size_pairs[pair_index] = com_size_pairs[pair_index-1];
            pair_index--;
        }
        com_size_pairs[pair_index] = new_pair;
    }
    // Count the number of non-empty communities.
    size_t num_nonempty = 0;
    while ( (num_nonempty < num_coms) && (com_size_pairs[num_nonempty].size > 0) ) {
        num_nonempty++;
    }
    std::cout << "# communities: " << num_nonempty << "\n";
    // Now assign each non-zero community a new number based on its rank.
    // Make the numbers start at 1 to follow the convention of the Infomap communities.
    std::vector<size_t> ranks ( num_coms, 0 );
    for (size_t pair_index = 0; pair_index < num_nonempty; pair_index++) {
        ranks[ com_size_pairs[pair_index].community ] = pair_index + 1;
        std::cout << "# size(" << ranks[ com_size_pairs[pair_index].community ] << "): " << com_size_pairs[pair_index].size << "\n";
        // std::cout << "com " << com_size_pairs[pair_index].community << ", size " << com_size_pairs[pair_index].size << ", rank " << ranks[ com_size_pairs[pair_index].community ] << "\n";
    }
    // Change each node's community assignment to reflect the new numbering.
    for (auto it = community_assignments.begin(); it != community_assignments.end(); ++it) {
        // std::cout << "node " << it->first << ", com " << (it->second).com_index;
        (it->second).com_index = ranks[(it->second).com_index];
        // std::cout << "->" << (it->second).com_index << "\n";
    }
}

// To make communities with entropies:
// 1. Start with the original graph.
// 2. Remove all the edges where the normed entropies in both directions are above a certain threshold.
//    (That is, keep only the edges with normed entropy at or below a certain threshold in at least one direction.)
// 3. The connected components are the communities.
// Since we have already sorted the transitions in order of increasing entropy,
// we can construct the partitionings for all possible thresholds as follows:
// 1. Start with each node assigned to a singleton community.
// 2. For each transition, R, starting from the lowest entropy to the highest,
//    merge the two communities connected by the edge by
//    re-assigning all nodes in the community of the target to the community of the source.
//    This is the partition for thresholds, H, such that
// edge_entropies[R].normed_entropy <= H < edge_entropies[R+1].normed_entropy
void get_communities(std::vector<EdgeEntropy> & edge_entropies) {
    // First, just make a set of all the unique nodes.
    std::unordered_set< std::string > nodes;
    for (size_t edge_index = 0; edge_index < edge_entropies.size(); edge_index++) {
        nodes.insert( edge_entropies[edge_index].source );
        nodes.insert( edge_entropies[edge_index].target );
    }
    // Make a list of the community assignments of all the nodes.
    // Start with each node in a singleton community.
    CommunityAssignments community_assignments ( nodes.size() );
    // Also keep track of the size of each community.
    std::vector<size_t> community_sizes( nodes.size() );
    size_t com_index = 0;
    for (auto it = nodes.begin(); it != nodes.end(); ++it) {
        community_assignments[ *it ].com_index = com_index;
        // Initialize degree to 0 so we can increment it in the next step.
        community_assignments[ *it ].degree = 0;
        // Initialize the size of each community to 1.
        community_sizes[com_index] = 1;
        com_index++;
    }
    // Compute the degree of each node.
    for (size_t edge_index = 0; edge_index < edge_entropies.size(); edge_index++) {
        // We increment the degree of the target when we get to the directed edge in the opposite direction.
        community_assignments[ edge_entropies[edge_index].source ].degree++;
    }
    // We start with as many communities as we have nodes.
    size_t num_coms = community_assignments.size();
    // If no edges are intra-edges, then modularity is 0.
    double modularity = 0.0;
    // Initially, all are singletons.
    size_t num_singletons = nodes.size();
    // Thus, the maximum size is 1.
    size_t max_size = 1;
    // Keep track of the partitioning that gives us the best modularity.
    double max_modularity = 0.0;
    CommunityAssignments best_coms;
    std::vector<size_t> best_coms_sizes;
    // Go through the edges from lowest to highest normed entropy until everything is grouped into one community.
    for (size_t edge_index = 0; edge_index < edge_entropies.size(); edge_index++) {
        size_t source_com_index = community_assignments[ edge_entropies[edge_index].source ].com_index;
        size_t target_com_index = community_assignments[ edge_entropies[edge_index].target ].com_index;
        // The communities to which these nodes belong may already have been merged via another edge.
        // If so, then we do not need to update the number of communities of modularity.
        // If not, merge them now.
        if (target_com_index != source_com_index) {
            // Merge the communities connected by this edge.
            for (auto it = community_assignments.begin(); it != community_assignments.end(); ++it) {
                if ( (it->second).com_index == target_com_index ) {
                    (it->second).com_index = source_com_index;
                }
            }
            // std::cout << "intra-edge: " << edge_entropies[edge_index].source << "-" << edge_index << "->" << edge_entropies[edge_index].target << "\n";
            // for (auto it = community_assignments.begin(); it != community_assignments.end(); ++it) {
            //     std::cout << (it->first) << ":\t" << (it->second).com_index << "\n";
            // }
            // We no longer have the target community as a distinct community.
            num_coms--;
            // If either one was a singleton, it no longer is.
            num_singletons -= (community_sizes[source_com_index] == 1);
            num_singletons -= (community_sizes[target_com_index] == 1);
            // Adjust the community sizes to reflect the change.
            community_sizes[source_com_index] += community_sizes[target_com_index];
            community_sizes[target_com_index] = 0;
            // The source community is the only one to increase in size this step,
            // so it is the only candidate for a new maximum size.
            if (community_sizes[source_com_index] > max_size) {
                max_size = community_sizes[source_com_index];
            }
            // Recompute modularity for the new community assignments.
            modularity = get_modularity(edge_entropies, community_assignments);
            // If we find a new best partitioning, save a copy of it.
            if (modularity > max_modularity) {
                max_modularity = modularity;
                best_coms = CommunityAssignments (community_assignments);
                best_coms_sizes = std::vector<size_t> (community_sizes);
            }
        }
        // Record the number of communities remaining and the modularity of this partitioning.
        edge_entropies[edge_index].num_coms = num_coms;
        edge_entropies[edge_index].num_singletons = num_singletons;
        edge_entropies[edge_index].max_size = max_size;
        edge_entropies[edge_index].modularity = modularity;
    }
    // std::cout << "# max modularity: " << max_modularity << "\n";
    // for (auto it = best_coms.begin(); it != best_coms.end(); ++it) {
    //    std::cout << (it->first) << " " << (it->second).com_index << "\n";
    // }
    // Fix the numbering of the best community partitioning.
    standardize_community_nums(best_coms, best_coms_sizes);
    // Print the best community partitioning.
    std::cout << "# modularity: " << max_modularity << "\n";
    for (auto it = best_coms.begin(); it != best_coms.end(); ++it) {
        std::cout << (it->first) << " " << (it->second).com_index << "\n";
    }
}

bool write_entropies(std::string file_name, std::vector<EdgeEntropy> edge_entropies) {
    bool is_valid;
    std::ofstream entropy_file(file_name);
    if ( entropy_file.is_open() ) {
        // entropy_file << "source,target,entropy,intervals,normed,communities,singletons,maxsize,modularity\n";
        entropy_file << "source,target,entropy,intervals,normed\n";
        for (size_t edge_index = 0; edge_index < edge_entropies.size(); edge_index++) {
            // entropy = (  sum(xi)*log(N) - sum( xi*log(xi) )  )/N
            entropy_file << edge_entropies[edge_index].source << ","
            << edge_entropies[edge_index].target << ","
            << edge_entropies[edge_index].entropy << ","
            << edge_entropies[edge_index].num_intervals << ","
            << edge_entropies[edge_index].normed_entropy << "\n";// ","
            // << edge_entropies[edge_index].num_coms << ","
            // << edge_entropies[edge_index].num_singletons << ","
            // << edge_entropies[edge_index].max_size << ","
            // << edge_entropies[edge_index].modularity << "\n";
        }
        if ( entropy_file.fail() ) {
            std::cout << strerror(errno);
            is_valid = false;
        }
        else {
            is_valid = true;
        }
        entropy_file.close();
    }
    else {
        std::cout << "Error: failed to open file " << file_name << "\n";
        is_valid = false;
    }
    return is_valid;
}

// Calculate the entropy of time intervals between visits to some node.
// Let the intervals be x1, x2, ..., xn and the total number of steps be N.
// entropy = -sum( (xi/N)*log(xi/N) )
//         = -sum(  (xi/N)*( log(xi)-log(N) )  )
//         = -sum( (xi/N)*log(xi)-(xi/N)*log(N) )
//         = -sum( (xi/N)*log(xi) ) + -sum( -(xi/N)*log(N) )
//         = -sum( (xi/N)*log(xi) ) + sum( (xi/N)*log(N) )
//         = -sum( xi*log(xi) )/N + ( sum(xi)/N )*log(N)
//         = -sum( xi*log(xi) )/N + sum(xi)*log(N)/N
//         = sum(xi)*log(N)/N - sum( xi*log(xi) )/N
//         = (  sum(xi)*log(N) - sum( xi*log(xi) )  )/N
// Store sum(xi) in interval_sums, xi*log(xi) in entropy_sums, N in num_steps.
// Normalize the entropies so that we can compare them.
// The entropy is at its maximum when all intervals are of equal size,
// that is, when visits are the least clumpy,
// or, when a randomly selected step is equally likely to be in any interval.
// In this case, if there are n intervals and N steps, every xi = N/n, so xi/N = 1/n.
// entropy = -sum( (xi/N)*log(xi/N) ) = -sum( (1/n)*log(1/n) ) = -n*(1/n)*log(1/n) = -log(1/n) = log(n)
// Thus, normed entropy = (  sum(xi)*log(N) - sum( xi*log(xi) )  )/( N*log(n) )
std::vector<double> get_normed_entropies(std::vector<double> & interval_sums,
                                         std::vector<double> & entropy_sums,
                                         std::vector<size_t> & num_intervals,
                                         size_t num_steps) {
    std::vector<double> normed_entropies ( interval_sums.size() );
    double log_num_steps = log2( (double) num_steps );
    for (size_t edge_index = 0; edge_index < interval_sums.size(); edge_index++) {
        normed_entropies[edge_index] = ( interval_sums[edge_index]*log_num_steps - entropy_sums[edge_index] )
        / (  (double) num_steps * log2( (double) num_intervals[edge_index] )  );
    }
    return normed_entropies;
}

// Output a matrix M as a CSV file so that M(i,j) is sum( masses(where length=2^i).^orders(j) )
bool write_matrix_to_csv(std::string file_name, std::vector<double> matrix_data, size_t num_rows, size_t num_cols) {
    bool is_valid;
    std::ofstream csv_file(file_name);
    if ( csv_file.is_open() ) {
        // The last segment mass is for a length with no complete segments and so is not valid.
        size_t last_row_offset = (num_rows-1)*num_cols;
        size_t num_cols_minus_1 = num_cols-1;
        for (size_t row_offset = 0; row_offset <= last_row_offset; row_offset += num_cols) {
            // entropy = (  sum(xi)*log(N) - sum( xi*log(xi) )  )/N
            for (size_t col_index = 0; col_index < num_cols; col_index++) {
                csv_file << matrix_data[row_offset + col_index];
                if (col_index < num_cols_minus_1) {
                    csv_file << ",";
                }
                else {
                    csv_file << "\n";
                }
            }
        }
        if ( csv_file.fail() ) {
            std::cerr << strerror(errno);
            is_valid = false;
        }
        else {
            is_valid = true;
        }
        csv_file.close();
    }
    else {
        std::cerr << "Error: failed to open file " << file_name << "\n";
        is_valid = false;
    }
    return is_valid;
}

int main(int argc, const char * argv[]) {
    
    auto start = std::chrono::system_clock::now();
    
    std::vector<double> orders = { -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
    const size_t num_orders = orders.size();
    
    // std::cout << "MBRW for clusters\n";
    run_params_t params = parse_args(argc, argv);
    if (!params.is_valid) {
        std::cout << "Reading comand line args failed.\n";
        return 1;
    }
    std::cout << "MBRW entropy:\n"
    << "input file: " << params.input_file << "\n"
    << "output file: " << params.output_file << "\n"
    << "bias: " << params.bias << "\n"
    << "memory length: " << params.memory_length << "\n"
    << "rand seed: " << params.rand_seed << "\n"
    << "max doublings: " << params.max_doublings << "\n";
    // Convert from a char * to a std::string so we can append file names.
    std::vector<Edge> edges;
    bool read_succeeded = read_edges( params.input_file, edges );
    if (!read_succeeded) {
        std::cout << "Reading in edge list failed.\n";
        return 1;
    }
    // for (size_t edge_i = 0; edge_i < edges.size(); edge_i++) {
    //     std::cout << edges[edge_i].get_source() << "\t"
    //     << edges[edge_i].get_target() << "\n";
    // }
    // Initialize the walk agent.
    // std::cout << "initializing walk agent...\n";
    
    // Update for each segment length for which we have started a segment.
    // All segment lengths are powers of 2.
    // Visual:
    // 01
    //  |=segment start step, -=other step in a segment
    //  *=step at which we visited node A, the node visited on the last step
    //               *                                           *       *
    //  len=1  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| => y
    //  len=2  |-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| => y
    //  len=4  |---|---|---|---|---|---|---|---|---|---|---|---|---|---|-- => y
    //  len=8  |-------|-------|-------|-------|-------|-------|-------|-- => y
    //  len=16 |---------------|---------------|---------------|---------- => n
    //  len=32 |-------------------------------|-------------------------- => n
    //  len=64 |---------------------------------------------------------- => n
    // y: last step does contribute to segment mass
    // n: last step does not contribute to segment mass
    // steps_taken-1 is a multiple of 2^n => steps_taken-1 is a multiple of 2^(n-1).
    // That is, if this step is the start of a new segment for a given segment length,
    // then it must be the start of a new segment for all the shorter segment lengths.
    // Thus, we only need to check up to the first segment length
    // of which this is not the start of a new segment.
    // Example 1: steps_taken = 9, so (9-1) = 8 = 2^3
    // index=0, length=1: Trivially, every segment is the first step of a new 1-step segment.
    // index=1, length=2: Every odd-numbered segment is the first step of a 2-step segment.
    // index=2, length=4: We have completed two 4-step segments and are at step 1 of a third.
    // index=3, length=8: We have completed a single 8-step segment and are at step 1 of a second.
    // index=4, length=16: We have not completed any 16-step segments yet.
    // Example 1: steps_taken = 19, so (19-1) = 18 = 2^4 + 2
    // index=0, length=1: Trivially, every segment is the first step of a new 1-step segment.
    // index=1, length=2: Every odd-numbered segment is the first step of a 2-step segment.
    // index=2, length=4: We have completed 4 4-step segments and are at step 2 of a third.
    // index=3, length=8: We have completed 2 8-step segments and are at step 2 of a second.
    // index=4, length=16: We have completed 1 16-step segment and are at step 2 of a third.
    // index=5, length=32: We have not completed any 32-step segments yet.
    
    MBRWAgent agent = MBRWAgent(edges, params.bias, params.memory_length, params.start_vertex, params.rand_seed);
    // We do not know a priori how many segment lengths we will need to exhaust the whole network.
    // We keep doubling until we get a segment that has mass equal to num_nodes.
    // The constructor for the agent takes the first transition so that
    // last_edge and last_step_last_taken always have valid values.
    // This means we now have two nodes visited: the starting node and the one reached after 1 step.
    // First, update everything for the starting node:
    // The selection of the starting node counts as step 1 for our purposes here,
    // and moving to the next node counts as step 2.
    size_t steps_taken = 2;
    // step_last_visited[n] = step at which we last visited node n.
    // Initialize these to 0 so that
    // we count the first visit to each node as before the start of the first interval of any size
    // so that each first visit counts toward each segment mass.
    const size_t num_nodes = agent.get_num_nodes();
    std::vector<size_t> step_last_visited( num_nodes, 0 );
    step_last_visited[ agent.get_previous_node() ] = 1;
    step_last_visited[ agent.get_current_node() ] = 2;
    // We now have 3 segment lengths with indices 0 -> 2^0=1, 1 -> 2^1=2, and 2 -> 2^2=4.
    size_t last_length_index = 2;
    // segment_length[i] = 2^i for all i such that 2^i <= steps_taken.
    std::vector<size_t> segment_length = { 1, 2, 4 };
    // segment_mass_sum[i] = number of complete segments of length 2^i
    // We have completed two segments of length 1, 1 segment of length 2, and no segments of length 4.
    std::vector<size_t> num_segments = { 2, 1, 0 };
    // segment_start_step[i] = start step of the current segment in progress of length 2^i
    // The two length-1 segments and the 1 length-2 segment are complete.
    // The segments of each size will begin on step 3.
    // We are still building the first segment of length 4, which started on step 1.
    std::vector<size_t> segment_start_step = { 3, 3, 1 };
    // next_segment_start_step[i] = start step of the next segment of length 2^i
    std::vector<size_t> segment_end_step = { 3, 4, 4 };
    // partial_segment_mass[i] = mass so far of current segment in progress of length 2^i
    // segment_mass_sum[i] = sum of all the masses of completed segments of length 2^i
    // The two 1-step segments have masses of 1, which sum to 2.
    // We have not completed the 4-step segment, so it has a partial segment mass of 2.
    std::vector<size_t> partial_segment_mass = { 0, 0, 2 };
    // For segment lengths of 1 and 2, minimum and maximum are both 1 and 2 respectively.
    // For segment length 4, the theoretical minimum is 2 with backtracking, 3 without it,
    // while the theoretical maximum is 4.
    // Put the theoretical maximum in place of the minimum and theoretical minimum in place of the maximum
    // so that they will get replaced once we have actual values from complete segments.
    std::vector<size_t> min_mass = { 1, 2, 4 };
    std::vector<size_t> max_mass = { 1, 2, 2 };
    std::vector<double> segment_mass_power_sum = {};
    std::vector<double> log_mass_mean = {};
    // For q != 0, we are tracking the sum of segment masses raised to the power of q in segment_mass_power_sum
    // and then log of the sum over num segments divided by q (equivalent to log of the q-root of this sum over num segments) in log_mass_mean.
    // For q == 0, we instead track the sum of the logs of the segment masses masses in segment_mass_power_sum
    // and then this sum divided by the number of segments in log_mass_mean.
    // This is equivalent to the log of the geometric mean of the segment masses.
    // See https://en.wikipedia.org/wiki/Generalized_mean
    for (size_t order_index = 0; order_index < num_orders; order_index++) {
        double q = orders[order_index];
        if (q == 0) {
            double sum_1 = 2*std::log2( 1.0 );
            segment_mass_power_sum.push_back( sum_1 );
            log_mass_mean.push_back( sum_1/2.0 );
        } else {
            double sum_1 = 2*std::pow(1.0, q);
            segment_mass_power_sum.push_back( sum_1 );
            log_mass_mean.push_back( std::log2(sum_1/2.0)/q );
        }
    }
    // Because the graph cannot have self-loops, the first two nodes must be distinct.
    // Thus, the single 2-step segment has mass 2.
    for (size_t order_index = 0; order_index < num_orders; order_index++) {
        double q = orders[order_index];
        if (q == 0) {
            double sum_2 = std::log2( 2.0 );
            segment_mass_power_sum.push_back( sum_2 );
            log_mass_mean.push_back( sum_2 );
        } else {
            double sum_2 = std::pow( 2.0, q );
            segment_mass_power_sum.push_back( sum_2 );
            log_mass_mean.push_back( std::log2(sum_2)/q );
        }
    }
    // Fill in 0s for when we have a complete length-4 segment.
    for (size_t order_index = 0; order_index < num_orders; order_index++) {
        segment_mass_power_sum.push_back(0);
        log_mass_mean.push_back(0);
    }
    // Since we have only done trivial cases thus far, keep going.
    size_t new_segment_length = 4;
    size_t next_checkpoint = new_segment_length;
    bool has_not_saturated = true;
    bool has_not_converged = true;
    // Keep track of the total number of distinct nodes we have encountered.
    // We use this to initialize the partial segment mass for each new segment length,
    // since the first segment runs from the beginning of the walk up to some future step.
    size_t distinct_nodes_so_far = 2;
    std::cout << "walking...\n";
    // We have filled segments of lengths 1 (2^0) and 2 (2^1), so 1 doubling filled so far.
    for (size_t num_doublings = 1; num_doublings < params.max_doublings; num_doublings++) {
        
        // When we reach a given power of 2, start filling in a partial segment for the next power of 2.
        while ( steps_taken < next_checkpoint ) {
            
            agent.take_step();
            steps_taken++;
            size_t node = agent.get_current_node();
            
            distinct_nodes_so_far += (step_last_visited[node] == 0);
            size_t length_index;
            // Increment the length in progress until we get a segment in progress long enough that
            // it stretches back to when we last visited this node.
            for ( length_index = 0; (length_index <= last_length_index) && (step_last_visited[node] < segment_start_step[length_index]); length_index++ ) {
                partial_segment_mass[length_index]++;
            }
            // For each segment length for which this step completes a segment,
            size_t length_offset = 0;
            for (length_index = 0; (length_index <= last_length_index) && (steps_taken == segment_end_step[length_index]); length_index++) {
                // Keep track of the minimum and maximum masses for each segment length.
                size_t mass = partial_segment_mass[length_index];
                if (mass < min_mass[length_index]) {
                    min_mass[length_index] = mass;
                }
                if (mass > max_mass[length_index]) {
                    max_mass[length_index] = mass;
                }
                // For q != 0, we are tracking the sum of segment masses raised to the power of q in segment_mass_power_sum
                // and then log of the sum over num segments divided by q (equivalent to log of the q-root of this sum over num segments) in log_mass_mean.
                // For q == 0, we instead track the sum of the logs of the segment masses masses in segment_mass_power_sum
                // and then this sum divided by the number of segments in log_mass_mean.
                // This is equivalent to the log of the geometric mean of the segment masses.
                double mass_double = (double) mass;
                for (size_t order_index = 0; order_index < num_orders; order_index++) {
                    double q = orders[order_index];
                    if (q == 0.0) {
                        segment_mass_power_sum[ length_offset + order_index ] += std::log2( mass_double );
                    } else {
                        segment_mass_power_sum[ length_offset + order_index ] += std::pow( mass_double, q );
                    }
                }
                // increment the number of complete segments,
                num_segments[length_index]++;
                // shift the start and end points to the locations for the next segment,
                segment_start_step[length_index] += segment_length[length_index];
                segment_end_step[length_index] += segment_length[length_index];
                // and reset partial segment mass.
                partial_segment_mass[length_index] = 0;
                // Increment length_offset by num_orders, because we just filled in num_orders values.
                length_offset += num_orders;
            }
            
            step_last_visited[node] = steps_taken;
            
        }
        // std::cout << "completed " << steps_taken << " steps\n";
        
        // Check for the maximum change in log mean segment mass for any segment length.
        double change;
        double fractional_change;
        double max_change = 0.0;
        double max_fractional_change = 0.0;
        // Only check the change for cases where we have at least two complete segments,
        // since change is only meaningful for those cases.
        size_t length_offset = 0;
        double double_num_segments;
        for (size_t length_index = 0; length_index <= last_length_index; length_index++) {
            double_num_segments = (double) num_segments[length_index];
            for (size_t order_index = 0; order_index < num_orders; order_index++) {
                size_t mean_index = length_offset + order_index;
                double old_mean = log_mass_mean[mean_index];
                double q = orders[order_index];
                double new_mean;
                double new_sum = segment_mass_power_sum[mean_index];
                // For q != 0, we are tracking the sum of segment masses raised to the power of q in segment_mass_power_sum
                // and then log of the sum over num segments divided by q (equivalent to log of the q-root of this sum over num segments) in log_mass_mean.
                // For q == 0, we instead track the sum of the logs of the segment masses masses in segment_mass_power_sum
                // and then this sum divided by the number of segments in log_mass_mean.
                // This is equivalent to the log of the geometric mean of the segment masses.
                if (q == 0.0) {
                    new_mean = new_sum/double_num_segments;
                } else {
                    new_mean = std::log2( segment_mass_power_sum[mean_index]/double_num_segments )/q;
                }
                // We can only calculate the change if we have completed more than one segment of this length.
                if (num_segments[length_index] > 1) {
                    change = std::abs( new_mean - old_mean );
                    fractional_change = change/old_mean;
                    // std::cout << "length: " << segment_length[length_index] << ": " << num_segments[length_index] << " segments, old mean mass " << mean_segment_mass[length_index] << ", new mean mass " << new_mean_mass << ", change " << change << ", fractional change " << fractional_change << "\n";
                    if (max_change < change) {
                        max_change = change;
                    }
                    if (max_fractional_change < fractional_change) {
                        max_fractional_change = fractional_change;
                    }
                }
                // Now that we have calculated the change, we can replace the old value.
                log_mass_mean[mean_index] = new_mean;
            }
            length_offset += num_orders;
        }
        
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        
        std::cout << "step/max segment length = " << new_segment_length << " steps, min segment mass= " << min_mass[last_length_index] << " nodes, max mean segment mass change: " << max_change << ", max fractional change: " << max_fractional_change << ", seconds: " << elapsed_seconds.count() << "\n";
        // std::cout << "o=\t-Inf";
        // for (size_t order_index = 0; order_index < num_orders; order_index++) {
        //     std::cout << "," << orders[order_index];
        // }
        // std::cout << ",Inf\n";
        // for (size_t length_index = 0; length_index <= last_length_index; length_index++) {
        //     std::cout << "l=" << segment_length[length_index] << "\tm=" << std::log2( (double) min_mass[length_index] );
        //     length_offset = num_orders*length_index;
        //     for (size_t order_index = 0; order_index < num_orders; order_index++) {
        //         std::cout << "," << log_mass_mean[length_offset+order_index];
        //     }
        //     std::cout << "," << std::log2( (double) max_mass[length_index] ) << "\n";
        // }
        
        // Add a new segment length twice the previous longest segment length.
        new_segment_length = 2*new_segment_length;
        // std::cout << "adding segment length " << new_segment_length << "\n";
        segment_length.push_back(new_segment_length);
        // The first segment of a given length starts at step 1.
        segment_start_step.push_back(1);
        // The next segment of this length will start when we reach segment length + 1.
        segment_end_step.push_back(new_segment_length);
        // We start the new segment of the new length at the beginning of the walk,
        // so all distinct vertices encountered so far are included.
        partial_segment_mass.push_back(distinct_nodes_so_far);
        // Set the min segment mass to num_nodes so that it will get replaced when we have a complete segment.
        min_mass.push_back(num_nodes);
        // Set the max segment mass to 0 so that it will get replaced when we have a complete segment.
        max_mass.push_back(0);
        // We have not completed any segments of this length yet.
        num_segments.push_back(0);
        for (size_t order_index = 0; order_index < num_orders; order_index++) {
            segment_mass_power_sum.push_back(0);
            log_mass_mean.push_back(0);
        }
        // We now have a new last element in each array.
        last_length_index++;
        // If we have not yet saturated the network,
        // check again when we have completed a segment of the new length.
        next_checkpoint = new_segment_length;
        
        
    }
    // Make a new matrix that includes the logs of the min and max masses on the left and right of the other q-mean masses, respectively.
    size_t num_rows = last_length_index;
    size_t num_cols = 1+num_orders+1;
    std::vector<double> final_log_masses ( num_rows*num_cols );
    size_t length_offset = 0;
    size_t final_length_offset = 0;
    for (size_t length_index = 0; length_index < last_length_index; length_index++) {
        // Log the minimum.
        final_log_masses[final_length_offset] = std::log2( (double) min_mass[length_index] );
        // Copy over the log q-mean massses.
        for (size_t order_index = 0; order_index < num_orders; order_index++) {
            final_log_masses[final_length_offset + order_index + 1] = log_mass_mean[length_offset + order_index];
        }
        // Log the maximum.
        final_log_masses[final_length_offset + num_orders + 1] = std::log2( (double) max_mass[length_index] );
        length_offset += num_orders;
        final_length_offset += num_cols;
    }
    // std::cout << "saving results...\n";
    bool write_succeeded = write_matrix_to_csv(params.output_file, final_log_masses, num_rows, num_cols);
    if (!write_succeeded) {
        std::cout << "Writing segment mass generalized means failed.\n";
        return 1;
    }
    std::cout << "done\n";
    return 0;
}
