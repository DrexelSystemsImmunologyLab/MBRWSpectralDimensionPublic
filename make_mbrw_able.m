function [new_edge_file_names, new_node_name_file_names] = make_mbrw_able(edge_file_name, node_name_file_name)
%MAKEMBRWABLE Create versions of a network usable with MBRW.
%   edge_file_name: the name of the text file with the original network
%   represented as edges with endpoints indices separated by tab characters
%   Endpoint indices should be positive integers.
%   node_name_file_name: the name of the text file with node names
%   new_edge_file_names: a cell array of character vectors representing
%   names of files in the same folder with _cc_[connected component number]
%   inserted before the file extension.
%   Components are in descending order of node count.
%   Each file stores a connected component of the 2-core
%   (largest subgraph such that all nodes have degree at least 2)
%   of the simple undirected graph version of the original graph.
%   new_edge_file_names: a cell array of character vectors representing
%   printout:
%   We also print out information about
%   the time each section of the function takes
%   and the size of the network at that stage.
%   If the 2-core of the graph is empty, we do not create any files,
%   and we return an empty cell array for newfilenames.

disp('reading in graph...')
tic
G = read_graph(edge_file_name,node_name_file_name);
toc
fprintf( 'has %u nodes and %u edges\n', ...
    numnodes(G), numedges(G) )

disp('simplifying...')
tic
G_simple = simplify(G);
toc
fprintf( '%u nodes, %u edges remain\n', ...
    numnodes(G_simple), numedges(G_simple) )

disp('taking 2-core...')
tic
G_2core = kcore(G_simple,2);
toc
fprintf( '%u nodes, %u edges remain\n', ...
    numnodes(G_2core), numedges(G_2core) )

% If the 2-core is empty, stop here.
% Do not make any more files.
if numnodes(G_2core) == 0
    new_edge_file_names = {};
    return
end

disp('finding connected components...')
tic
[bins, binsizes] = conncomp(G_2core);
ubins = unique(bins);
numbins = numel(ubins);
% Put the components in descending order of size.
[~, sort_order] = sort(binsizes,'descend');
ubins_sorted = ubins(sort_order);
toc
fprintf('found %u connected components\n', numbins)

% Separate the file extension from the rest of the file name.
% We will insert connected component numbers between them.
[base_edge_file_name, edge_file_extension] = ...
    rmfileextension(edge_file_name);
[base_node_name_file_name, node_name_file_extension] = ...
    rmfileextension(node_name_file_name);
new_edge_file_names = cell(numbins,1);
new_node_name_file_names = cell(numbins,1);
tic
for cc = 1:numbins
    
    fprintf('retrieving connected component %u...\n', cc)
    G_cc = subgraph( G_2core,  bins == ubins_sorted(cc) );
    fprintf('%u nodes, %u edges\n', ...
        numnodes(G_cc), numedges(G_cc) )
    
    [s, t] = findedge(G_cc);
    node_names = G_cc.Nodes.Name;
    cc_edge_file_name = [base_edge_file_name ...
        sprintf('_cc_%u',cc) edge_file_extension];
    cc_node_name_file_name = [base_node_name_file_name ...
        sprintf('_cc_%u',cc) node_name_file_extension];
    fprintf('writing to file %s...\n', cc_edge_file_name)
    writematrix([s,t],cc_edge_file_name,'FileType','text','Delimiter','\t')
    fprintf('writing to file %s...\n', cc_node_name_file_name)
    writecell(node_names,cc_node_name_file_name,'FileType','text')
    new_edge_file_names{cc} = cc_edge_file_name;
    new_node_name_file_names{cc} = cc_node_name_file_name;
    
end
disp('done')
toc

end

