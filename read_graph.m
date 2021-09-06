function G = read_graph(edge_file_name,node_name_file_name)
%READ_GRAPH Read in the edges and, optionally the node names, of a graph.
%   edge_file_name: name of a text file with
%   each line representing an edge
%   as a pair of node indices separated by tabs
%   node_name_file_name: name of a text file with
%   each line representing the name of an edge

edges = readmatrix(edge_file_name,'FileType','text','Delimiter','\t');
s = edges(:,1);
t = edges(:,2);
if exist('node_name_file_name','var')
    node_names = readcell(node_name_file_name,'FileType','text');
    % MATLAB's graph() function expects node weights before node names,
    % but we are only working with unweighted graphs.
    weights = ones( size(s) );
    G = graph(s,t,weights,node_names);
else
    G = graph(s,t);
end

