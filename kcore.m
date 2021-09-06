function G_kcore = kcore(G,k)
%KCORE Find the largest subgraph such that all degrees >= k.
%   G: a MATLAB graph object
%   k: a scalar, integer representing the minimum node degree to allow
%   G_kcore: the largest subgraph of G such that
%   all nodes have degree at least k.
%   It may be empty if no such subgraph exists.

G_kcore = G;
findlessnodes = @(g) find( degree(g) < k );
lessnodes = findlessnodes(G_kcore);
i = 0;
while ~isempty(lessnodes)
    i = i+1;
    fprintf('iteration %u: %u nodes, %u edges remaining, removing %u...\n', ...
        i, numnodes(G_kcore), numedges(G_kcore), numel(lessnodes) )
    G_kcore = rmnode(G_kcore,lessnodes);
    lessnodes = findlessnodes(G_kcore);
end
disp('done')

end

