function [edges, gene_names] = extract_full_gunsalus_et_al_2005(raw_file_path,edges_file_path,gene_names_file_path)
%EXTRACT_FULL_GUNSALUS_ET_AL_2005 Put interaction edges in a file.
%   Author: Adam Craig, Last Updated: 2021-09-06
%   The network file downloaded from Gunsalus et al. 2005
%   represents the network with three columns:
%   Gene A, abreviation for an interaction type, Gene B
%   pp: direct protein-protein interactions,
%       the type of interactions we select.
%   The output network uses indices (numbers) instead of gene names,
%   but we save the gene names in the correct order to another file.

edges_with_type = readcell(raw_file_path,'FileType','text','Delimiter','\t');
is_right_type = contains( edges_with_type(:,2), 'pp' );
edges_by_name = edges_with_type( is_right_type, [1 3] );
[gene_names, ~, stubs] = unique(edges_by_name);
edges = reshape(stubs, [], 2);
if exist('edges_file_path','var')
    writematrix(edges,edges_file_path,'FileType','text','Delimiter','\t')
end
if exist('gene_names_file_path','var')
    writecell(gene_names,gene_names_file_path,'FileType','text','Delimiter','\t')
end

end

