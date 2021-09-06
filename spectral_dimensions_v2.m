function Dq = spectral_dimensions_v2(log2_generalized_means, cutoff_mode, cutoff_value, log2_lengths)
%SPECTRAL_DIMENSIONS_V2 Compute generalized spectral dimensions
%   log2_lengths: a column vector of log2(walk segment length) values
%   differs from spectral_dimensions() in order of arguments
%   lets you use default values for log2_lengths
%   log2_generalized_means: a matrix where log2_generalized_means(i,:)
%   are log2 of the generalized means of the segment masses for segments
%   of length log2_lengths(i).
%   We do not need to know which generalized mean was used for each column
%   for the purposes of this calculation.
%   Dq: a row vector of the spectral dimensions
%       calculated as twice the slopes of lines fit to the q-means
%   We also choose some cutoff at which finite size effects mean
%   that the linear model no longer applies.
%   This can be set with either cutoff_index or cutoff_mass.
%   cutoff_length should be either a scalar integer or a vector of integers
%   of the same length as size(log2_generalized_means,2).
%   We then use log2_lengths(1:cutoff_index(j))
%   and log2_generalized_means(1:cutoff_index(j),j)
%   to calculate Dq(j).
%   If a value for cutoff_mass is passed in,
%   we use only (length,generalized mean mass) pairs for which
%   the mass is less than cutoff_mass.
%   If neither one is provided,
%   we use the maximum mass in the vector as the cutoff_mass.

[num_l, num_q] = size(log2_generalized_means);

% The masses output by segment_mass_log_multimeans are for these lengths.
if ~exist('log2_lengths','var')
    log2_lengths = ( 0:num_l-1 )';
end
% Make sure it is a column vector.
log2_lengths = log2_lengths(:);
if numel(log2_lengths) ~= num_l
    error('Size of length vector and number of rows in mass matrix must match.')
end

% By default, only look at segment lengths less than total number of nodes.
% If this length is not provided, assume we have saturated the network,
% so that the max segment mass achieved is the number of nodes.
if ~exist('cutoff_mode','var')
    cutoff_mode = 'length';
    cutoff_value = 2^max(log2_generalized_means,[],'all');
end
if strcmpi(cutoff_mode,'length')
    % If we receive a length, cut off lengths above it.
    ci = find( log2_lengths <= log2(cutoff_value), 1, 'last' );
    if isempty(ci)
        ci = num_l;
    end
    cutoff_index = repmat( ci, 1, num_q );
elseif strcmpi(cutoff_mode,'mass')
    % If we receive a mass, cut off any segment masses larger than it.
    cutoff_index = NaN(1,num_q);
    cutoff_log_mass = log2(cutoff_value);
    for j = 1:num_q
        ci = find( log2_generalized_means(:,j) <= cutoff_log_mass, 1, 'last' );
        if isempty(ci)
            cutoff_index(j) = num_l;
        else
            cutoff_index(j) = ci;
        end
    end
elseif strcmpi(cutoff_mode,'index')
    % For more fine-grained control,
    % allow the caller to pass in cutoff indices for indibidual q-values.
    cutoff_index = cutoff_value;
end
% Start with segments of length 4 or more.
% The no-backtracking rule constrains segment masses to be at least 3.
min_index = find( log2_lengths >= 2, 1 );
Dq = NaN(1,num_q);
for j = 1:num_q
    indices = min_index:cutoff_index(j);
    x = log2_lengths( indices );
    y = log2_generalized_means( indices, j );
    p = polyfit( x, y, 1 );
    Dq(j) = 2*p(1);
end

end

