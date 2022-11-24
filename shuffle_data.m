function [subset_out,idx] = shuffle_data(subset_in)

% This function reorders in a random sequence the rows of the input matrix.
%
% Input: 
% subset        = observations
%
% Output: 
% subset_out    = shuffled observations
% idx           = indexes used in the permutation


% create a random permutation
[N,M] = size(subset_in);
idx   = randperm(N);
idx   = idx';

% initialize the output vector
subset_out = nan(N,M);

% shuffle
for j = 1:N
    subset_out(j,:) = subset_in(idx(j),:);
end


% This code has been written by Stefano Galelli.