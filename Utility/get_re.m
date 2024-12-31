function RE = get_re(W, V)
%% Reconstruction Error (RE)
% Computes the mean absolute reconstruction error between V and its projection onto W.
%
% Inputs:
%   W - A matrix of size (items × dimension), representing the latent space or basis.
%   V - A matrix of size (items × subject), representing the target data.
%
% Output:
%   RE - The mean absolute reconstruction error.

% Validate input dimensions
if size(W, 1) ~= size(V, 1)
    error('The number of rows in W and V must match.');
end

% Compute the projection of V onto W
projection = W * (W' * V);

% Compute the absolute reconstruction error
error_matrix = abs(round(projection) - V);

% Sum the error across rows (for each subject) and then take the mean
RE = mean(sum(error_matrix, 1));

end