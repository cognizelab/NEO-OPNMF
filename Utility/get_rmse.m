function RMSE = get_rmse(W, V)
% RMSE: Root Mean Square Error between V and the projection of V onto W.
%
% Inputs:
%   W - A matrix of size (items × dimension), representing the basis or subspace.
%   V - A matrix of size (items × subject), representing the target data.
%
% Output:
%   RMSE - The root mean square error between V and its projection onto W.

% Check that the number of rows in W and V match
if size(W, 1) ~= size(V, 1)
    error('The number of rows in W and V must match.');
end

% Calculate the projection of V onto W
% V' * W results in a (subject × dimension) matrix
% Taking the transpose of (V' * W) results in a (dimension × subject) matrix
% W * ((V' * W)') gives a (items × subject) projection
projection = W * ((V' * W)');

% Compute the error matrix between V and its projection
error_matrix = V - projection;

% Compute the Frobenius norm of the error matrix
frobenius_norm = norm(error_matrix, 'fro');

% Normalize the Frobenius norm by the total number of elements in V
% sqrt(numel(V)) ensures RMSE is normalized by the total size of V
RMSE = frobenius_norm / sqrt(numel(V));

end