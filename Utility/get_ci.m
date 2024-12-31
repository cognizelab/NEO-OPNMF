function CI = get_ci(H1, H2)
%% Concordance Index (CI)
% This function calculates the Concordance Index (CI) between two matrices H1 and H2.
% The CI measures the similarity of the pairwise cosine similarity matrices
% derived from H1 and H2, normalized to account for the number of items.
%
% Inputs:
%   H1 - A matrix of size (NumItems × NumDimensions), representing the first set of data.
%   H2 - A matrix of size (NumItems × NumDimensions), representing the second set of data.
%
% Output:
%   CI - The Concordance Index, a scalar value between 0 and 1.

% Step 1: Validate input dimensions
[NumItems1, NumDims1] = size(H1);
[NumItems2, NumDims2] = size(H2);
if NumItems1 ~= NumItems2 || NumDims1 ~= NumDims2
    error('H1 and H2 must have the same dimensions.');
end
NumItems = NumItems1; % Number of items (rows)

% Step 2: Compute the Euclidean norms of rows for H1 and H2
% Vectorized computation of row norms
rowNormsH1 = sqrt(sum(H1.^2, 2)); % Norms of each row in H1
rowNormsH2 = sqrt(sum(H2.^2, 2)); % Norms of each row in H2

% Step 3: Normalize rows of H1 and H2 by their respective norms
% This step converts rows into unit vectors (cosine normalization)
H1_normalized = H1 ./ rowNormsH1; % Normalize rows in H1
H2_normalized = H2 ./ rowNormsH2; % Normalize rows in H2

% Step 4: Compute pairwise cosine similarity matrices
% S1 and S2 are symmetric matrices where S(i,j) represents the cosine
% similarity between the i-th and j-th rows of the normalized matrices.
S1 = H1_normalized * H1_normalized'; % Similarity matrix for H1
S2 = H2_normalized * H2_normalized'; % Similarity matrix for H2

% Step 5: Compute the Frobenius norm of the difference between S1 and S2
similarityDifference = norm(S1 - S2, 'fro')^2;

% Step 6: Calculate the Concordance Index (CI)
% The CI formula normalizes the similarity difference by the number of
% pairwise comparisons (NumItems * NumItems - NumItems).
CI = 1 - (similarityDifference / (NumItems * NumItems - NumItems));
end