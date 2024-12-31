function PS = get_ps(p1,p2,varargin)
%% Partition Similarity (PS)
% This function computes various similarity metrics between two partitions.
%
% Inputs:
%   p1 - Partition 1: A vector specifying the community assignment for each element.
%   p2 - Partition 2: A vector specifying the community assignment for each element.
%   varargin - Methods to compute similarity measures (optional):
%       'general'   - general similarity
%       'item'      - item-level similarity
%
% Output:
%   PS - A structure containing the calculated similarity metrics:
%       PS.zRand   - Z-score of Rand Similarity Coefficient
%       PS.SR      - Rand Similarity Coefficient
%       PS.SAR     - Adjusted Rand Similarity Coefficient
%       PS.VI      - Variation of Information
%       PS.NMI     - Normalized Mutual Information
%       PS.Jaccard - Jaccard Index
%       PS.Dice    - Dice Coefficient
%       PS.MV      - Modular Variability

idx = 1;
while idx <= (nargin - 2)
    switch varargin{idx}
        case {'general'}
            [PS.zRand,PS.SR,PS.SAR,PS.VI] = zrand(p1,p2);     
            PS.NMI = normalized_mutual_information(p1,p2);
            PS.Jacard = getJaccard(p1,p2);  
            PS.Dice = getDiceCoeff(p1,p2);    
        case {'item'}
            PS.MV = mv(p1,p2);   
    end
    idx = idx + 1;
end

end

function Jacard = getJaccard(A,B)
%% Jaccard Index    
A = agreement(A);
B = agreement(B);
Jacard = sum(A.*B)/(sum(A+B) - sum(A.*B));
end

function DC = getDiceCoeff(A,B)
%% Dice Coefficient    
A = agreement(A);
B = agreement(B);
DC = 2*(sum(A.*B))/sum(A + B);
end

function [zRand,SR,SAR,VI] = zrand(part1,part2)
%% ZRAND Calculates the z-Rand score and Variation of Information
%  distance between a pair of partitions.
%
%   [zRand,SR,SAR,VI] = ZRAND(part1,part2) calculates the z-score of the
%   Rand similarity coefficient between partitions part1 and part2. The
%   Rand similarity coefficient is an index of the similarity between the
%   partitions, corresponding to the fraction of node pairs identified the
%   same way by both partitions (either together in both or separate in
%   both)
%
%   NOTE: This code requires genlouvain.m to be on the MATLAB path
%
%   Inputs:     part1,  | Partitions that are being
%               part2,  | compared with one another
%
%   Outputs:    zRand,  z-score of the Rand similarity coefficient
%               SR,     Rand similarity coefficient
%               SAR,    Adjusted Rand similarity coefficient
%               VI,     Variation of information
%               
%
%   Amanda L. Traud, Eric D. Kelsic, Peter J. Mucha, and Mason A. Porter,
%   "Comparing Community Structure to Characteristics in Online Collegiate
%   Social Networks," SIAM Review 53, 526-543 (2011).
if size(part1, 1) == 1
    part1 = part1';
end
if size(part2, 1)==1
    part2 = part2';
end
if length(part1) ~= length(part2)
    disp('ERROR: partitions not of equal length')
    return
end
nij = sparse(part1+1, part2+1, 1);
ni = sum(nij, 2);
nj = sum(nij, 1);
n = length(part1);
M = n*(n-1)/2;
M1 = sum(ni.^2-ni)/2;
M2 = sum(nj.^2-nj)/2;
a = full(sum(sum(nij.^2-nij)))/2; %same in both
b = M1-a;                         %same in 1, diff in 2
c = M2-a;                         %same in 2, diff in 1
d = M-(a+b+c);                    %diff in both
SR = (a+d)/(a+b+c+d);
meana = M1*M2/M;
SAR = (a-meana)/((M1+M2)/2-meana);
C1 = 4*sum(ni.^3)-8*(n+1)*M1+n*(n^2-3*n-2);
C2 = 4*sum(nj.^3)-8*(n+1)*M2+n*(n^2-3*n-2);
vara = M/16 - (4*M1-2*M)^2*(4*M2-2*M)^2/(256*M^2) + C1*C2/(16*n*(n-1)*(n-2)) + ...
    ((4*M1-2*M)^2-4*C1-4*M)*((4*M2-2*M)^2-4*C2-4*M)/(64*n*(n-1)*(n-2)*(n-3));
zRand = (a-meana)/sqrt(vara);

c1 = unique(part1);
c2 = unique(part2);
H1 = 0; H2 = 0; I = 0;
for i = c1'
    pi = ni(i+1)/n;
    H1 = H1-pi*log(pi);
    for j = c2'
        if nij(i+1,j+1)
            pj = nj(j+1)/n;
            pij = nij(i+1,j+1)/n;
            I = I+pij*log(pij/pi/pj);
        end
    end
end
for j = c2'
    pj = nj(j+1)/n;
    H2 = H2-pj*log(pj);
end
VI = (H1+H2-2*I);

end

function [nmi] = normalized_mutual_information(partition1,partition2,varargin)
%% NORMALIZED_MUTUAL_INFORMATION - This function computes the normalized mutual
%         information (normalization taken from from Danon 2005, note other normalizations
%         exist) which is a measure of similarity between two different partitions of 
%         community structure.  The measure is bounded in [0 1] where 1 implies perfect
%         agreement between partitions.  In cases where the number of nodes
%         divided by the number of clusters <~ 100, the adjusted value should be
%         used to correct for chance (see Vinh et al 2010).
% 
% Inputs: 
%         partition1 - vector containing community assignment for each node in first partition
% 
%         partition2 - vector containing community assignment for each node in second partition
% 
%         varargin - set equal to 'adjusted' if want to correct for chance and compute the AMI 
%               (use this option if the number of nodes diveded by the number of communities <~ 100.  
%               Default is the unadjusted NMI (set varaible equal to 'unadjusted' or leave empty).
% 
% Outputs: 
%         nmi - normalized mutual information.  Gives measure of
%               similarity between partitons.  can be adjusted for chance 
%               (AMI- adjusted normalized mutual information - see Vihn 2010)
% 
% Written by Sarah Feldt Muldoon

% References: Danon L, Diaz-Guilera A, Duch J, Arenas A (2005) Comparing
%         community structure identification. J Stat Mech:P09008.; Vinh NX, Epps J, 
%         Bailey J (2010) Information Theoretic Measures for Clusterings Comparison: 
%         Variants, Properties, Normalization and Correction for Chance. The Journal
%         of Machine Learning Research.

%if no keyword, calculate the unadjusted NMI
if isempty(varargin)
    varargin{1} = 'unadjusted';
end

num_nodes = length(partition1);
comm1_array = unique(partition1);
comm2_array = unique(partition2);
num_comm1 = length(comm1_array);
num_comm2 = length(comm2_array);
contingency_matrix = zeros(num_comm1,num_comm2);

%make contingency table (also called confusion_matrix)
for i = 1:num_comm1
    comm1 = comm1_array(i);
    for j = 1:num_comm2
        comm2 = comm2_array(j);
        contingency_matrix(i, j) = sum(partition1 == comm1 & partition2 == comm2);
    end
end

sum_over_i = sum(contingency_matrix,2);
sum_over_j = sum(contingency_matrix,1);
total_sum = sum(sum(contingency_matrix));

%compute mutual information
mi_matrix = contingency_matrix .* log((contingency_matrix .* total_sum) ./ (sum_over_i * sum_over_j));
mi = sum(mi_matrix(isfinite(mi_matrix))) / total_sum;

%compute normalization terms (we normalize by (h1+h2)/2 as in Danon 2005)
h1 = sum_over_i .* log(sum_over_i ./ total_sum);
h1 = -sum(h1(isfinite(h1))) / total_sum;
h2 = sum_over_j .* log(sum_over_j ./ total_sum);
h2 = -sum(h2(isfinite(h2))) / total_sum;

if strcmp(varargin{1}, 'adjusted') == 1
    %compute expected mutual information (Vihn 2010)
    expected_mi = 0;
    for i = 1:num_comm1
        for j=1:num_comm2
            k_min = max([sum_over_i(i) + sum_over_j(j) - total_sum, 1]);
            k_max = min([sum_over_i(i), sum_over_j(j)]);
            for k = k_min:k_max
                term1 = k / total_sum;
                term2 = log(total_sum * k / (sum_over_i(i) * sum_over_j(j)));
                term3_numerator = factorial(sum_over_i(i)) * factorial(sum_over_j(j)) * factorial(total_sum - sum_over_i(i)) * factorial(total_sum - sum_over_j(j));
                term3_denominator = factorial(total_sum) * factorial(k) * factorial(sum_over_i(i) - k) * factorial(sum_over_j(j) - k) * factorial(total_sum - sum_over_i(i) - sum_over_j(j) + k);
                emi_term = term1 * term2 * term3_numerator / term3_denominator;
                if isfinite(emi_term)
                    expected_mi = expected_mi + emi_term;
                end
            end
        end
    end
    %compute AMI_sum as in Vihn 2010
    nmi = (mi - expected_mi) / (.5 * (h1 + h2) - expected_mi);
elseif strcmp(varargin{1}, 'unadjusted') == 1
    %compute NMI as in Danon 2005 (NMI_sum in Vihn 2010)
    nmi = 2 * mi / (h1 + h2);
end

end

function [MV,cMV] = mv(N1,N2,ID)
%% Calculation of Modular Variability (MV) for A Given Node
%
%  Z.K.X. 2017/12/15
%-----------------------------------------------------------------------------------%
%   MV: For each node, its module affiliation variability is calculated 
%       between two modular partitions through a metric of modular variability (MV) 
%
%                      |Xk(i)∩Xk(j)|   |Xk(i)∩Xk(j)|
%       MVk(i,j) = 1 - -------------- * --------------
%                          |Xk(i)|          |Xk(j)|
%       where Xk (i) and Xk (j) denote the module labels to which node k belongs in 
%       modular partitions i and j, respectively. Xk (i)∩Xk (j) represents the 
%       common node set between modules Xk (i) and Xk (j), and |Xk(i)∩Xk(j)| 
%       denotes the number of nodes in the common node set. For node k between 
%       modular structures i and j, a small overlap between the two modules Xk (i) 
%       and Xk (j) indicates large module affiliation variability.
%  cMV: For each node, the MV values are calculated for every networks identified by 
%       ID. The cMV indicates how possible the community pattern from a given node to 
%       a given network would be changed. 
%-----------------------------------------------------------------------------------%
if (nargin < 3)
	ID = [];
end
if size(N1,2) == 1 
    N1 = agreement(N1);
    N2 = agreement(N2);
end
for i = 1:length(N1)
    M1 = sum(N1(i,:))+1;
    M2 = sum(N2(i,:))+1;
    cover = find(N1(i,:)==1&N2(i,:)==1);
    com = length(cover)+1;
    MV(i,1) = 1-(com/M1)*(com/M2);
end
if ~isempty(ID)
    for i = 1:length(N1)
        for j = 1:max(ID)
            F = find(ID==j);
            M1 = sum(N1(i,F))+1;
            M2 = sum(N2(i,F))+1;
            cover = find(N1(i,F)==1&N2(i,F)==1);
            com = length(cover)+1;
            cMV(i,j) = 1-(com/M1)*(com/M2);
        end
    end
end

end

function D = agreement(ci,buffsz)
%% AGREEMENT Agreement matrix from clusters
%
%   D = AGREEMENT(CI) takes as input a set of vertex partitions CI of
%   dimensions [vertex x partition]. Each column in CI contains the
%   assignments of each vertex to a class/community/module. This function
%   aggregates the partitions in CI into a square [vertex x vertex]
%   agreement matrix D, whose elements indicate the number of times any two
%   vertices were assigned to the same class.
%
%   In the case that the number of nodes and partitions in CI is large
%   (greater than ~1000 nodes or greater than ~1000 partitions), the script
%   can be made faster by computing D in pieces. The optional input BUFFSZ
%   determines the size of each piece. Trial and error has found that
%   BUFFSZ ~ 150 works well.
%
%   Inputs,     CI,     set of (possibly) degenerate partitions
%               BUFFSZ, optional second argument to set buffer size
%
%   Outputs:    D,      agreement matrix
%
%   Richard Betzel, Indiana University, 2012

%modification history
%09.24.2012 - added loop for big N that makes the function slower but also
% prevents it from maxing out memory.

n = size(ci,2);

if nargin < 2
    buffsz = 1000;
end

if n <= buffsz
    
    ind = dummyvar(ci);
    D = ind*ind';
    
else
    
    a = 1:buffsz:n;
    b = buffsz:buffsz:n;
    
    if length(a) ~= length(b)
        b = [b, n];
    end
    
    x = [a' b'];
    nbuff = size(x,1);
    
    D = zeros(size(ci,1));
    for i = 1:nbuff
       y = ci(:,x(i,1):x(i,2));
       ind = dummyvar(y);
       D = D + ind*ind';
    end
    
end

D = D.*~eye(length(D));
end