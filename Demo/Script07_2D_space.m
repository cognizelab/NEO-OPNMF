%% --------- Projecting 2 Clusters onto A 2D space --------- 
% This script is used to display the volatility of the same items under different allocation schemes.

%% Toolkit configuration

% Set up path for core codes.

clc; clear;  
base_dir = pwd;

addpath(genpath(fullfile(base_dir,'Code')));

%% Data loading

% "data" are the residurals of Meta-Two traits after controling for sex and age.
% "loadings" are the Pearson's correlations between Big-Five and Meta-Two traits.
% "T" involves the row data of Meta-Two traits in HCP-YAD.

load(fullfile(base_dir,'Data','data_HCP_meta2traits.mat')); 

%% Fuzzy C-Means Clustering  
N = 2; % set the number of clusters

exp = 2; % set the fuzziness exponent, which determines the level of cluster fuzziness
maxIter = 100; % set the maximum number of iterations for the FCM algorithm
minImprove = 0.00001; % set the minimum amount of objective function improvement between two consecutive iterations
displayObjective = false; % set whether to display the objective function value at each iteration 
options = [exp maxIter minImprove displayObjective]; % combine all the options into a single array

scores = zscore(data); % standardize the data by subtracting the mean and dividing by the standard deviation

[C,U] = fcm(scores,N,options); % run the Fuzzy C-Means (FCM) clustering algorithm on the standardized data

% C is the cluster centers
% U is the fuzzy partition matrix
 
[~, group] = max(U); % identify groups base on 'U'

%% Biplot

S.group = group'; % define different groups
S.group_center = C; % S.group_center = [];
S.density_plot = 1; % define the transparency by the data density within the range where the data point is located  
S.line_label = {'Conscientiousness','Agreeableness','Extrovesion','Openness','Neuroticism'}; % labels for the Big-Five traits

S.dot_degree_color{1} = [0.552941176470588,0.627450980392157,0.796078431372549]; % RGB color for group 1 (filled color)
S.dot_degree_color{2} = [0.988235294117647,0.552941176470588,0.384313725490196]; % RGB color for group 2 (filled color)
S.dot_color = [S.dot_degree_color{1};S.dot_degree_color{2}]; % RGB color for the two groups (outline color)

Y = plot_biplot(loadings,scores,S); % display figure

% Please try entering 'doc plot_biplot' for more details about the input options.