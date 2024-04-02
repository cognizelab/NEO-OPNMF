%% --------- Basic Evaluation --------- 
% This script is used to conduct the basic version of the model evaluation 
% as well as the visualization of primary outcomes.
% The evaluation pipleine has been tested in MatLab 2018a and later.

%% Toolkit configuration

% Set up path for core codes.

clc; clear;  
base_dir = pwd;

addpath(genpath(fullfile(base_dir,'Code')));

%% Data loading

% Load demo data.
% "data" is a subjects (rows) * items (colums) matrix.
% "sex" and "age" are corresponding demographic variables.
% "datasets" defines which dataset each data point belongs to.

load(fullfile(base_dir,'Data','data_evaluation.mat')); 

%% Score setting

% Recode reverse scoring items (range: 0-4).
% "reverse" marks items that are scored in reverse.
% Items that require reverse scoring are marked as 1 and other items are marked as 0.

f = find(reverse==1); % reverse-scored items
data(:,f) = max(max(data)) - data(:,f); % turn all items into the same direction

%% Model evaluation

% Execute the evaluation procedure with the default configurations.
[Model,EVA,Log,Data] = fmatrix(data); % execute the evaluation procedure with default parameters

% Set parameters for basic version of model evaluation.

model.fmatrix = 'OPNMF'; % choose a model ('OPNMF' or 'FA'; default --> 'OPNMF')
model.factor = [2:8]; % specify the range of factor decomposition (default --> [1:9])
model.group_level = 1; % if generate a group-level model
                       % 0 --> do not generate the group-level model
                       % 1 --> generate group-level models and perform cross-validation evaluations (default)
                       % 2 --> only generate group-level models and skip the other steps
model.outlier = []; % if conduct outlier filtering (2 vectors double)
                    % [3 0.5] --> subjects are defined as outliers and subsequently removed if their fitting errors 
                    %             surpassed 3 scaled MAD in more than half of the matrix decompositions
                    % [ ] --> do not execute this operation (default)

DM.missing = 'mode'; % addressing missing values ('remove' or 'mode'; default --> 'remove')
                     % 'remove' --> remove data with any missing value
                     % 'mode' --> replace missing values with mode of a given item
      
CV.n = 100; % choose the CV number (default --> CV.n = 5000)
CV.method = 'K-fold'; % choose a resampling method ('K-fold' or 'Bootstrap'; default --> 'K-fold')
CV.K = 5; % choose a K value for K-fold cross-valudation (default --> CV.K = 5)
CV.direction = 'one-way'; % specify the direction of model generalization ('one-way' or 'two-way'; default --> 'one-way')
                          % 'one-way' --> generalize the models trained on a large portion of data (e.g., 4/5 in 5-folds CV)
                          %                to a small portion of data (e.g., 1/5 in 5-folds CV)
                          % 'two-way' --> when using a split-half resampling strategy (CV.K = 2), 
                          %               since the data on both sides are equal (i.e., 1/2 vs. 1/2), 
                          %               the same operations can be performed separatel

OS.filename = 'test';     % set a name for the output file (default --> 'FMATRIX')
                          % all output data will be saved to this file

[Model,EVA,Log,Data] = fmatrix(data,model,DM,CV,OS); % execute the evaluation procedure 

%% Important inputs for OPNMF

% data: subjects * items matrix  
% model: model configuration   
% DM: data manipulation   
% CV: cross validation setting

% Detailed information about the input parameters can be found in the fmatrix.m file.

%% Important outputs for OPNMF

% Model.group_level_W: coefficients of the basis matrix that are trained using the entire set of data and can be easily generalized to new data
% Model.group_level_H: individual scores in the loading matrix that are available for follow-up analysis (e.g., clustering and individual difference)

% EVA.aRand: Adjusted Rand Similarity Coefficient
% EVA.CI: Concordance Index
% EVA.VI: Variation of Information
% EVA.increased_SED: Increased Reconfiguration Error
% EVA.MV: Item Variability

% Log: log of the model evaluation parameters 

% Data: pre-processed data that is fed into the evaluation pipeline

%% Visualization

% Display the results of model evaluation and the basis matrices for selected models, based on the outputs of the previous stage.
% The following visualization requires MatLab 2021a or later to call the function of "boxchart".
% "T" recodes the questionnaire information of NEO-FFI.
% "Color_NEOAC" indicates n (number of theoretical dimensions) * 3 RGB colors for theoretical dimensions.

parameter.coefficient = [2 5]; % select the models you wish to display (e.g., [2 5] means "displaying 2F and 5F models"; "coefficient" handles the heatmap)
parameter.items = T.Description; % set the label for each item
parameter.ID = T.ID; % set theoretical dimensions for each item
parameter.dimlabel = {'Neuroticism','Extrovesion','Openness','Agreeableness','Conscientiousness'}; % set labels for theoretical dimensions
parameter.itemcolor = Color_NEOAC; % set RGB colors for theoretical dimensions

plot_fmatrix(Model,EVA,Log,Data,parameter); % display figures

% The HCP-based model can be demonstrated by loading 'Model_NEO_FFI_60items_HCP_DYA.mat'. 

clearvars -except parameter base_dir
load(fullfile(base_dir,'Models','Model_NEO_FFI_60items_HCP_DYA.mat'));

plot_fmatrix(Model,EVA,Log,Data,parameter); % display figures based on HCP data  