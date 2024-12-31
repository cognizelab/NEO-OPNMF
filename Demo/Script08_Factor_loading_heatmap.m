%% --------- Factor Loading Heatmap --------- 
% This script is used to demonstrate the factor structure based on 
% classical dimensionality reduction models (e.g., factor analysis or principal component analysis).
% We will skip the model evaluation phase and use a more flexible approach to display the loading heatmap.

%% Toolkit configuration

% Set up path for core codes.

clc; clear;  
base_dir = pwd; % root directory

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

%% Questionnaire basic information configuration
item = T.Description;
ID = T.ID;
dl = {'Neuroticism','Extrovesion','Openness','Agreeableness','Conscientiousness'};

%% Factor analysis

% Set parameters for factor analysis.

model.fmatrix = 'PAF'; % select a dimensionality reduction method
                       % 'OPNMF' --> Orthogonal Projective Non-negative Matrix Factorization  
                       % 'MLE' --> Maximum Likelihood Estimate (Factor Analysis)
                       % 'PAF' --> Principal Axis Factoring (Factor Analysis)
                       % 'PCA' --> Principal Component Analysis

model.factor = [2 5]; % only train two-factor and five-factor models (based on prior assumptions)
model.group_level = 2; % 2 --> only generate group-level models and skip the cross-validation steps 

DM.missing = 'mode'; % 'mode' --> replace missing values with mode of a given item
      
CV = []; % no parameters need to be configured

OS.filename = 'test'; % set a name for the output file (default --> 'FMATRIX')   
OS.fa_rotate = 'promax'; % perform oblique rotation based on the promax factor rotation method in MATLAB
 
[Model,EVA,Log,Data] = fmatrix(data,model,DM,CV,OS); % execute the factor analysis

%% Five-factor model

% Display the loading matrix of the five-factor model.

[newmatrix,showitem,ml] = plot_factor_loading(Model.group_level_W{2},'ID',ID,'itemlabel',item,'dimlabel',dl,'itemcolor',Color_NEOAC);

% Adjust the presentation scheme based on the generated heatmap.

colorlimits = [-0.6 0.6]; % restrict the color bar to a specific numerical range
revdim = [1 2 5]; % reverse the loading direction of certain dimensions

% Generate the adjusted loading matrix.

[newmatrix,showitem,ml] = plot_factor_loading(Model.group_level_W{2},'ID',ID,'itemlabel',item,'dimlabel',dl,'itemcolor',Color_NEOAC,'colorlimits',colorlimits,'revdim',revdim);

%% Two-factor model

% Display the loading matrix of the two-factor model.

[newmatrix,showitem,ml] = plot_factor_loading(Model.group_level_W{1},'ID',ID,'itemlabel',item,'dimlabel',dl,'itemcolor',Color_NEOAC);

% Adjust the presentation scheme based on the generated heatmap.

colorlimits = [-0.6 0.6]; % restrict the color bar to a specific numerical range
revdim = [ ];  

[newmatrix,item_all,ml] = plot_factor_loading(Model.group_level_W{1},'ID',ID,'itemlabel',item,'dimlabel',dl,'itemcolor',Color_NEOAC,'colorlimits',colorlimits,'revdim',revdim);

%% Two-factor model with reordered loadings

% Display the loading coefficients for each dimension, ordered by the magnitude of the loading coefficients.

[newmatrix,item_d1,ml] = plot_factor_loading(Model.group_level_W{1}(:,1),'ID',ID,'itemlabel',item,'dimlabel',dl,'itemcolor',Color_NEOAC,'direction','positive','itemsize',10);
[newmatrix,item_d2,ml] = plot_factor_loading(Model.group_level_W{1}(:,2),'ID',ID,'itemlabel',item,'dimlabel',dl,'itemcolor',Color_NEOAC,'direction','positive','itemsize',10);

% Refer to the documentation by entering doc plot_factor_loading to access more detailed information about the input options and their functionalities.