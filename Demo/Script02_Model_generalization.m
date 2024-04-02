%% --------- Model Generalization --------- 
% The script is used to quickly apply the pre-trained model to other independent samples.

%% Toolkit configuration

% Set up path for core codes.

clc; clear;  
base_dir = pwd;

addpath(genpath(fullfile(base_dir,'Code')));

%% Data loading

% Load a pre-trained model.

load(fullfile(base_dir,'Models','Model_NEO_FFI_60items_HCP_DYA.mat')); 

% Load data from an independent sample.

load(fullfile(base_dir,'Data','data_evaluation.mat')); 

%% Model generalization

model = Model.group_level_W{1}; % specify a pre-trained 2F model (basis matrix)
newdata = fmatrix(data(1:100,:),model); % calculate the dot product of data and pre-trained basis matrix

%% Important outputs for OPNMF

% newdata.data_dr: dimensionality reduced data based on the pre-trained model
% newdata.group_level_SED: group-level Reconfiguration Error that can be used to 
%                          assess the quality of model generalization (lower is better)
% newdata.individual_level_SED: individual-level Reconfiguration Error that can be used to 
%                               identify the individuals whose data are not well described by the model 