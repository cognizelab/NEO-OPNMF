%% --------- Consolidation of Evaluation Results --------- 

% This script is used to combine results from two OPNMF outputs. 
% This will facilitate the parallel execution of the evaluation process, especially for large "CV.n".
% The evaluation pipleine has been tested in MatLab 2018a and later.

%% Toolkit configuration

% Set up path for core codes.

clc; clear;  
base_dir = pwd;

addpath(genpath(fullfile(base_dir,'Code')));

%% Data loading

% Load demo data.

load(fullfile(base_dir,'Data','data_evaluation.mat')); 

%% Score setting

% Recode reverse scoring items (range: 0-4).
% "reverse" marks items that are scored in reverse.
% Items that require reverse scoring are marked as 1 and other items are marked as 0.

f = find(reverse==1); % reverse-scored items
data(:,f) = max(max(data)) - data(:,f); % turn all items into the same direction

%% Model evaluation

CV.n = 10;

% Run the first evaluation.

OS.filename = 'test1';      
[Model,EVA,Log,Data] = fmatrix(data,[],[],CV,OS); 

% Run the second evaluation.

OS.filename = 'test2';      
[Model,EVA,Log,Data] = fmatrix(data,[],[],CV,OS); 

%% Result combination

% Merge the results of the two evaluations.

output = 'test'; % set a name for the output file (default --> 'FMATRIX')

[Model,EVA,Log,Data] = combine_fmatrix('test1','test2',output); % merge the results from 'test1' and 'test2'