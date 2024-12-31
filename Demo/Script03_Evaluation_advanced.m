%% --------- Advanced Evaluation --------- 
% This script is used to conduct the advanced versions of the model evaluation 
% as well as the visualization of primary outcomes.
% The evaluation pipleine has been tested in MatLab 2018a and later.

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

%% Model evaluation based on Bootstrap

% In the following example, we will apply the OPNMF model trained in dataset 1 to datasets 2 and 3. 
% Here, datasets 2 and 3 are subjected to bootstrap resampling in each run.

% Set parameters for model evaluation.

model.factor = [2:5];  

DM.group = datasets; % a column vector of classified (e.g., family relationship, study batch, datasets, gender, culture\country, etc) 
                    % or continuous (e.g., age) data labels
DM.resample_subset = [2 3]; % group ID to be resampled (e.g., [2 3] --> resample data points numbered 2 and 3)
DM.resample_method = 'Bootstrap'; % resampling method within each dataset
                                  % 'random' --> randomly select N data points
                                  % 'K-fold' --> randomly select 1/K data points  
                                  % 'Bootstrap' --> generate an equal number of datasets, allowing for duplication of elements  
      
CV.n = 100;  
CV.direction = 'mixed'; % perform cross-sample or cross-dataset cross validation
CV.st = [1 2; 1 3]; % source and target datasets of the generalization direction if the 'mixed' strategy is selected
                    % (e.g., [1 2; 1 3] --> train a model in dataset 1 then generalize it to dataset 2 and 3) 

OS.filename = 'test';      
 
[Model,EVA,Log,Data] = fmatrix(data,model,DM,CV,OS);  

% Display the results.
parameter.coefficient = [2 5];  
parameter.items = T.Description;  
parameter.ID = T.ID;  
parameter.dimlabel = {'Neuroticism','Extrovesion','Openness','Agreeableness','Conscientiousness'};  
parameter.itemcolor = Color_NEOAC; 

plot_fmatrix(Model,EVA,Log,Data,parameter);  

%% Model evaluation based on cross-validation with group information considered

% In the following example, we will perform five-fold cross-validation, 
% controlling the grouping information in the dataset, so as to never split data with the same ID between folds.

% Generate group labels.

datasets = repmat([1:20],15,1); datasets = datasets(:); % generate ID numbers for this evaluation
clear model DM CV % clear previous configuration parameters

% Set parameters for model evaluation.

model.factor = [2:5];  

DM.group = datasets;  
       
CV.n = 100;  
CV.method = 'K-fold';
CV.group_based = 1; % perform cross-validation with group information considered

OS.filename = 'test';      
 
[Model,EVA,Log,Data] = fmatrix(data,model,DM,CV,OS);  

% Display the results.
parameter = [ ]; % display evaluation metrics only

plot_fmatrix(Model,EVA,Log,Data,parameter);  

% Based on the above configuration, resample groups numbered 1 and 2. 
% In each loop, draw 5 data points from each of the two groups. 
% In practice, if the sample size of some samples is significantly 
% larger than other samples, this operation can be performed.

DM.resample_subset = [1 2]; % specify the groups to be resampled (e.g., [1 2] --> resample groups numbered 1 and 2); 
                            % the group information is in DM.group.
DM.resample_method = 'random'; % specify the sampling method
                                  % 'random' --> randomly select N data points
                                  % 'K-fold' --> randomly select 1/K data points  
                                  % 'Bootstrap' --> generate an equal number of datasets, allowing for duplication of elements  

DM.resample_N = [5 5]; % randomly select several data points from each group (e.g., [5 5] --> draw 5 data points from each of the two groups);
                       % this parameter is only meaningful (randomly select N data points) when DM.resample_method = 'random'

[Model,EVA,Log,Data] = fmatrix(data,model,DM,CV,OS);  