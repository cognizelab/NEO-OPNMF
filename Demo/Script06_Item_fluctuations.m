%% --------- Fluctuation of Items' Affiliations --------- 
% This script is used to display the volatility of the same items under different allocation schemes.

%% Toolkit configuration

% Set up path for core codes.

clc; clear;  
base_dir = pwd;

addpath(genpath(fullfile(base_dir,'Code')));

%% Data loading

% Load demo data.

load(fullfile(base_dir,'Data','data_alluvial.mat')); 

%% Visualization

% Generate item alluvial plots.

S = plot_alluvial_factors(pos_data,colors_map,factors,x_labels,y_labels,line_width,line_alpha,...
    'figure_width',1000,'figure_high',1300, ...
    'line_style',1,'dot_style','s','dot_size',80,'category_spacing',3,'y_FontSize',6);

% Refer to the documentation by entering doc plot_alluvial_factors to access more detailed information about the input options and their functionalities.