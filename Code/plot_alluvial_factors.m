function S = plot_alluvial_factors(pos_data,colors_map,factors,x_labels,y_labels,line_width,line_alpha,varargin)
%%  
%------------------------------------------------------------------------------------------------------------------------------%
% - Z.K.X. 2023/03/28
%------------------------------------------------------------------------------------------------------------------------------%
%% Input
%  (1) pos_data: N of items * N of variants matrix
%                  <specify the relative position of items in different variants>
%  (2) colors_map: N of items * 3 RGB color schemes for coloring different items
%                  <default = [0.5 0.5 0.5] for all items>
%  (3) factors: factor ID （N of items * N of variants matrix）
%                  <specify the factors/dimensions for different items across variants>
%  (4) x_labels: cell array of labels for variants (X-axis) 
%  (5) y_labels: cell array of labels for items (Y-axis) 
%  (6) line_width: line width
%                  <default = 2, for all lines>
%  (7) line_alpha: line transparency
%                  <default = 1, for all lines, i.e., non-transparent>
%% Output
%  (1) S: setting parameters
%% Extra Input
%  (1) figure_width: figure width 
%                  <default = [600 800]>
%  (2) figure_high: figure high 
%                  <default = [600 800]>
%  (3) top_bottom: order of items
%                  <'y', from top to bottom (default)>
%  (4) line_style: line style
%                  <1, curve lines (default); 2, straight line>
%  (5) dot_style: dot style
%                  <default = 'o'>
%  (6) dot_size: dot size
%                  <default = 150>
%  (7) x_FontSize: font size for X-axis 
%  (8) y_FontSize: font size for Y-axis 
%------------------------------------------------------------------------------------------------------------------------------%
%% Example   
% clc; clear;
% pos_data = [[1:10]',[1:10]',[4 1 7 6 2 5 3 9 8 10]'];
% colors_map = [repmat([0.5, 0.5, 0.5],3,1);repmat([0.2, 0.6, 0.3],3,1);repmat([0.8, 0.2, 0.7],4,1)];
% 
% factors = [[1;1;1;2;2;2;3;3;3;3],[1;1;1;1;2;2;2;3;3;3],[1;1;1;2;2;2;2;3;3;3]];
% 
% x_labels = {'S1','S2','S3'};
% y_labels = {'Item1','Item2','Item3','Item4','Item5','Item6','Item7','Item8','Item9','Item10'};
% 
% line_width = 5;
% line_alpha = 1
%;
% S = plot_alluvial_factors(pos_data);
% S = plot_alluvial_factors(pos_data,colors_map,factors,x_labels,y_labels,line_width,line_alpha,'line_style',1,'dot_style','s','dot_size',200,'category_spacing',2);
% S = plot_alluvial_factors(pos_data,colors_map,factors,x_labels,y_labels,line_width,line_alpha,'top_bottom','n','line_style',2,'dot_style','s','category_spacing',2);
% 
% factors = [1 1 1 2 2 2 3 3 3 3]';
% S = plot_alluvial_factors(pos_data,colors_map,factors,x_labels,y_labels,'line_style',1,'dot_size',200,'category_spacing',1,'line_width',5);
%------------------------------------------------------------------------------------------------------------------------------%

%% Default Setting
if nargin < 2 || isempty(colors_map) || size(colors_map,1) ~= size(pos_data,1)
    colors_map = repmat(0.5,size(pos_data,1),3);
end
if nargin < 3 || isempty(factors) || size(factors,1) ~= size(pos_data,1)  
    factors = ones(size(pos_data,1),size(pos_data,2));
elseif size(factors,2) ~= size(pos_data,2)  
    factors = repmat(factors,1,size(pos_data,2));
end

if (nargin >= 8)
    idx = 1;
    while idx <= (nargin - 7)
        switch varargin{idx}
            case {'figure_width'}
                figure_width = varargin{idx+1}; 
            case {'figure_high'}
                figure_high = varargin{idx+1};                 
            case {'top_bottom'}
                top_bottom = varargin{idx+1};                    
            case {'line_style'}
                line_style = varargin{idx+1};
            case {'dot_style'}
                dot_style = varargin{idx+1};   
            case {'dot_size'}
                dot_size = varargin{idx+1};                   
            case {'category_spacing'}
                category_spacing = varargin{idx+1};                 
            case {'x_FontSize'}
                x_FontSize = varargin{idx+1};                              
            case {'y_FontSize'}
                y_FontSize = varargin{idx+1};                    
        end
        idx = idx + 1;
    end
else
    figure_width = 600;
    figure_high = 800;
    top_bottom = 'y';
    line_style = 1;
    dot_style = 'o';
    dot_size = 150;
    category_spacing = 1;
    x_FontSize = 12;
    y_FontSize = 12;
end

if ~exist ('figure_width','var')  
    figure_width = 600;
end
if ~exist ('figure_high','var') 
    figure_high = 800;
end
if ~exist ('top_bottom','var')
    top_bottom = 'y';
end
if ~exist ('line_style','var')
    line_style = 1;
end
if ~exist ('dot_style','var')
    dot_style = 'o';
end
if ~exist ('dot_size','var')
    dot_size = 100;
end
if ~exist ('category_spacing','var')
    category_spacing = 1;
end
if ~exist ('x_FontSize','var')
    x_FontSize = 12;
end
if ~exist ('y_FontSize','var')
    y_FontSize = 12;
end
if ~exist ('line_width','var') || isempty(line_width)
    line_width = ones(size(pos_data, 1),1)*2;
else 
    if numel(line_width) == 1
        line_width = repmat(line_width,size(pos_data, 1),1);
    end
end
if ~exist ('line_alpha','var') || isempty(line_width)
    line_alpha = ones(size(pos_data, 1),1);
else 
    if numel(line_alpha) == 1
        line_alpha = repmat(line_alpha,size(pos_data, 1),1);
    end
end

S.figure_wh = [figure_width,figure_high]; S.top_bottom = top_bottom; S.line_style = line_style; S.dot_style = dot_style; S.dot_size = dot_size; 
S.category_spacing = category_spacing; S.line_width = line_width; S.line_alpha = line_alpha;

if exist ('top_bottom','var') & strcmp(top_bottom,'y')
    pos_data = max(pos_data) + 1 - pos_data; 
    factors = flip(factors); 
end

%% Data Manipulation
num_stages = size(pos_data, 2);

for j = 1:size(factors,2)
    gap_row = [];
    for i = 2:size(factors,1)
        if factors(i,j) ~= factors(i-1,j)
            gap_row = [gap_row;i];
        end
    end
    GAP{j} = gap_row;
end

pos_data2 = pos_data;

for j = 1:size(factors,2)
    gap_row = GAP{j};
    for i = 1:numel(gap_row)
        if i == numel(gap_row)
            f = find(pos_data(:,j)>=gap_row(i));
            pos_data2(f,j) = pos_data2(f,j)+i*category_spacing;            
        else
            f = find(pos_data(:,j)>=gap_row(i)&pos_data(:,j)<gap_row(i+1));
            pos_data2(f,j) = pos_data2(f,j)+i*category_spacing;
        end
    end
end

%% Main Figure
figure('Position',[50,50,figure_width,figure_high ]);

hold on

for i = 1:size(pos_data2, 1)
    color = colors_map(i,:);           
    for j = 1:num_stages - 1
        if line_style == 1
            x = [j, j + 0.5, j + 1];
            y = [pos_data2(i, j), mean([pos_data2(i, j), pos_data2(i, j + 1)]) + 0.5 * abs(pos_data2(i, j) - pos_data2(i, j + 1)), pos_data2(i, j + 1)];
            plot_quadratic_bezier(x, y, color, line_width(i), line_alpha(i));
        elseif line_style == 2
            control_points = [j, pos_data2(i, j); (j + j + 1) / 2, (pos_data2(i, j) + pos_data2(i, j + 1)) / 2; j + 1, pos_data2(i, j + 1)];
            plot_bezier(control_points, color, line_width(i), line_alpha(i));
        end
        scatter(j, pos_data2(i, j), dot_size, color, 'filled', dot_style, 'MarkerEdgeColor', 'white');
    end
    scatter(num_stages, pos_data2(i, num_stages), dot_size, color, 'filled', dot_style, 'MarkerEdgeColor', 'white');
end
hold off

%% Labels
set(gca, 'XTick', 1:num_stages);

if exist('x_labels') & ~isempty(x_labels)
    set(gca, 'XTickLabel',x_labels,'FontSize',x_FontSize);
    S.x_labels = x_labels;
else
    set(gca, 'XTickLabel', arrayfun(@(x) sprintf('Variant %d', x), 1:num_stages, 'UniformOutput', false),'FontSize',x_FontSize);
end

ylim([0 max(max(pos_data)) + 1]);

if exist ('y_labels','var') & iscell(y_labels)
    if size(y_labels,1) ~= 1
        y_labels = y_labels';
    end   
    S.y_labels = y_labels;
    if exist ('top_bottom','var') & strcmp(top_bottom,'y')
        y_labels = fliplr(y_labels);
    end    
    gap_row = GAP{1};
    if ~isempty(gap_row)
        for i = 1:numel(gap_row)
            cells_to_insert{i} = repmat({''},1,category_spacing);
        end            
        y_labels = insert_cells(y_labels, cells_to_insert, gap_row - 1);
    end
    set(gca, 'YTick', 1:length(y_labels), 'YTickLabel', y_labels); ax = gca; set(ax.YAxis, 'FontSize', y_FontSize);
end    

box off

end

%% Quadratic Bessel Curves
function plot_quadratic_bezier(x, y, color, line_width, line_alpha)

t = linspace(0, 1, 100);
Bx = (1 - t).^2 * x(1) + 2 * (1 - t) .* t * x(2) + t.^2 * x(3);
By = (1 - t).^2 * y(1) + 2 * (1 - t) .* t * y(2) + t.^2 * y(3);
 
plot(Bx, By, 'LineWidth', line_width, 'Color', [color,line_alpha]);

end

%% Bessel Curves
function plot_bezier(points, color, line_width, line_alpha)

n = size(points, 1);
t = linspace(0, 1);
B = zeros(length(t), 2);

for i = 1:length(t)
    for j = 1:n
        B(i, :) = B(i, :) + points(j, :) * nchoosek(n-1, j-1) * t(i)^(j-1) * (1 - t(i))^(n - j);
    end
end
 
plot(B(:, 1), B(:, 2), 'LineWidth', line_width, 'Color', [color,line_alpha]);

end

%% Cell Interpolation 
function new_cell_array = insert_cells(original_cell_array, cells_to_insert, insert_positions)

num_inserts = numel(insert_positions);
assert(numel(cells_to_insert) == num_inserts, 'Number of cells to insert must match number of insert positions.');
new_cell_array = original_cell_array;
 
for i = 1:num_inserts
    num_cells_inserted_before = sum(cellfun(@numel, cells_to_insert(1:i-1)));
    insert_pos = insert_positions(i) + num_cells_inserted_before;
    if insert_pos >= numel(new_cell_array)
        new_cell_array = [new_cell_array, cells_to_insert{i}];
    else
        new_cell_array = [new_cell_array(1:insert_pos), cells_to_insert{i}, new_cell_array(insert_pos+1:end)];
    end
end

end