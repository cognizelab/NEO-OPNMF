function S = plot_biplot(loadings,scores,S)
%% Generate a biplot of the coefficients in the matrix coefs.
%---------------------------------------------------------------------------------------------------------------------------------------------------%
% - Z.K.X. 2022/11/09
%---------------------------------------------------------------------------------------------------------------------------------------------------%
%% Input
%  (1) loadings: N * 2 loadings for the latent variables to be represented as vectors in the biplot               
%  (2) scores: N * 2 scores for each observation to be represented as dots in the biplot               
%  (3) S: Extra settings (optional).
%         --- S.line_color, line colors (Custom RGB colormap as N * 3 matrix, default = [0.3 0.3 0.3]) 
%         --- S.line_label, line labels 
%         --- S.group, vector with a group number for each observation (default, all observations in the same group)
%         --- S.group_center, centers for each group (n group * 2)
%         --- S.group_center_form, cell array of cells with formatting options for each group center
%         --- S.dot_size, dot size (default = 100;)
%         --- S.dot_color, line colors (Custom RGB colormap as Nx3 matrix, default = cbrewer2('Set2',n group)) 
%         --- S.dot_degree, N * n group (0-1), which indicates the degree of membership of each data point within each group
%         --- S.dot_degree_color, custom RGB colormap as N * 3 matrix in different cells
%         --- S.dot_form, cell array of cells with formatting options for each group of dots
%         --- S.ellipses, add 95% confidence ellipses (default, 1)   
%         --- S.density_plot, density plot (default, 0)   
%         --- S.edge_color, edge color (default, [0.5 0.5 0.5])   
%% Example
% -------------------------- Case 1 --------------------------
% clc; clear;
% load fisheriris;
% data = meas;
% labels = {'Sepal length' 'Sepal width' 'Petal length' 'Petal width'};
% group = species;
% [loadings, scores] = pca(zscore(data));
% 
% loadings = loadings(:,1:2);
% scores = scores(:,1:2);
% S.group = group;
% S.line_label = labels;
% 
% plot_biplot(loadings,scores,S);
% -------------------------- Case 2 --------------------------
% % Data
% clc; clear;
% load fisheriris;
% data = meas;
% [loadings, scores] = pca(zscore(data));
% loadings = loadings(:,1:2);
% scores = scores(:,1:2);
% 
% % Fuzzy C-Means Clustering
% N = 3;
% exp = 2;
% maxIter = 100;
% minImprove = 0.00001;
% displayObjective = false;
% options = [exp maxIter minImprove displayObjective];
% 
% [C,U] = fcm(scores,N,options);
% 
% maxU = max(U);
% for i=1:N
%     index = find(U(i,:) == maxU);
%     group(index,1) = i;
% end
% 
% % Plot
% S.group = group;
% S.group_center = C;
% S.dot_degree = U';
% S.line_label = {'A' 'B' 'C' 'D'};
% 
% plot_biplot(loadings,scores,S);
%---------------------------------------------------------------------------------------------------------------------------------------------------%

%% default settings
if nargin < 2, scores = []; end
if nargin < 3 || isempty(S)
    S = struct;
end

if ~isfield(S,'density_plot') || isempty(S.density_plot)
    if ~isempty(scores) & size(scores,1) >= 10000
        S.density_plot = 1;
    else
        S.density_plot = 0;
    end
end
if ~isfield(S,'color') || isempty(S.color)
    S.line_color = repmat([0.3 0.3 0.3],size(loadings,1),1);
end
if ~isfield(S,'ellipses') || isempty(S.ellipses)
    S.ellipses = 1;
end
if ~isfield(S,'group') || isempty(S.group)
    try
        S.group = ones(1,size(scores,1));
    catch
        S.group = ones(1,size(loadings,1));
    end
end
if ~isfield(S,'dot_size') || isempty(S.dot_size)
    S.dot_size = 100;
end
if ~isfield(S,'dot_color') || isempty(S.dot_color)
    S.dot_color = cbrewer2('Set2',numel(unique(S.group)));
end
if ~isfield(S,'edge_color') || isempty(S.edge_color)
    S.edge_color = [0.5 0.5 0.5];
end
if ~isfield(S,'dot_form') || isempty(S.dot_form)
    for i = 1:numel(unique(S.group))
        S.dot_form{i} = {'o','filled','MarkerFaceColor',S.dot_color(i,:),'MarkerEdgeColor',[0.9 .9 .9],'LineWidth',0.3};
    end
end
if ~isfield(S,'group_center_form') || isempty(S.group_center_form)
    for i = 1:numel(unique(S.group))
        S.group_center_form{i} = {'+','filled','MarkerFaceColor',S.dot_color(i,:),'MarkerEdgeColor',[0.3 .3 .3],'LineWidth',2};
    end
end
if ~isfield(S,'line_label') || isempty(S.line_label)
    for i = 1:size(loadings,1)
        S.line_label{i} = num2str(i);
    end
end
if isfield(S,'dot_degree') & ~isempty(S.dot_degree)
    if ~isfield(S,'dot_degree_color') | isempty(S.dot_degree_color)
        k = cbrewer2('Set2',8);
        for i = 1:size(k,1)
            S.dot_degree_color{i} = k(i,:);
        end
    end
end

%% axis
if ~isempty(get(gca,'children'))
    figure();
end
axis square; 
hold on;

% set axis
maxload = max(max(abs(loadings)));
[hax, h0] = formataxis(maxload);
hlab = labelaxis();
    
%% scatter plots
if ~isempty(scores)   
    scores = zscore(scores);
    maxlen = max(sqrt(loadings(:,1).^2 + loadings(:,2).^2));
    if ~isfield(S,'group_center') || isempty(S.group_center)
        scores = scores./max(max(abs(scores)))*maxlen;
    else
        lc = size(S.group_center,1);
        curr = [S.group_center;scores];
        curr = curr./max(max(abs(curr)))*maxlen;
        S.group_center = curr(1:lc,:); curr(1:lc,:) = [];
        scores = curr;        
    end
    groupID = unique(S.group);  
    for igroup = 1:length(groupID)
        if iscell(groupID)
            idx = strcmp(groupID(igroup),S.group);
        else
            idx = (S.group == groupID(igroup));
        end
        % plot dots
        if isfield(S,'density_plot') && S.density_plot == 1
            density_plot(scores,idx,S,igroup);
        else
            if isfield(S,'dot_degree') && ~isempty(S.dot_degree)
                f = find(idx==1);
                x = color2value(S.dot_degree_color{igroup},S.dot_degree(f,igroup),10);
                scatter(scores(idx,1),scores(idx,2),S.dot_size,x,'filled','o','MarkerEdgeColor',S.edge_color,'LineWidth',0.5);
            else
                scatter(scores(idx,1),scores(idx,2),S.dot_size,S.dot_form{igroup}{:},'MarkerEdgeColor',S.edge_color,'LineWidth',0.5);
            end
        end
        % plot ellipses
        if S.ellipses == 1
            %# substract mean
            Mu = mean(scores(idx,:) );
            X0 = bsxfun(@minus,scores(idx,:),Mu);
            %# eigen decomposition [sorted by eigen values]
            [V D] = eig( X0'*X0 ./ (sum(idx)-1) );  %#' cov(X0)
            [D order] = sort(diag(D), 'descend');
            D = diag(D);
            V = V(:, order);
            t = linspace(0,2*pi,100);
            e = [cos(t) ; sin(t)];        %# unit circle
            VV = V*sqrt(D);               %# scale eigenvectors
            e = bsxfun(@plus, VV*e, Mu'); %#' project circle back to orig space
            plot(e(1,:), e(2,:),'Color',[0.5 0.5 0.5],'LineWidth',2);
            fill(e(1,:),e(2,:),S.dot_color(igroup,:),'FaceAlpha',0.3);
        end 
        if isfield(S,'group_center') && ~isempty(S.group_center)
            scatter(S.group_center(igroup,1),S.group_center(igroup,2),100,S.group_center_form{igroup}{:});
        end
    end
end

%% line plots
if ~isempty(loadings)
    % line
    for i = 1:size(loadings,1)
        ht = draw_arrow([0 0],loadings(i,:),0.15,S.line_color(i,:));  
    end
    % label
    labelpoints(loadings(:,1), loadings(:,2),S.line_label);
end

box off
 
end
%% ----------------------------------- END -----------------------------------
%%
function out = draw_arrow(startpoint,endpoint,headsize,color)
% by Ryan Molecke
% accepts two [x y] coords and one double headsize

   v1 = headsize*(startpoint-endpoint)/2.5;

   theta = 22.5*pi/180;
   theta1 = -1*22.5*pi/180;
   rotMatrix = [cos(theta)  -sin(theta) ; sin(theta)  cos(theta)];
   rotMatrix1 = [cos(theta1)  -sin(theta1) ; sin(theta1)  cos(theta1)];
   
   v2 = v1*rotMatrix;
   v3 = v1*rotMatrix1;
   x1 = endpoint;
   x2 = x1 + v2;
   x3 = x1 + v3;
   hold on;
   fill([x1(1) x2(1) x3(1)],[x1(2) x2(2) x3(2)],color);     % this fills the arrowhead (black)
   out = plot([startpoint(1) endpoint(1)],[startpoint(2) endpoint(2)],'linewidth',2,'color',color);
end

%%
function h = abline(m,b,varargin) 
% plots y=m*x+b line behind other plots 
% example: abline(1,0,'Color',[.8 .8 .8]) % reference y = x line
xlim = get(gca,'Xlim');
hExist = get(gca,'children');
hLine = line(xlim, b+m*xlim, varargin{:});
uistack(hExist,'top');
if (nargout>0), h = hLine; end
end

function [hax, h0] = formataxis(l)
% setup axis for biplot
% round limit up (5% margin)
l = l*1.05;
% set boundaries
axis([-l l -l l]) % axis('equal')
hax = gca;
% add origin lines
hx = abline(0,0, 'Color',[.1 .1 .1], 'LineStyle', ':');
hy = vline(0, 'Color',[.1 .1 .1], 'LineStyle', ':');
h0 = [hx hy];
% prettify 
set(hax, ... 
    'FontName'    , 'Helvetica', ... 
    'FontSize'    , 18  , ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'     , ...
    'YMinorTick'  , 'off'     , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...
    'LineWidth'   , 1.2         );
%  'YGrid'       , 'on'      , ...
box off

end

function h = labelaxis()
% label axis 
hx = xlabel('Component 1');
hy = ylabel('Component 2');
hlab = [hx hy];
set(hlab, ...
    'FontName'   , 'Helvetica', ...
    'FontSize'   , 20  );
if (nargout>0), h = hlab; end
end

function h = labelpoints(x, y, labels, color) 
% add labels to points in a scatter plot
if nargin < 4 || isempty(color), color = [.1 .1 .1]; end
% if labels = 0 or [] or {} do nothing
if ~ (isempty(labels) || (~iscell(labels) && isequal(labels, 0)))   
    % if labels are numbers convert to cell array of strings
    if isnumeric(labels), labels = cellstr(num2str(labels')); end
    % add labels to plot
    hold on
    htxt = text(x, y, labels);
    set(htxt, ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
        'FontSize', 18, 'FontName', 'Helvetica', 'Color', color )
    hold off
else
    htxt = [];
end
if (nargout>0), h = htxt; end
end

function format = parseformat(input, default)
% check if input is cell of cells
% if only one group is given it could be a simple cell e.g. {'Color','r'}
% but the parser expects cell of cells, so fix it
if ~iscell(input{1}) % simple cell
    input = {input};
end

% go over input (if less groups are given in format than in groups, the
% remaining will keep default formatting)
format = default;
for igroup = 1:numel(input)
    for pair = reshape(input{igroup}, 2, []) % pair is {Name;Value}
        % is there a default value?
        idx = find(strcmp(pair{1}, default{igroup}));
        if isempty(idx) % add value
            format{igroup} = {format{igroup}{:} pair{1} pair{2}};
        else % replace default value
            format{igroup}{idx+1} = pair{2};
        end
    end
end
end

function opts = parseopts(input, defaults)
% input:    cell array of name-value pairs
% defaults: struct with default values
% opts:     struct with parsed options

opts = defaults;
validNames = fieldnames(defaults);

% go over input pairs
for pair = reshape(input, 2, []) % pair is {Name;Value}
    
    % find option position (case insensitive, accepts substrings)
    name = pair{1};
    match = strncmpi(name, validNames, length(name));
    
    % replace defaults with supplied inputs
    if any(match)
        opts.(validNames{match}) = pair{2};
   else
      error('%s is not a recognized parameter name', name)
   end
end
end

function h = vline(x,varargin) 
% plots vertical line behind other plots crossing x-axis at x 
% example: vline(0,'Color',[.8 .8 .8]) % y axis line
ylim = get(gca,'Ylim');
hExist = get(gca,'children');
hLine = line([x x], ylim, varargin{:});
uistack(hExist,'top');
if (nargout>0), h = hLine; end
end

%%
function density_plot(scores,idx,S,igroup)

disp(['Density mapping for group ',num2str(igroup)]);

K = -1:0.02:1;

f = find(idx==1);
data = scores(f,:);
m = zeros(numel(K) - 1,numel(K) - 1);

for i = 1:numel(K) - 1
    for j = 1:numel(K) - 1
        if i == 1 | j == 1
            f = find(data(:,1) >= K(i) & data(:,1) <= K(i+1) & data(:,2) >= K(j) & data(:,2) <= K(j+1));
        else
            f = find(data(:,1) > K(i) & data(:,1) <= K(i+1) & data(:,2) > K(j) & data(:,2) <= K(j+1));
        end
        m(j,i) = numel(f); % disp([i j]);
    end
end
        
m = rescale(flipud(m));

[a,b] = find(m>0);

k = K; k(k==0) = [];

for i = 1:numel(a)
    value(i,1) = m(a(i),b(i));
end

x = color2value(S.dot_degree_color{igroup},value,10);
scatter(k(b),k(100-a),100,x,'o','filled');  

end

%%
function [value_color,n] = color2value(color,value,n)

% color: RGB N * 3 OR 1 * 3
% value: value to be mapped
% n: number of sections

if size(color,1) > 1
    n = size(color,1);
    a = [];
else
    if nargin < 3 | isempty(n)
        n = 10;
    end        
    a = 0:1/(n):1; a(1) = [];
end

if isempty(a)
    color = color;
else
    p = color;
    for i = 1:length(a)
        color(i,:) = p + ([1 1 1] - p)*( 1 - a(i)/1); 
    end
end

k = 100/(n);
k = [k:k:100];
for i = 1:numel(k)
    edge(i) = prctile(value,k(i));
end

edge = [0,edge];

for i = 1:length(edge) - 1
    f = find(value > (edge(i)) & value <= (edge(i+1)));        
    value_color(f,:) = repmat(color(i,:),numel(f),1);
end

f = find(mean(value,2)==0 | isnan(mean(value,2)));
value_color(f,:) = repmat([1 1 1],numel(f),1);
end 