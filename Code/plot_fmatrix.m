function out = plot_fmatrix(Model,EVA,Log,Data,parameter)
%% Visualization of Matrix Factorization & Evaluation
%% -------------------------------------------------------------------------------------------------------------------------------
%% Input
%  (1-4) Model, EVA, Data & Log: outputs of 'fmatrix.m'
%  (5) paratmter: parameters for displaying factor loading
%      [Figure 01: coefficient in the basis matrix]
%      - parameter.coefficient - itentify specific factors, whose loading matrices will be displayed
%      - parameter.direction - itentify the criteria for item categorization ('abs', based on the maximal absolute value of loadings;'positive', based on the maximal value of loadings)
%      - parameter.items - cell or table array of strings to describe each item
%      - parameter.itemcolor - RGB colors to identify dimensions of items
%      - parameter.itemsize - font size of itemms (default = 10)
%      - parameter.dimlabel - cell array string to identify dimensions of items
%      - parameter.flabel - cell array string to identify factors
%      - parameter.MV - item variation index
%      - parameter.ID - dimension ID number, which is equal to value in "Log.DM.dimension"
%      - parameter.rorder - numerical value (double) for the adjustment of the display order of factors
%      - parameter.cat - if display factor loadings based on different categories (default = 'n')
%      [Figure 02: Evaluation Measures]
%      - parameter.evaluation - itentify specific factors, whose valuation indexes will be displayed
%      - parameter.layout - layout style (1 = <2 * 2>; 2 = <1 * 4>)
%      - parameter.outlier - if show outliers (1 = yes; 0 = no)
%      - parameter.kdata - display the first k data points (default = whole data points)
%% Output
%  (1) out
%         - out.NM: reorganized factor loading
%         - out.NI: reorganized items (serial number)
%         - out.ML: maximal loading
%         - out.parameter: log for parameters
%% Dependency
%  cbrewer2 (https://www.mathworks.com/matlabcentral/fileexchange/58350-cbrewer2)
%% -------------------------------------------------------------------------------------------------------------------------------
% - Z.K.X. 2021/06/12 (MATLAB R2021a)
%% -------------------------------------------------------------------------------------------------------------------------------

%% Factor Loading
if nargin < 2
    EVA = [];
end
if nargin < 3
    Log = [];
end
if nargin < 4
    Data = [];
end

if nargin < 5
    parameter.coefficient = [];
    parameter.show_factor2 = [];
    parameter.direction = 'abs';
    parameter.evaluation = Log.model.factor;
    parameter.show_factor4 = [];
    parameter.ID = [];
end

if ~isfield(parameter,'coefficient')
    parameter.coefficient = [];
else
    if isempty(Log)
        Log.model.factor = [min(parameter.coefficient):max(parameter.coefficient)];
        Log.DM.dimension = parameter.ID;
    end
end

if ~isfield(parameter,'direction')
    parameter.direction = 'abs';
end
if ~isfield(parameter,'show_factor2')
    parameter.show_factor2 = [];
end    
if ~isfield(parameter,'evaluation') 
    parameter.evaluation = Log.model.factor;
end
if ~isfield(parameter,'show_factor4')
    parameter.show_factor4 = [];
end        
if ~isfield(parameter,'itemcolor')
    parameter.itemcolor = [];
end
if ~isfield(parameter,'show_type')
    parameter.show_type = 1;
end     

%% Figure 1: Loading Matrix (coefficient)
if ~isempty(parameter.coefficient) && ~isfield(Model,'group_level_RMSE_permutation_average') 
    if isfield(Model,'group_level_W_cat') 
        if isfield(parameter,'cat') & strcmp(parameter.cat,'y')
            Model.group_level_W = Model.group_level_W_cat;
        end
    end
    for show = 1:length(parameter.coefficient)
        for g = 1:size(Model.group_level_W,1)
            f = find(Log.model.factor==parameter.coefficient(show));
            fac = parameter.coefficient(show);
            FL = Model.group_level_W{g,f}; 
            ID = Log.DM.dimension;
            if isempty(ID)
                if ~isfield(parameter,'ID')
                    error("parameters.ID is required!");
                else
                    ID = parameter.ID;
                end
            end
            if isfield(parameter,'items') & ~isempty(parameter.items)
                items = parameter.items;
            else
                for i = 1:size(Model.group_level_W{1},1)
                    items{i,1} = ['Item-',num2str(i,'%02d')];
                end
            end
            if isfield(parameter,'reorder') & ~isempty(parameter.reorder)
                rorder = parameter.reorder;
            else
                rorder = [];
            end
            if isfield(parameter,'itemsize') & ~isempty(parameter.itemsize)
                itemsize = parameter.itemsize;
            else
                itemsize = 10;
            end                      
            if isfield(parameter,'flabel') & ~isempty(parameter.flabel)
                flabel = parameter.flabel;
            else
                flabel = [];
            end            
            itemcolor = parameter.itemcolor;
            if ~isempty(EVA) & isfield(EVA,'MV')
                if isfield(parameter,'MV') & ~isempty(parameter.MV)
                    SI = parameter.MV(:,fac);
                else
                    c = mean(EVA.MV,3);
                    SI = c(:,fac);
                end
            else
                if isfield(parameter,'MV') & ~isempty(parameter.MV)
                    SI = parameter.MV(:,fac);
                else
                    SI = ones(size(FL,1),1);
                end                    
            end
            if isfield(parameter,'dimlabel')
                dimlabel = parameter.dimlabel;
            else
                dimlabel = [];
            end
            s = ['Group',num2str(g),'_F',num2str(fac)]; 
            figure('name',s);
            [newmatrix,showitem,ml] = plot_factor_loading(FL,SI,ID,items,rorder,flabel,itemcolor,dimlabel,itemsize,parameter.direction);
            out.NM{g,show} = newmatrix; out.NI{g,show} = showitem; out.ML{g,show} = ml;
        end
    end
end

%% Figure 2: Cross-validation Evaluation (Box Plots)
if isfield(EVA,'aRand') && ~isempty(parameter.evaluation) && length(parameter.evaluation) < 9 
    show_factor = parameter.evaluation;
    for i = 1:numel(show_factor)
        LABEL{i} = num2str(show_factor(i));
    end
    if ~isfield(parameter,'layout') || isempty(parameter.layout)
        parameter.layout = 1;
    end    
    if ~isfield(parameter,'outlier') || isempty(parameter.outlier)
        parameter.outlier = 0;
    end 
    if ~isfield(parameter,'kdata') || isempty(parameter.kdata)
        parameter.kdata = Log.CV.n;
    end    
    if parameter.layout == 1
       figure('Position',[100,100,1100,900]);
    elseif parameter.layout == 2
        figure('Position',[100,100,1900,430]);
    end
    kdata = parameter.kdata;    
    if parameter.layout == 1
        subplot(2,2,1);
    elseif parameter.layout == 2
        subplot(1,4,1);
    end
    if parameter.outlier == 1
        X = repmat(show_factor,kdata,1); X = X(:);
        Y = EVA.aRand([1:kdata],show_factor); Y = Y(:);
    else
        Y = []; X =[];
        for i = 1:numel(show_factor)
            y = rmoutliers(EVA.aRand([1:kdata],show_factor(i)),'quartiles');
            Y = [Y;y];
            X = [X;ones(numel(y),1)*i];
        end
    end
    boxchart(categorical(X),Y,'Notch','on');
    set(gca,'XTickLabel',LABEL,'FontSize',12); 
    hold on;
    plot(median(EVA.aRand([1:kdata],show_factor)),'-o','LineWidth',1.5,'MarkerSize',8);
    legend(["Data","Median"],'Location','southeast'); legend('boxoff');
    xlabel('Factors','FontSize',15);
    ylabel('Adjusted Rand Similarity Coefficient','FontSize',15);
    if parameter.layout == 1
        subplot(2,2,2);
    elseif parameter.layout == 2
        subplot(1,4,2);
    end
    if parameter.outlier == 1
        X = repmat(show_factor,kdata,1); X = X(:);
        Y = EVA.CI([1:kdata],show_factor); Y = Y(:);
    else
        Y = []; X =[];
        for i = 1:numel(show_factor)
            y = rmoutliers(EVA.CI([1:kdata],show_factor(i)),'quartiles');
            Y = [Y;y];
            X = [X;ones(numel(y),1)*i];
        end        
    end
    boxchart(categorical(X),Y,'Notch','on');
    set(gca,'XTickLabel',LABEL,'FontSize',12);
    hold on;
    plot(median(EVA.CI([1:kdata],show_factor)),'-o','LineWidth',1.5,'MarkerSize',8);
    legend(["Data","Median"],'Location','southeast'); legend('boxoff');
    xlabel('Factors','FontSize',15);
    ylabel('Concordance Index','FontSize',15);
    if parameter.layout == 1
        subplot(2,2,3);
    elseif parameter.layout == 2
        subplot(1,4,3);
    end
    if parameter.outlier == 1
        X = repmat(show_factor,kdata,1); X = X(:);
        Y = EVA.VI([1:kdata],show_factor); Y = Y(:);
    else
        Y = []; X =[];
        for i = 1:numel(show_factor)
            y = rmoutliers(EVA.VI([1:kdata],show_factor(i)),'quartiles');
            Y = [Y;y];
            X = [X;ones(numel(y),1)*i];
        end        
    end
    boxchart(categorical(X),Y,'Notch','on');
    set(gca,'XTickLabel',LABEL,'FontSize',12);
    hold on;
    plot(median(EVA.VI([1:kdata],show_factor)),'-o','LineWidth',1.5,'MarkerSize',8);
    legend(["Data","Median"],'Location','southeast'); legend('boxoff');
    xlabel('Factors','FontSize',15);
    ylabel('Variation of Information','FontSize',15);
    if parameter.layout == 1
        subplot(2,2,4);
    elseif parameter.layout == 2
        subplot(1,4,4);
    end
    if parameter.outlier == 1
        X = repmat(show_factor,kdata,1); X = X(:);
        Y = EVA.increased_SED([1:kdata],show_factor); Y = Y(:);
    else
        Y = []; X =[];
        for i = 1:numel(show_factor)
            y = rmoutliers(EVA.increased_SED([1:kdata],show_factor(i)),'quartiles');
            Y = [Y;y];
            X = [X;ones(numel(y),1)*i];
        end        
    end
    boxchart(categorical(X),Y,'Notch','on');
    set(gca,'XTickLabel',LABEL,'FontSize',12);
    hold on;
    plot(median(EVA.increased_SED([1:kdata],show_factor)),'-o','LineWidth',1.5,'MarkerSize',8);
    legend(["Data","Median"],'Location','southeast'); legend('boxoff');
    xlabel('Factors','FontSize',15);
    ylabel('Increased Reconfiguration Error','FontSize',15);
    h = sgtitle('Model Stability & Generalizability');
    set(h,'FontSize',20,'FontWeight','normal');    
end
 
%
out.parameter = parameter;

end
%% coefficient
function [newmatrix,showitem,ml] = plot_factor_loading(FL,SI,ID,items,rorder,flabel,itemcolor,dimlabel,itemsize,direction)
%% Input
%  (1) FL: factor loading matrix (items * factors)
%  (2) SI: secondary indicators:  
%  (3) ID: item dimension (serial number)
%  (4) items: item labels
%  (5) rorder: rearrangement of factors (serial number)
%  (6) flabel: factor label
%  (7) itemcolor: item color
%  (8) dimlabel: labels for theoretic item dimensions 
%  (9) itemsize: front size of tiem labels
%  (10) direction: criteria for item categorization 
 
if ~iscell(items)
    label = table2cell(items);
else
    label = items;
end

if ~isempty(rorder)
    newmatrix = FL(:,rorder);
else
    newmatrix = FL;
end

for i = 1:size(newmatrix,1)
    if strcmp(direction,'abs')
        ml(i,1) = find(abs(newmatrix(i,:))==max(abs(newmatrix(i,:))));
    elseif strcmp(direction,'positive')
        ml(i,1) = find(newmatrix(i,:)==max(newmatrix(i,:)));
    end
end

X = []; showitem = [];
for i = 1:length(unique(ml))
    f = find(ml==i); x = newmatrix(f,:);
    if strcmp(direction,'abs')
        [a,b] = sort(abs(x(:,i)),'descend');
    elseif strcmp(direction,'positive')
        [a,b] = sort(x(:,i),'descend');
    end
    x = x(b,:);
    X = [X;x];
    k = f(b);
    showitem = [showitem;k];
end

%%
ytick = 1:size(newmatrix,1);
xtick = 1:size(newmatrix,2);

if ~isempty(flabel)
    xticklabel = flabel;
else
    for i = 1:length(xtick)
        xticklabel{i} = num2str(xtick(i));
    end
end

if isempty(itemcolor)
    try
        itemcolor = cbrewer2('Set2',length(unique(ID)));
    catch
        itemcolor = generate_color_cat;
    end
end

if isempty(itemsize)
    itemsize = 10;
end

if nargin < 8 || isempty(dimlabel)
    for i = 1:length(unique(ID))
        dimlabel{i,1} = num2str(i);
    end
end
    
f1 = subplot('Position',[0.48,0.05,0.02,0.9]);
imagesc(ID(showitem));
colormap(f1,itemcolor); 
set(f1,'YTick',ytick,'YTickLabel',label(showitem),'XTick',[],'XTickLabel',[],'FontSize',itemsize);
xlabel('Items','FontSize',10,'Rotation',30);

f2 = subplot('Position',[0.5,0.05,0.05,0.9]);
stem(SI(showitem),'filled'); xlim([0.5 size(newmatrix,1)+0.5]); ylim([0 1]);
set(gca,'XTick',[],'FontSize',10);
view(90,90);
ylabel('Variation');

f3 = subplot('Position',[0.55,0.05,0.2,0.9]);
imagesc(X); 
set(gca,'YTick',[],'XTick',xtick,'XTickLabel',xticklabel,'FontSize',10); 
colorbar; 
xlabel('Factor');

try
    colormap(f3,flipud(cbrewer2('Spectral',50))); 
catch
    colormap(f3,generate_color_50); 
end

F = 0;
for i = 1:size(newmatrix,2) - 1
    line([i+0.5 i+0.5],[0 size(newmatrix,1)+0.5],'LineWidth',2,'LineStyle','-','Color','black'); 
    f = find(ml==i); F = F+numel(f);
    line([0 size(newmatrix,1)],[F+0.5 F+0.5],'LineWidth',2,'LineStyle','-','Color','black'); 
end

f4 = subplot('Position',[0.8,0.4,0.015,0.15]);
imagesc([1:length(unique(ID))]');
set(gca,'YTick',[1:length(unique(ID))],'YTickLabel',dimlabel,'XTick',[],'FontSize',12,'YAxisLocation','right'); 
colormap(f4,itemcolor);

newmatrix = X; [~,ml] = max(X, [], 2);

end

%% Color
function color_curr = generate_color_50()
color_curr = [0.368627450980392	0.309803921568628	0.635294117647059
            0.332902962548228	0.357721391577151	0.662056358721288
            0.297410156499750	0.403625043100944	0.685496577686323
            0.262190498809226	0.448560808248045	0.706300992018127
            0.227253292838168	0.493158051920390	0.725190871451502
            0.192935898134102	0.537843600479998	0.742750968131035
            0.194917463625426	0.586164422285955	0.742074164445842
            0.238534540078604	0.637581524678404	0.719736058071074
            0.294417778161078	0.687928657669767	0.687785439646511
            0.352244600735673	0.732580099971520	0.658881765806444
            0.412141038529684	0.766623612699109	0.646137917746237
            0.471160671311987	0.792663020566297	0.643958573071747
            0.527897955864135	0.814598814961382	0.644445362723125
            0.582579112579655	0.834093012647500	0.645610790696595
            0.635147044405028	0.852907942384945	0.645350435367962
            0.685818965271568	0.872862404052358	0.641133513206310
            0.739150480235605	0.894283510896355	0.628193274137407
            0.792378256834256	0.915760968253365	0.610733399555541
            0.841269745626912	0.935771299381019	0.595787791521905
            0.882384514228659	0.952734869482390	0.590910411630564
            0.913164829929038	0.965639048892307	0.604010373240172
            0.938253226604679	0.978184038872055	0.637118155506859
            0.959436484971286	0.989677093394554	0.679555277371867
            0.977267622337847	0.998033622257167	0.719449292662263
            0.992758944259555	1	0.745131436748527
            1	0.996091982777108	0.744767771139698
            1	0.979123685064294	0.714366179891001
            1	0.953289252191960	0.665187646744921
            1	0.922535647727882	0.609689085101958
            0.997985926415667	0.890669761508102	0.560320574599493
            0.994481038739963	0.859385127175737	0.524456694930908
            0.994186645738301	0.823329455530606	0.488802813479127
            0.995036348602311	0.783089803491765	0.453116664779688
            0.995233286539529	0.740063303478644	0.419150215181151
            0.993321362908456	0.695696526020084	0.388626398684196
            0.989844055548408	0.649301715379478	0.361676483457733
            0.987508070014149	0.596770827105191	0.335201127201337
            0.983499347950103	0.540993028012321	0.309935250193174
            0.975251549955637	0.485835612321005	0.286760910295060
            0.960691141463968	0.436191051972805	0.266418034331519
            0.941289731179181	0.394297645467696	0.256777872871590
            0.921700868432271	0.354745948681751	0.266622680528681
            0.899673120062352	0.316920673415453	0.285103654804381
            0.872832685669552	0.280744705060415	0.302401525669884
            0.839113445393712	0.246446678861223	0.309905351180865
            0.799644534080146	0.212027150558578	0.306652708472718
            0.757946714564209	0.173816603741874	0.299847504690544
            0.714035973978013	0.130623234728162	0.289581176133295
            0.667914767073927	0.0784367485504265	0.275887896758510
            0.619607843137255	0.00392156862745120	0.258823529411765];
end

function color_curr = generate_color_cat()
            color_curr = [222 236 246;175 200 225;226 242 205;182 218 167;249 213 213;239 152 161;251 227 192;251 201 154;...
                232 224 239; 194 177 215]/255;
end