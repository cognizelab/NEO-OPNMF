function [out,question1,question2] = plot_permutation(Model,EVA,Log,Data,parameter)
%% Visualization of Matrix Factorization & Evaluation
%% -------------------------------------------------------------------------------------------------------------------------------
%% Input
%  (1-4) Model, EVA, Data & Log: outputs of 'fmatrix.m'
%  (5) paratmter: parameters for displaying factor loading
%      - parameter.wh - width and height (default = [500 500]) 
%      - parameter.factor - itentify specific factors, whose permutation results will be displayed (e.g., [2 5]) 
%      - parameter.category - basis for grouping (e.g., 'Age' or 'Gender') 
%      - parameter.separate - separate the two figures of two questions ('y' = yes) 
%      - parameter.color - RGB colors for the two questions
%% Output
%  (1) out
%  (2) question1
%  (3) question2
%% -------------------------------------------------------------------------------------------------------------------------------
% - Z.K.X. 2023/03/20 (MATLAB R2022b)
%% -------------------------------------------------------------------------------------------------------------------------------

%% Factor Loading
if nargin < 5
    parameter.wh = [500 500];
    parameter.factor = 2;
    parameter.category = [];
    parameter.color = [0.988235294117647,0.552941176470588,0.384313725490196;0.552941176470588,0.627450980392157,0.796078431372549];
else
    if ~isfield(parameter,'wh') || isempty(parameter.wh) || numel(parameter.wh) < 2
        parameter.wh = [500 500];
    end
    if ~isfield(parameter,'factor') || isempty(parameter.factor)  
        parameter.factor = 2;
    end    
    if ~isfield(parameter,'category')  
        parameter.category = [];
    end     
    if ~isfield(parameter,'color') || isempty(parameter.color)  
        parameter.color = [0.988235294117647,0.552941176470588,0.384313725490196;0.552941176470588,0.627450980392157,0.796078431372549];
    end        
end

if ~isempty(Log.group_info)
    category = unique(Log.group_info.Group);
else
    category = unique(Log.DM.group);
end

Q1 = Model.group_level_SED_wg_cat - Model.group_level_SED_wg_permutation;
Q2 = Model.group_level_SED_cat - Model.group_level_SED_permutation;

for m = 1:numel(parameter.factor)
    % Q1
    question1(m).Factor = parameter.factor(m);

    figure('Position',[100,100,parameter.wh]);
    f = find(Log.model.factor==(parameter.factor(m)));
    for c = 1:numel(category)
        X1(:,c) = Q1(c,f,:);
    end
    for c = 1:numel(category)
        real_SED(:,c) = Model.group_level_SED_wg_cat(c,f,:);
    end
    for c = 1:numel(category)
        null_SED(:,c) = Model.group_level_SED_wg_permutation(c,f,:);
    end
    question1(m).real_SED = real_SED;
    question1(m).null_SED = null_SED;
    question1(m).real_null_SED = X1;

    bx1 = boxchart(X1,'Notch','on');
    bx1.BoxFaceColor = parameter.color(1,:); bx1.WhiskerLineColor = parameter.color(1,:); bx1.MarkerColor = parameter.color(1,:);   

    line([1-0.5,numel(category)+0.5],[0 0],'LineWidth',1.5,'LineStyle','--','Color',[0 0 0]); hold on;   
    if isfield(Log.DM,'group_label') & ~isempty(Log.DM.group_label)
        set(gca,'xticklabel',Log.DM.group_label,'FontSize',12);
    end
    ylabel('Reconfiguration Error (real-null)','FontSize',15);

    for i = 1:size(X1,2)
        [h,p,ci,stat] = ttest(X1(:,i)); 
        t1(i,1) = stat.tstat; p1(i,1) = p;
    end
    out(m).t1 = t1;
    out(m).p1 = p1;
    sig1 = zeros(1,size(X1,2));
    [a,b] = makeFDR(p1,0.05);
    sig1(p1<=a) = 1;
    out(m).bad1 = find(sig1==1&t1'>0)'; 
    out(m).good1 = find(sig1==1&t1'<0)';  

    % Q2
    question2(m).Factor = parameter.factor(m);

    if isfield(parameter,'separate') & strcmp(parameter.separate,'y')
        title(['Group-level ',num2str(parameter.factor(m)),'-factor Model'],'FontSize',20);  
        figure('Position',[100,100,parameter.wh]);      
    else        
        hold on
    end
    for c = 1:numel(category)
        X2(:,c) = Q2(c,f,:);
    end
    for c = 1:numel(category)
        real_SED(:,c) = Model.group_level_SED_cat(c,f,:);
    end
    for c = 1:numel(category)
        null_SED(:,c) = Model.group_level_SED_permutation(c,f,:);
    end
    question2(m).real_SED = real_SED;
    question2(m).null_SED = null_SED;
    question2(m).real_null_SED = X2;    
    
    bx2 = boxchart(X2,'Notch','on');
    bx2.BoxFaceColor = parameter.color(2,:); bx2.WhiskerLineColor = parameter.color(2,:); bx2.MarkerColor = parameter.color(2,:);   

    if isfield(parameter,'separate') & strcmp(parameter.separate,'y')
        line([1-0.5,numel(category)+0.5],[0 0],'LineWidth',1.5,'LineStyle','--','Color',[0 0 0]);  
        if isfield(Log.DM,'group_label') & ~isempty(Log.DM.group_label)
            set(gca,'xticklabel',Log.DM.group_label,'FontSize',12);
        end
        ylabel('Reconfiguration Error (real-null)','FontSize',15);
        if isempty(parameter.category)
            title(['Category-specific ',num2str(parameter.factor(m)),'-factor Model'],'FontSize',20);  
        else
            title([parameter.category,'-specific ',num2str(parameter.factor(m)),'-factor Model'],'FontSize',20);  
        end
    else
        title([num2str(parameter.factor(m)),'-factor Model'],'FontSize',20);  
    end
    for i = 1:size(X2,2)
        [h,p,ci,stat] = ttest(X2(:,i)); 
        t2(i,1) = stat.tstat; p2(i,1) = p;
    end    
    out(m).t2 = t2;
    out(m).p2 = p2;    
    sig2 = zeros(1,size(X2,2));
    [a,b] = makeFDR(p2,0.05);
    sig2(p2<=a) = 1;
    out(m).bad2 = find(sig2==1&t2'>0)'; 
    out(m).good2 = find(sig2==1&t2'<0)';  
end

end
%%
function [pID,pN,h] = makeFDR(p,q)
pp = p;
p = sort(p(:));
V = length(p);
I = (1:V)';

cVID = 1;
cVN = sum(1./(1:V));

pID = p(find(p<=I/V*q/cVID, 1, 'last' ));
pN  = p(find(p<=I/V*q/cVN,  1, 'last' ));

if isempty(pID) ~= 1
    h = find(pp<=pID);
else
    h = [];
end
return

end