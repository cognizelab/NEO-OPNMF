function [Model,EVA,Log,Data] = fmatrix(data,model,DM,CV,OS)
%% Matrix Factorization & Evaluation
%% -------------------------------------------------------------------------------------------------------------------------------
%% Input
%  (1) data: subjects * items matrix  
%  (2) model: model configuration   
%               ---  model.group_level <whole group-level model generation>
%                                 [1] 0 = NO
%                                 [2] 1 = YES (default)
%                                 [3] 2 = only execute the group-level model generation
%               ---  model.fmatrix <matrix factorization method>
%                                 [1] 'OPNMF' - Orthogonal Projective Non-negative Matrix Factorization (default)
%                                 [2] 'MLE' - Maximum Likelihood Estimate (Factor Analysis)
%                                 [3] 'PAF' - Principal Axis Factoring (Factor Analysis)
%                                 [4] 'PCA' - Principal Component Analysis
%               ---  model.factor <factor setting>
%                                one row of vector array which identifies the number of factos (default = factor = [2:1:9])     
%               ---  model.outlier <outlier filtering>
%                                  model.outlier(1) = the number of scaled MAD (default = 3)
%                                  model.outlier(2) = identify a proportion X, that subjects who are detected as 
%                                                     outliers in X proportion of factors would be removed (default = 0.5)       
%                                  (model.outlier = [] means non-execution)
%      ---------------------------------------------------------------------------------------------------------------------------
%      ||Note, 'model' setting can also simply be a pre-defiend model, which means you'd like to apply this model to new sample.|| 
%      ||Note, if 'model.outlier' is set to the empty set, the program will skip this step. ||
%      ---------------------------------------------------------------------------------------------------------------------------
%  (3) DM: data manipulation   
%               ---  DM.subID <subject ID>
%               ---  DM.missing <handling data with NaN values>
%                                 [1] 'remove' (default)
%                                 [2] 'mode'
%               ---  DM.group <class identification> 
%                              a column vector of classified (e.g., family relationship, study batch, study site, gender, culture\country, etc) 
%                              or continuous (e.g., age) data labels
%               ---  DM.group_label <class identification, which is useful only when DM.group is double array data> 
%               ---  DM.demographics <demographic information corresponding to the data> 
%               ---  DM.dimension <theoretical dimensions>
%               ---  DM.filter <data filtering>
%                                 [1] 'std'  
%                                 [2] value of upper MAD £¨should be > 1£©
%                                 [3] value of percentiles £¨should be > 0 & < 1£©
%                                 [4] empty ([]) - non-execution (default)
%               ---  DM.class <data classification, i.e., the grouping basis for conducting permutation tests>
%                              'none'\'cat'\'bin','window'
%               ---  DM.binB <boundaries or edges of bins for data classification, e.g., [20:5:40]>
%               ---  DM.windowL <length of windows for data classification>
%               ---  DM.windowW <step width of windows for data classification>
%               ---  DM.permutation <number of permutations for class, 0 = non-execution> 
%               ---  DM.permutation_save <'y' save random models; 'n' not save randonm models (default)>   
%               ---  DM.resample_subset <subset ID to be resampled>  
%               ---  DM.resample_method <resampling method>
%                                 [1] 'random'  
%                                 [2] 'K-fold'
%                                 [3] 'Bootstrap'
%               ---  DM.resample_N <number of subjects to be saved in different resampled subset>
%               ---  DM.resample_retain <number of runs a single resampling to be retained>
%  (4) CV: cross validation setting (no execution, if it is set to the empty set)
%               ---  CV.method <resampling strategy>  
%                                 [1] 'K-fold'  
%                                 [2] 'Bootstrap'
%               ---  CV.K <K folds>    
%               ---  CV.n <number of cross validation>  
%               ---  CV.direction <Generalization Direction>  
%                                 [1] 'one-way' (within-sample or within-dataset cross validation) 
%                                 [2] 'two-way' (within-sample or within-dataset cross validation) 
%                                 [3] 'mixed'  (cross-sample or cross-dataset cross validation)   
%               ---  CV.group_based <if execute group-based cross validation> 
%                                 [1] 1 = YES
%                                 [2] 0 = NO
%               ---  CV.st <source and target groups of generalization direction if 'mixed' strategy is selected>   
%                                || This setting allows a n * 2 matrix, indicating the mixed generalization directions
%                                   (e.g., [1 2; 1 3] means training a model in group 1 then generalize it to group 2 and 3).||
%               ---  CV.permutation <number of cross validation>    
%  (5) OS: other optional setting
%               || setting for OPNMF ||
%               ---  OS.opnmf_max_iter <maximum number of iterations in OPNMF (default = 50000)>
%               ---  OS.opnmf_tol <convergence tolerance in OPNMF (default = 1e-6)>
%               || setting for MLE (same as parameters used in matlab function 'factoran')||
%               ---  OS.fa_scores <'wls' (default) | 'regression'>
%               ---  OS.fa_start <'Rsquared' (default) | 'random' | Positive integer | Matrix with d rows>
%               ---  OS.fa_rotate <'varimax' (default) | 'none' | 'quartimax' | 'equamax' | 'parsimax' | 'orthomax' | 'promax' | 'procrustes' | 'pattern' | function handle>
%               ---  OS.fa_delta <0-1, default = 0.005> 
%               || setting for PAF (same as parameters used in matlab function 'rotatefactors')||
%               ---  OS.fa_rotate <'varimax' (default)>
%               || setting for PCA (same as parameters used in matlab function 'pca')||
%               ---  OS.pca_algorithm <'svd' (default) | 'eig' | 'als'>
%               ---  OS.pca_centered <false (default) | true>
%               ---  OS.pca_options <[] (default) | struct input of 'Option' in 'pca'>
%               ---  OS.pca_VariableWeights <[] (default) | 'variance'> 
%               ---  OS.fa_rotate <'varimax' (default)>
%               || outputs||
%               ---  OS.filename <output filename> 
%% Output
%  (1) Model: model outputs
%               ---  Model.data_dr, dimensionality reduced data based on a pre-trained model (basis matrix)
%               ---  Model.group_level_W \ Model.group_level_W_cat 
%               ---  Model.group_level_H \ Model.group_level_H_cat
%               ---  Model.group_level_RMSE \ Model.group_level_RMSE 
%               ---  Model.group_level_SED \ Model.group_level_SED_cat
%               ---  Model.individual_level_RMSE \ Model.individual_level_RMSE_cat
%               ---  Model.individual_level_SED \ Model.individual_level_SED_cat
%               ---  Model.outliers
%               ---  Model.outliers_removed_ID
%               ---  Model.outliers_removed_position
%  (2) EVA: evaluation indexes
%               --- EVA.Dice, Dice Coefficient            
%               --- EVA.Jacard, Jaccard Index            
%               --- EVA.NMI, Normalized Mutual Information      
%               --- EVA.Rand, Rand Similarity Coefficient   
%               --- EVA.zRand, Z-score of Rand Similarity Coefficient   
%               --- EVA.aRand, adjusted Rand Similarity Coefficient   
%               --- EVA.CI, Concordance Index
%               --- EVA.VI, Variation of Information   
%               --- EVA.MV, Modular Variability (for detecting the unstable items in specific factor model) 
%               --- EVA.increased_RMSE, Increased Root Mean Square Error (RMSE)    
%               --- EVA.increased_SED, Increased Squared Euclidean distance  
%  (3) Log: log of the model evaluation parameters 
%  (4) Data:data fed into the model training
%% Dependency
%  Brainlets (https://github.com/asotiras/brainparts)
%% Reference
% - Chen, J., Patil, K. R., Weis, S., Sim, K., Nickl-Jockschat, T., Zhou, J., ... & Habel, U. (2020). 
%   Neurobiological Divergence of the Positive and Negative Schizophrenia Subtypes Identified on a New Factor Structure of 
%   Psychopathology Using Non-negative Factorization: An International Machine Learning Study. Biological psychiatry, 87(3), 282-293.
%% -------------------------------------------------------------------------------------------------------------------------------
% - Z.K.X. 2023/03/20 (MATLAB R2018a)
%% -------------------------------------------------------------------------------------------------------------------------------

%% Default Setting  
if nargin < 2 || isempty(model)  
    model.fmatrix = 'OPNMF';
    model.factor = [2:9];
    model.outlier = [];
    model.group_level = 1;
end
if  ~isa(model,'double') 
    if ~isfield(model,'group_level') 
        model.group_level = 1;
    end
    if ~isfield(model,'factor')  
        model.factor = [2:9];
    end
    if ~isfield(model,'fmatrix')  
        model.fmatrix = 'OPNMF';
    end    
    if ~isfield(model,'outlier')  
        model.outlier = [];
    elseif ~isempty(model.outlier)
        if length(model.outlier) < 2
            model.outlier(2) = 0.5;
        end
    end
end
% ---------------------------------
if nargin < 3 || isempty(DM)  
    DM.subID  = [1:size(data,1)]';
    DM.missing  = 'remove';
    DM.group  = [1:size(data,1)]';
    DM,dimension = [];
    DM.filter = [];
    DM.class = [];
    DM.permutation = 0;
    DM.resample_subset = [];
end
if ~isfield(DM,'subID') || isempty(DM.subID)
    DM.subID  = [1:size(data,1)]';
end
if ~isfield(DM,'resample_subset')
    DM.resample_subset  = [];
end
if ~isfield(DM,'missing')
    DM.missing  = 'remove';
end
if ~isfield(DM,'group')
    DM.group  = [1:size(data,1)]';
else
    if iscell(DM.group)
        g = unique(DM.group);
        DM.group_label = g;
        for i = 1:length(g)
            f = find(strcmp(DM.group,g{i}));
            c(f,1) = i;
        end
        DM.group = c;
        clear c g;
    end
end
if ~isfield(DM,'dimension')
    DM.dimension  = [];
end
if ~isfield(DM,'permutation')
    DM.permutation  = [];
end
if ~isfield(DM,'demographics')
    DM.demographics = [];
end   
if ~isempty(DM.resample_subset) 
    if ~isfield(DM,'resample_method') || isempty(DM.resample_method)
        DM.resample_method  = 'random';
    end
    if ~isfield(DM,'resample_N')
        DM.resample_N = [];
    end
    if ~isfield(DM,'resample_retain')
        DM.resample_retain = 1;
    end    
end
% ---------------------------------
if nargin < 4 
    CV.method = 'K-fold';
    CV.K = 5;    
    CV.n = 5000;
    CV.direction = 'one-way';
    CV.group_based = 1;
end
if ~isempty(CV)
    if ~isfield(CV,'method')
        CV.method = 'K-fold';
    end
    if ~isfield(CV,'K')
        CV.K = 5;
    end
    if ~isfield(CV,'K')
        CV.K = 5;
    end
    if ~isfield(CV,'n')
        CV.n = 1000;
    end    
    if ~isfield(CV,'direction')
        CV.direction = 'one-way';
    end
    if ~isfield(CV,'group_based')
        CV.group_based = 1;
    end
    if ~isfield(CV,'permutation')
        CV.permutation = 0;
    end
end
% ---------------------------------
if nargin < 5 || isempty(OS)
    OS.opnmf_max_iter = 5000;
    OS.opnmf_tol = 1e-6;
end
if ~isfield(OS,'opnmf_max_iter')
    OS.opnmf_max_iter = 5000;
end
if ~isfield(OS,'opnmf_tol')
    OS.opnmf_tol = 1e-6;
end
% ---------------------------------
%% Model Identification
if ~isa(model,'double') 
    if strcmp(model.fmatrix,'OPNMF')
        if find(data<0)
          error('The entered data cannot contain negative numbers.')
        end
        model_expression = "[W,H] = opnmf(X,K,[],[],OS.opnmf_max_iter,OS.opnmf_tol);";
        disp('Model: Orthogonal Projective Non-negative Matrix Factorization');
    elseif strcmp(model.fmatrix,'MLE')
        model_expression = "[W,H] = factor_analysis(X,K,OS);";
        disp('Model: Maximum Likelihood Estimate');
    elseif strcmp(model.fmatrix,'PAF')
        model_expression = "[W,H] = principal_axis_factoring(X,K);";
        disp('Model: Principal Axis Factoring');        
    elseif strcmp(model.fmatrix,'PCA')
        model_expression = "[W,H] = principal_component_analysis(X,K,OS);";
        disp('Model: Principal Component Analysis');
    end
end

%% Applying Pre-defined Model to New Data
if  isa(model,'double') 
    Model.data_dr = data*model;
    Model.group_level_RMSE = rmse(model,data');
    Model.group_level_SED = sed(model,data');
    for sub = 1:size(data,1)
        Model.individual_level_RMSE(sub,1) = rmse(model,data(sub,:)');
        Model.individual_level_SED(sub,1) = sed(model,data(sub,:)');
    end    
    return
end

%% Handling Data with Missing Values   
Log.data.inputN = size(data,1);

if ~isempty(DM.missing) 
    [a,b] = find(isnan(data));
    Log.data.missingN = length(unique(a));
    if strcmp(DM.missing,'remove')    
        if ~isempty(a)
            Log.removed_data_ID = DM.subID(a);
            Log.removed_data_Position = a;
            data(a,:) = [];
            DM.subID(a,:) = [];
            DM.group(a,:) = [];
            if ~isempty(DM.demographics)
                DM.demographics(a,:) = [];
            end
        end
    else
        if strcmp(DM.missing,'mode')
            v = mode(data);
        end
        for i = 1:length(a)
            data(a(i),b(i)) = v(b(i));
        end
    end
end

%% Data Filtering Based on Theoretical Dimension
if isfield(DM,'filter') && ~isempty(DM.dimension)      
    for i = 1:numel(unique(DM.dimension))
        f = find(DM.dimension==i);
        within(:,i) = std(data(:,f)');   
    end    
    v = mean(within,2)./std(data')';
    if isstr(DM.filter) & strcmp(DM.filter,'std')
        f = find(v>=1 | isnan(v));
    elseif DM.filter >= 1
        [TF,L,U,C] = isoutlier(v,'ThresholdFactor',DM.filter); 
        f = find(v>=U);
    elseif DM.filter > 0 & DM.filter < 1
        K = DM.filter * 100;
        K = prctile(v,K);
        f = find(v>=K);
    end
    Log.data.filteringN = length(f);
    Log.std_value = v;
    Log.poor_quality_ID = DM.subID(f);
    data(f,:) = [];
    DM.subID(f,:) = [];
    DM.group(f,:) = []; 
    if ~isempty(DM.demographics)
        DM.demographics(f,:) = [];
    end    
end

%% Group Level Model Building
if isvector(model.factor)
    if model.group_level > 0
        disp('--------- Group Level Model Building ---------');
        for j = 1:length(model.factor)
            K = model.factor(j);
            X = data';
            eval(model_expression);
            Model.group_level_W{j} = W;
            Model.group_level_H{j} = H;
            for k = 1:size(W,1)
                Model.group_level_maximal_loading{j}(k,1) = find(abs(W(k,:))==max(abs(W(k,:))));
            end
            Model.group_level_RMSE(j) = rmse(Model.group_level_W{j},X);
            Model.group_level_SED(j) = sed(Model.group_level_W{j},X);
            for sub = 1:size(data,1)
                Model.individual_level_RMSE{j}(sub,1) = rmse(W,data(sub,:)');
                Model.individual_level_SED{j}(sub,1) = sed(W,data(sub,:)');
            end 
        end
        if ~isempty(model.outlier)
            disp('--------- Outlier Filtering ---------');
            for j = 1:length(model.factor)
                [TF,L,U,C] = isoutlier(Model.individual_level_SED{j},'ThresholdFactor',model.outlier(1)); 
                temp = zeros(size(data,1),1);
                temp(Model.individual_level_SED{j}>=U) = 1;
                Model.outliers(:,j) = temp;
            end
            temp = sum(Model.outliers,2); f = find(temp>fix(model.outlier(2)*length(model.factor)));
            Model.outliers_removed_ID = DM.subID(f);
            Model.outliers_removed_position = f;
            Log.data.outliersN = length(f);
            data(f,:) = [];
            DM.subID(f,:) = [];
            DM.group(f,:) = [];
            if ~isempty(DM.demographics)
                DM.demographics(f,:) = [];
            end
            for j = 1:length(model.factor)
                K = model.factor(j);
                X = data';
                eval(model_expression);
                Model.group_level_W{j} = W;
                Model.group_level_H{j} = H;
                for k = 1:size(W,1)
                    Model.group_level_maximal_loading{j}(k,1) = find(abs(W(k,:))==max(abs(W(k,:))));
                end
                Model.group_level_RMSE(j) = rmse(Model.group_level_W{j},X);
                Model.group_level_SED(j) = sed(Model.group_level_W{j},X);
                for sub = 1:size(data,1)
                    Model.individual_level_RMSE{j}(sub,1) = rmse(W,data(sub,:)');
                    Model.individual_level_SED{j}(sub,1) = sed(W,data(sub,:)');
                end 
            end        
        end
    end
end

%% Data Classification
if isfield(DM,'class') && ~isempty(DM.class) 
    class = DM.class;
    group = DM.group;
    [dataset,Log.group_info] = classdata(data,class,group,DM);
    if model.group_level ~= 2
        disp('--------- Group Classification ---------');
        if length(find(group == min(unique(group)))) < min(model.factor)
            error('The amount of data in some groups is too small!');
        end
        for i = 1:length(dataset)
            for j = 1:length(model.factor)
                K = model.factor(j);
                X = dataset{i}';
                eval(model_expression);
                Model.group_level_W_cat{i,j} = W;
                Model.group_level_H_cat{i,j} = H;
                for k = 1:size(W,1)
                    Model.group_level_maximal_loading_cat{i,j}(k,1) = find(abs(W(k,:))==max(abs(W(k,:))));
                end
                Model.group_level_SED_cat(i,j) = sed(Model.group_level_W_cat{i,j},X);
                Model.group_level_SED_wg_cat(i,j) = sed(Model.group_level_W{j},X);
                Model.group_level_SED_wg_inc_cat(i,j) = ire(W,Model.group_level_W{j},X);
                for sub = 1:size(dataset{i},1)
                    Model.individual_level_SED_cat{i,j}(sub,1) = sed(W,dataset{i}(sub,:)');
                end 
            end
        end
    end
end

%% Permutation Test for Class\Group
if ~isempty(DM.permutation) & DM.permutation > 0 & model.group_level ~= 2
    disp('--------- Permutation Test for Class\Group ---------');
    for per = 1:DM.permutation
        tStart = tic;
        per_id = randperm(size(data,1));
        Log.permutation_ID(per,:) = per_id;
        per_group = DM.group(per_id);
        per_dataset = classdata(data,DM.class,per_group,DM);        
        for i = 1:length(dataset)
            for j = 1:length(model.factor)
                K = model.factor(j);
                X = per_dataset{i}';  
                eval(model_expression);
                Model.group_level_SED_permutation(i,j,per) = sed(W,per_dataset{i}');
                Model.group_level_SED_wg_permutation(i,j,per) = sed(Model.group_level_W{j},per_dataset{i}');
                Model.group_level_SED_wg_inc_permutation(i,j,per) = ire(W,Model.group_level_W{j},per_dataset{i}');
                if isfield(DM,'permutation_save') & strcmp(DM.permutation_save,'y') 
                    Model.group_level_W_cat_permutation{i,j}(:,:,per) = W;
                    Model.group_level_H_cat_permutation{i,j}(:,:,per) = H;
                end
            end
        end     
        disp(['Permutation Test for Class\Group: ',num2str(per)]);
        tEnd = toc(tStart);
        fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
    end
end

%% Within-sample Cross Validation
if ~isempty(CV) && model.group_level ~= 2
if strcmp(CV.direction,'one-way') || strcmp(CV.direction,'two-way')
    disp('--------- Within-sample Cross Validation ---------');
    if strcmp(CV.method,'K-fold')
        if ~isfield(CV,'K') || isempty(CV.K)         
            CV.K = 5;
        end
    end
    if CV.group_based == 1
        group = DM.group;
    else
        group = [1:size(data,1)];
    end
    seg = unique(group);  
    XX = data';   
    for i = 1:CV.n  
        tStart = tic;
        if ~isempty(DM.resample_subset)
            if i <= DM.resample_retain || mod(i,DM.resample_retain) == 0
                [rdata,group,order] = resample_subgroup(data,DM.group,DM.resample_method,DM.resample_subset,DM.resample_N);
                XX = rdata'; seg = unique(group); 
            end
        end
        for K = model.factor     
            if strcmp(CV.method,'K-fold')
               index = crossvalind('Kfold',numel(seg),CV.K);
               cross = unique(index);
            elseif strcmp(CV.method,'Bootstrap')
               [f,o] = bootstrap(seg);
               cross = 1;
            end          
            for j = 1:numel(cross)         
                if strcmp(CV.method,'K-fold')
                    f = find(index~=j); o = find(index==j);      
                end                
                [is,pos] = ismember(group,seg(f));
                f = find(pos~=0);                    
                [is,pos] = ismember(group,seg(o));
                o = find(pos~=0);  
                X = XX(:,f); V1 = X;
                eval(model_expression);
                W1 = W; H1 = H;
                X = XX(:,o); V2 = X;         
                eval(model_expression);
                W2 = W; H2 = H;
                % -----------------------------------------------------------------------------------------
               [~,P1(:,1)] = yael_kmax(abs(single(W1)'),1);
               [~,P2(:,1)] = yael_kmax(abs(single(W2)'),1);
               if strcmp(CV.direction,'one-way')
                   ire_temp(j) = ire(W2,W1,V2);
                   irmse_temp(j) = rmse(W1,V2) - rmse(W2,V2);
               elseif strcmp(CV.direction,'two-way')
                   ire_temp(j) = (ire(W2,W1,V2)+ire(W1,W2,V1))/2;
                   irmse_temp(j) = ((rmse(W1,V2) - rmse(W2,V2))+(rmse(W2,V1) - rmse(W1,V1)))/2;
               end
               temp = ps(P1,P2,'zrand','mv','jaccard');
               VI(j) = temp.VI;
               Jacard(j) = temp.Jacard;
               zRand(j) = temp.zRand;
               aRand(j) = temp.SAR;
               Rand(j) = temp.SR;
               MV(:,j) = temp.MV;
               CI(j) = ci(W1,W2); 
               % -----------------------------------------------------------------------------------------
               if CV.permutation > 0
                   for per = 1:CV.permutation
                       indx = randperm(size(data,2));
                       [~,P1(:,1)] = yael_kmax(abs(single(W1)'),1);
                       [~,P2(:,1)] = yael_kmax(abs(single(W2(indx,:))'),1);               
                       if strcmp(CV.direction,'one-way')
                           ire_temp_rand(per) = ire(W2(indx,:),W1,V2(indx,:));
                           irmse_temp_rand(per) = rmse(W1,V2(indx,:)) - rmse(W2(indx,:),V2(indx,:));
                       elseif strcmp(CV.direction,'two-way')
                           ire_temp_rand(per) = (ire(W2(indx,:),W1,V2(indx,:))+ire(W1,W2(indx,:),V1))/2;
                           irmse_temp_rand(per) = -1*(((rmse(W2(indx,:),V2(indx,:)) - rmse(W1,V2(indx,:)))+(rmse(W1,V1) - rmse(W2(indx,:),V1)))/2);
                       end
                       temp = ps(P1,P2,'zrand','mv','jaccard');
                       VI_rand(per) = temp.VI;
                       Jacard_rand(per) = temp.Jacard;
                       zRand_rand(per) = temp.zRand;
                       aRand_rand(per) = temp.SAR;
                       Rand_rand(per) = temp.SR;
                       MV_rand(:,per) = temp.MV;
                       CI_rand(per) = ci(W1,W2(indx,:));
                    end
                    ire_temp_rand2(j) = mean(ire_temp_rand);
                    irmse_temp_rand2(j) = mean(irmse_temp_rand);
                    CI_rand2(j) = mean(CI_rand);
                    VI_rand2(j) = mean(VI_rand);
                    Jacard_rand2(j) = mean(Jacard_rand);
                    zRand_rand2(j) = mean(zRand_rand);
                    aRand_rand2(j) = mean(aRand_rand);
                    Rand_rand2(j) = mean(Rand_rand);
                    MV_rand2(:,j) = mean(MV_rand,2);
               end
               % -----------------------------------------------------------------------------------------
            end         
            EVA.increased_SED(i,K) = mean(ire_temp);   
            EVA.increased_RMSE(i,K) = mean(irmse_temp);   
            EVA.VI(i,K) = mean(VI);
            EVA.Jacard(i,K) = mean(Jacard);
            EVA.zRand(i,K) = mean(zRand);
            EVA.aRand(i,K) = mean(aRand);
            EVA.Rand(i,K) = mean(Rand);
            EVA.CI(i,K) = mean(CI);
            EVA.MV(:,K,i) = mean(MV,2);     
            % ----------------------------------
            if CV.permutation > 0
                EVA.increased_SED_random(i,K) = mean(ire_temp_rand2);   
                EVA.increased_RMSE_random(i,K) = mean(irmse_temp_rand2);   
                EVA.VI_random(i,K) = mean(VI_rand2);
                EVA.Jacard_random(i,K) = mean(Jacard_rand2);
                EVA.zRand_random(i,K) = mean(zRand_rand2);
                EVA.aRand_random(i,K) = mean(aRand_rand2);
                EVA.Rand_random(i,K) = mean(Rand_rand2);
                EVA.CI_random(i,K) = mean(CI_rand2);
                EVA.MV_random(:,K,i) = mean(MV_rand2,2);      
            end
            % ----------------------------------
        end           
        disp(['Within-sample Cross Validation: ',num2str(i)]);
        tEnd = toc(tStart);
        fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
    end     
end
end

%% Cross-sample (Class\Group) Cross Validation
if ~isempty(CV) & model.group_level ~= 2
if strcmp(CV.direction,'mixed')
    disp('--------- Cross-sample (Class\Group) Cross Validation ---------');
    if ~isfield(CV,'st') || isempty(CV.st)
        n = numel(unique(DM.group));
        st = [];
        for i = 1:n
            for j = 1:n
                if i~=j
                    st = [st;i,j];
                end
            end
        end 
    else
        st = CV.st;
    end
    if  isempty(DM.resample_subset)
        DM.resample_subset = [1:length(unique(DM.group))];
    end
    for i = 1:CV.n      
        tStart = tic;
        [rdata,group,order] = resample_subgroup(data,DM.group,DM.resample_method,DM.resample_subset,DM.resample_N);
        G = unique(group);
        for n = 1:numel(G)
            if i <2 || (i >2 & ~isempty(find(DM.resample_subset==G(n))))
                % disp(['Models in resampled group-',num2str(G(n)),' data is trainning!']);
                f = find(group==G(n));
                X = rdata'; X = X(:,f); x{n} = X; 
                for K = model.factor
                    eval(model_expression);
                    setW{n,K} = W; setH{n,K} = H; 
                end
            end
        end     
        % -----------------------------------------------------------------------------------------
        for k = model.factor
            clear corrW;
            for j = 1:size(st,1)                
               [~,P1(:,1)] = yael_kmax(abs(single(setW{st(j,1),k})'),1);
               [~,P2(:,1)] = yael_kmax(abs(single(setW{st(j,2),k})'),1);               
               ire_temp(j) = ire(setW{st(j,2),k},setW{st(j,1),k},x{st(j,2)});       
               irmse_temp(j) = rmse(setW{st(j,1),k},x{st(j,2)}) - rmse(setW{st(j,2),k},x{st(j,2)});
               temp = ps(P1,P2,'zrand','mv','jaccard');
               VI(j) = temp.VI;
               Jacard(j) = temp.Jacard;
               zRand(j) = temp.zRand;
               Rand(j) = temp.SR;
               aRand(j) = temp.SAR;
               MV(:,j) = temp.MV;
               CI(j) = ci(setW{st(j,1),k},setW{st(j,2),k});                      
            end
            EVA.increased_SED(i,k) = mean(ire_temp);            
            EVA.increased_RMSE(i,k) = mean(irmse_temp);   
            EVA.VI(i,k) = mean(VI);
            EVA.Jacard(i,k) = mean(Jacard);
            EVA.zRand(i,k) = mean(zRand);
            EVA.aRand(i,k) = mean(aRand);
            EVA.Rand(i,k) = mean(Rand);
            EVA.CI(i,k) = mean(CI);
            EVA.MV(:,k,i) = mean(MV,2);     
        end     
        % -----------------------------------------------------------------------------------------
        if CV.permutation > 0
            for k = model.factor
                for j = 1:size(st,1)       
                    for per = 1:CV.permutation
                        indx = randperm(size(X,1));
                        S = setW{st(j,2),k};
                        S = S(indx,:);
                       [~,P1(:,1)] = yael_kmax(abs(single(setW{st(j,1),k})'),1);
                       [~,P2(:,1)] = yael_kmax(abs(single(S)'),1);    
                       SX = x{st(j,2)};
                       SX = SX(indx,:);
                       ire_temp_rand(j,per) = ire(S,setW{st(j,1),k},SX);       
                       irmse_temp_rand(j,per) = rmse(setW{st(j,1),k},SX) - rmse(S,SX);
                       CI_rand(j,per) = ci(setW{st(j,1),k},S);  
                       temp = ps(P1,P2,'zrand','mv','jaccard');
                       VI_rand(j,per) = temp.VI;
                       Jacard_rand(j,per) = temp.Jacard;
                       zRand_rand(j,per) = temp.zRand;
                       Rand_rand(j,per) = temp.SR;
                       aRand_rand(j,per) = temp.SAR;
                       MV_rand(:,j,per) = temp.MV;
                    end
                end
                ire_temp_rand2 = mean(ire_temp_rand,2);
                irmse_temp_rand2 = mean(irmse_temp_rand,2);
                CI_rand2 = mean(CI_rand,2);
                VI_rand2 = mean(VI_rand,2);
                Jacard_rand2 = mean(Jacard_rand,2);
                zRand_rand2 = mean(zRand_rand,2);
                aRand_rand2 = mean(aRand_rand,2);
                Rand_rand2 = mean(Rand_rand,2);
                MV_rand2 = mean(MV_rand,3);
                EVA.increased_SED_random(i,k) = mean(ire_temp_rand2);            
                EVA.increased_RMSE_random(i,k) = mean(irmse_temp_rand2);   
                EVA.VI_random(i,k) = mean(VI_rand2);
                EVA.Jacard_random(i,k) = mean(Jacard_rand2);
                EVA.zRand_random(i,k) = mean(zRand_rand2);
                EVA.aRand_random(i,k) = mean(aRand_rand2);
                EVA.Rand_random(i,k) = mean(Rand_rand2);
                EVA.CI_random(i,k) = mean(CI_rand2);
                EVA.MV_random(:,k,i) = mean(MV_rand2,2);     
            end                         
        end
        % -----------------------------------------------------------------------------------------
        disp(['Cross-sample (Class\Group) Cross Validation: ',num2str(i)]);
        tEnd = toc(tStart);
        fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
    end
end
end

%%
if ~exist('Model'); Model = []; end
if ~exist('EVA'); EVA = []; end
Log.model = model; Log. DM = DM; Log.CV = CV; Log.OS = OS; Data = data;

if ~isfield(OS,'filename')
    save FMATRIX Model EVA Log Data
else
    eval(['save ',OS.filename,' Model EVA Log Data']);
end
end
%% ----------------------------------------------------------------------------------------------------------------------------------------

%% ----------------------------------------------------------------------------------------------------------------------------------------
function RMSE = rmse(W,V)
%% Root mean square error (RMSE).
% W = items * dimension
% V = items * subjject
RMSE = norm(V-W*((V'*W)'),'fro')/sqrt(size(V,1)*size(V,2)); 
end

function SED = sed(W,V)
%% Reconstruction Error (RE).
% W = items * dimension
% V = items * subjject
SED = mean(abs(sum(abs(round(W*(W'*V))-V),1)));
end

function [rdata,rgroup,order] = resample_subgroup(data,group,method,ID,N)
%%
if length(unique(ID)) < length(ID)
    error('Subset ID to be resampled must be unique!');
end

if isempty(N)
    if strcmp(method,'K-fold')  
        N = 5;
    elseif strcmp(method,'random')   
        for i = 1:length(unique(group))
            k(i) = numel(find(group==i));
        end
        N = ones(1,numel(ID))*fix(mean(k));
        if find(k(ID)<N)
            error("Can not use 'DM.resample_method = random' without an appropriate N of resamples!");
        end
    end
end

if strcmp(method,'random')   
    if length(N) < length(ID)
        N = repmat(N,1,length(ID));
    end
end

order = [];
for i = 1:length(unique(group))
    f = find(group==i);
    [a,b] = ismember(ID,i);        
    if mean(b) > 0    
        b = find(b==1);
        if strcmp(method,'Bootstrap')
            f = bootstrap(f);
            if ~isempty(N)
                f = f(N(b));
            end
        elseif strcmp(method,'K-fold')
            index = crossvalind('Kfold',numel(f),N); 
            f = f(index~=1);     
        elseif strcmp(method,'random')
            index = randsample(numel(f),N(b));
            f = f(index);     
        end
    end
    order = [order;f];    
end

rdata = data(order,:);
rgroup = group(order);
end

function [f,o] = bootstrap(in)
%%
s = size(in);     
if min(s) == 1  
   f = in(ceil(max(s)*rand(max(s),1)));    
else         
   f = in(ceil(s(1)*s(2)*rand(s(1),s(2),1))); 
end
o = setdiff(in,f);
end

function IRE = ire(H1,H2,V1)
%% Increased reconstruction error.    
% H2: S; H1: T
IRE = mean(abs(sum(abs(round(H2*(H2'*V1))-V1),1)-sum(abs(round(H1*(H1'*V1))-V1),1)));  
end

function CI = ci(H1,H2)
%% Concordance index (CI).
NumItem = size(H1,1); % item number         
dim = size(H1,2); % dimension number    
lowloadHB1 = nan(NumItem,dim); lowloadHB2 = nan(NumItem,dim);
lowloadHB1(:,:) = H1; lowloadHB2(:,:) = H2;
for j=1:NumItem
    EDHB1=norm(lowloadHB1(j,:),2);
    EDHall=norm(lowloadHB2(j,:),2);
    EDHB11(j,1)=EDHB1;
    EDHall11(j,1)=EDHall;
end
HB1normalized=bsxfun(@rdivide,lowloadHB1(:,:),EDHB11);
Hallnormalized=bsxfun(@rdivide,lowloadHB2(:,:),EDHall11);
S1 = HB1normalized*HB1normalized';
S2 = Hallnormalized*Hallnormalized';
CI = (1-(norm(S1-S2,'fro')^2)./(NumItem*NumItem-NumItem));
end

%% ----------------------------------------------------------------------------------------------------------------------------------------
function PS = ps(p1,p2,varargin)
%% Partition similarity. 
%% Input
%   (1) p1: partition 1
%   (2) p2: partition 2
%   (3) algorithm: consistency-based algorithms
%% Output
%  ---  dice: Dice Coefficient   
%  ---  jaccard: Jaccard Index    
%  ---  NMI: Normalized Mutual Information
%  ---  zRand: Z-score of Rand Similarity Coefficient
%  ---  SAR: Adjusted Rand Similarity Coefficient    
%  ---  VI: Variation of Information 
%  ---  MV: Modular Variability   
%% -------------------------------------------------------------------------------------------------------------------------------
idx = 1;
while idx <= (nargin - 2)
    switch varargin{idx}
        case {'zrand'}
            [PS.zRand,PS.SR,PS.SAR,PS.VI] = zrand(p1,p2);
        case {'nmi'}
            PS.NMI = normalized_mutual_information(p1,p2);
        case {'jaccard'}
            PS.Jacard = getJaccard(p1,p2);  
        case {'dice'}
            PS.Dice = getDiceCoeff(p1,p2);    
        case {'mv'}
            PS.MV = mv(p1,p2);   
    end
    idx = idx + 1;
end

end

function Jacard = getJaccard(A,B)
%% Jaccard Index    
A = agreement(A);
B = agreement(B);
Jacard = sum(A.*B)/(sum(A+B) - sum(A.*B));
end

function DC = getDiceCoeff(A,B)
%% Dice Coefficient    
A = agreement(A);
B = agreement(B);
DC = 2*(sum(A.*B))/sum(A + B);
end

function [zRand,SR,SAR,VI] = zrand(part1,part2)
%% ZRAND Calculates the z-Rand score and Variation of Information
%  distance between a pair of partitions.
%
%   [zRand,SR,SAR,VI] = ZRAND(part1,part2) calculates the z-score of the
%   Rand similarity coefficient between partitions part1 and part2. The
%   Rand similarity coefficient is an index of the similarity between the
%   partitions, corresponding to the fraction of node pairs identified the
%   same way by both partitions (either together in both or separate in
%   both)
%
%   NOTE: This code requires genlouvain.m to be on the MATLAB path
%
%   Inputs:     part1,  | Partitions that are being
%               part2,  | compared with one another
%
%   Outputs:    zRand,  z-score of the Rand similarity coefficient
%               SR,     Rand similarity coefficient
%               SAR,    Adjusted Rand similarity coefficient
%               VI,     Variation of information
%               
%
%   Amanda L. Traud, Eric D. Kelsic, Peter J. Mucha, and Mason A. Porter,
%   "Comparing Community Structure to Characteristics in Online Collegiate
%   Social Networks," SIAM Review 53, 526-543 (2011).
if size(part1, 1) == 1
    part1 = part1';
end
if size(part2, 1)==1
    part2 = part2';
end
if length(part1) ~= length(part2)
    disp('ERROR: partitions not of equal length')
    return
end
nij = sparse(part1+1, part2+1, 1);
ni = sum(nij, 2);
nj = sum(nij, 1);
n = length(part1);
M = n*(n-1)/2;
M1 = sum(ni.^2-ni)/2;
M2 = sum(nj.^2-nj)/2;
a = full(sum(sum(nij.^2-nij)))/2; %same in both
b = M1-a;                         %same in 1, diff in 2
c = M2-a;                         %same in 2, diff in 1
d = M-(a+b+c);                    %diff in both
SR = (a+d)/(a+b+c+d);
meana = M1*M2/M;
SAR = (a-meana)/((M1+M2)/2-meana);
C1 = 4*sum(ni.^3)-8*(n+1)*M1+n*(n^2-3*n-2);
C2 = 4*sum(nj.^3)-8*(n+1)*M2+n*(n^2-3*n-2);
vara = M/16 - (4*M1-2*M)^2*(4*M2-2*M)^2/(256*M^2) + C1*C2/(16*n*(n-1)*(n-2)) + ...
    ((4*M1-2*M)^2-4*C1-4*M)*((4*M2-2*M)^2-4*C2-4*M)/(64*n*(n-1)*(n-2)*(n-3));
zRand = (a-meana)/sqrt(vara);

c1 = unique(part1);
c2 = unique(part2);
H1 = 0; H2 = 0; I = 0;
for i = c1'
    pi = ni(i+1)/n;
    H1 = H1-pi*log(pi);
    for j = c2'
        if nij(i+1,j+1)
            pj = nj(j+1)/n;
            pij = nij(i+1,j+1)/n;
            I = I+pij*log(pij/pi/pj);
        end
    end
end
for j = c2'
    pj = nj(j+1)/n;
    H2 = H2-pj*log(pj);
end
VI = (H1+H2-2*I);

end

function [nmi] = normalized_mutual_information(partition1,partition2,varargin)
%% NORMALIZED_MUTUAL_INFORMATION - This function computes the normalized mutual
%         information (normalization taken from from Danon 2005, note other normalizations
%         exist) which is a measure of similarity between two different partitions of 
%         community structure.  The measure is bounded in [0 1] where 1 implies perfect
%         agreement between partitions.  In cases where the number of nodes
%         divided by the number of clusters <~ 100, the adjusted value should be
%         used to correct for chance (see Vinh et al 2010).
% 
% Inputs: 
%         partition1 - vector containing community assignment for each node in first partition
% 
%         partition2 - vector containing community assignment for each node in second partition
% 
%         varargin - set equal to 'adjusted' if want to correct for chance and compute the AMI 
%               (use this option if the number of nodes diveded by the number of communities <~ 100.  
%               Default is the unadjusted NMI (set varaible equal to 'unadjusted' or leave empty).
% 
% Outputs: 
%         nmi - normalized mutual information.  Gives measure of
%               similarity between partitons.  can be adjusted for chance 
%               (AMI- adjusted normalized mutual information - see Vihn 2010)
% 
% Written by Sarah Feldt Muldoon

% References: Danon L, Diaz-Guilera A, Duch J, Arenas A (2005) Comparing
%         community structure identification. J Stat Mech:P09008.; Vinh NX, Epps J, 
%         Bailey J (2010) Information Theoretic Measures for Clusterings Comparison: 
%         Variants, Properties, Normalization and Correction for Chance. The Journal
%         of Machine Learning Research.


%if no keyword, calculate the unadjusted NMI
if isempty(varargin)
    varargin{1} = 'unadjusted';
end

num_nodes = length(partition1);
comm1_array = unique(partition1);
comm2_array = unique(partition2);
num_comm1 = length(comm1_array);
num_comm2 = length(comm2_array);
contingency_matrix = zeros(num_comm1,num_comm2);

%make contingency table (also called confusion_matrix)
for i = 1:num_comm1
    comm1 = comm1_array(i);
    for j = 1:num_comm2
        comm2 = comm2_array(j);
        contingency_matrix(i, j) = sum(partition1 == comm1 & partition2 == comm2);
    end
end

sum_over_i = sum(contingency_matrix,2);
sum_over_j = sum(contingency_matrix,1);
total_sum = sum(sum(contingency_matrix));

%compute mutual information
mi_matrix = contingency_matrix .* log((contingency_matrix .* total_sum) ./ (sum_over_i * sum_over_j));
mi = sum(mi_matrix(isfinite(mi_matrix))) / total_sum;

%compute normalization terms (we normalize by (h1+h2)/2 as in Danon 2005)
h1 = sum_over_i .* log(sum_over_i ./ total_sum);
h1 = -sum(h1(isfinite(h1))) / total_sum;
h2 = sum_over_j .* log(sum_over_j ./ total_sum);
h2 = -sum(h2(isfinite(h2))) / total_sum;

if strcmp(varargin{1}, 'adjusted') == 1
    %compute expected mutual information (Vihn 2010)
    expected_mi = 0;
    for i = 1:num_comm1
        for j=1:num_comm2
            k_min = max([sum_over_i(i) + sum_over_j(j) - total_sum, 1]);
            k_max = min([sum_over_i(i), sum_over_j(j)]);
            for k = k_min:k_max
                term1 = k / total_sum;
                term2 = log(total_sum * k / (sum_over_i(i) * sum_over_j(j)));
                term3_numerator = factorial(sum_over_i(i)) * factorial(sum_over_j(j)) * factorial(total_sum - sum_over_i(i)) * factorial(total_sum - sum_over_j(j));
                term3_denominator = factorial(total_sum) * factorial(k) * factorial(sum_over_i(i) - k) * factorial(sum_over_j(j) - k) * factorial(total_sum - sum_over_i(i) - sum_over_j(j) + k);
                emi_term = term1 * term2 * term3_numerator / term3_denominator;
                if isfinite(emi_term)
                    expected_mi = expected_mi + emi_term;
                end
            end
        end
    end
    %compute AMI_sum as in Vihn 2010
    nmi = (mi - expected_mi) / (.5 * (h1 + h2) - expected_mi);
elseif strcmp(varargin{1}, 'unadjusted') == 1
    %compute NMI as in Danon 2005 (NMI_sum in Vihn 2010)
    nmi = 2 * mi / (h1 + h2);
end

end

function [MV,cMV] = mv(N1,N2,ID)
%% Calculation of Modular Variability (MV) for A Given Node
%
%  Z.K.X. 2017/12/15
%-----------------------------------------------------------------------------------%
%   MV: For each node, its module affiliation variability is calculated 
%       between two modular partitions through a metric of modular variability (MV) 
%
%                      |Xk(i)¡ÉXk(j)|   |Xk(i)¡ÉXk(j)|
%       MVk(i,j) = 1 - -------------- * --------------
%                          |Xk(i)|          |Xk(j)|
%       where Xk (i) and Xk (j) denote the module labels to which node k belongs in 
%       modular partitions i and j, respectively. Xk (i)¡ÉXk (j) represents the 
%       common node set between modules Xk (i) and Xk (j), and |Xk(i)¡ÉXk(j)| 
%       denotes the number of nodes in the common node set. For node k between 
%       modular structures i and j, a small overlap between the two modules Xk (i) 
%       and Xk (j) indicates large module affiliation variability.
%  cMV: For each node, the MV values are calculated for every networks identified by 
%       ID. The cMV indicates how possible the community pattern from a given node to 
%       a given network would be changed. 
%-----------------------------------------------------------------------------------%
if (nargin < 3)
	ID = [];
end
if size(N1,2) == 1 
    N1 = agreement(N1);
    N2 = agreement(N2);
end
for i = 1:length(N1)
    M1 = sum(N1(i,:))+1;
    M2 = sum(N2(i,:))+1;
    cover = find(N1(i,:)==1&N2(i,:)==1);
    com = length(cover)+1;
    MV(i,1) = 1-(com/M1)*(com/M2);
end
if ~isempty(ID)
    for i = 1:length(N1)
        for j = 1:max(ID)
            F = find(ID==j);
            M1 = sum(N1(i,F))+1;
            M2 = sum(N2(i,F))+1;
            cover = find(N1(i,F)==1&N2(i,F)==1);
            com = length(cover)+1;
            cMV(i,j) = 1-(com/M1)*(com/M2);
        end
    end
end

end

function D = agreement(ci,buffsz)
%% AGREEMENT Agreement matrix from clusters
%
%   D = AGREEMENT(CI) takes as input a set of vertex partitions CI of
%   dimensions [vertex x partition]. Each column in CI contains the
%   assignments of each vertex to a class/community/module. This function
%   aggregates the partitions in CI into a square [vertex x vertex]
%   agreement matrix D, whose elements indicate the number of times any two
%   vertices were assigned to the same class.
%
%   In the case that the number of nodes and partitions in CI is large
%   (greater than ~1000 nodes or greater than ~1000 partitions), the script
%   can be made faster by computing D in pieces. The optional input BUFFSZ
%   determines the size of each piece. Trial and error has found that
%   BUFFSZ ~ 150 works well.
%
%   Inputs,     CI,     set of (possibly) degenerate partitions
%               BUFFSZ, optional second argument to set buffer size
%
%   Outputs:    D,      agreement matrix
%
%   Richard Betzel, Indiana University, 2012

%modification history
%09.24.2012 - added loop for big N that makes the function slower but also
% prevents it from maxing out memory.

n = size(ci,2);

if nargin < 2
    buffsz = 1000;
end

if n <= buffsz
    
    ind = dummyvar(ci);
    D = ind*ind';
    
else
    
    a = 1:buffsz:n;
    b = buffsz:buffsz:n;
    
    if length(a) ~= length(b)
        b = [b, n];
    end
    
    x = [a' b'];
    nbuff = size(x,1);
    
    D = zeros(size(ci,1));
    for i = 1:nbuff
       y = ci(:,x(i,1):x(i,2));
       ind = dummyvar(y);
       D = D + ind*ind';
    end
    
end

D = D.*~eye(length(D));
end

function [W H] = opnmf(X, K, w0, initMeth, max_iter, tol,iter0, save_step, outputdir)
%% Orthogonal Projective NonNegative Matrix Factorization
% input:
%   X          nonnegative data input (D times N)
%   K          number of components
%   w0         given initialization of W (optional)
%   initMeth   determines which initialization method will be used if no w0
%              is provided
%   max_iter   maximum number of iterations (default 50000)
%   tol        convergence tolerance (default 1e-5)
%   iter0      initial iteration (used when resuming optimization after
%              possible failure - use in combination with saved 
%              intermediate results)
%   save_step  save intermediate results every # number of steps
%   outputdir  directory where intermediate results are saved
%
% output:
%   W          the factorizing matrix (D times K)
%   H          expansion coefficients   
%
% Relevant references:
%
% Sotiras, A., Resnick, S. M., & Davatzikos, C. (2015). Finding imaging 
% patterns of structural covariance via Non-Negative Matrix Factorization. 
% NeuroImage, 108, 1-16.
%
% Yang, Z., & Oja, E. (2010). Linear and nonlinear projective nonnegative 
% matrix factorization. Neural Networks, IEEE Transactions on, 21(5), 
% 734-749.
[Dinit,N] = size(X);
% Basic argument check
if ~isscalar(K) || ~isnumeric(K) || K<1 || K>min(Dinit,N) || K~=round(K)
    error('opnmf:badK','K should be positive integer no larger than the number of rows or columns in X');
end
if ~ismatrix(X) || ~isnumeric(X) || ~isreal(X) || any(any(X<0)) || any(any(~isfinite(X)))
    error('opnmf:badX','opnmf:X must be a matrix of non-negative values.')
end
% get rid of the background in order to boost the computational efficiency
% assumption: background corresponds to positions with stackwise mean value
% equal to zero
mean_im = mean(X,2);
data_matrix_nz = X((mean_im>0),:) ; 
X = data_matrix_nz ; clear data_matrix_nz ;
% variables
check_step = 100;
D = size(X,1);
% initialize w0
if ~exist('w0','var') || isempty(w0)
     %%disp('Initializing w0: ');
    if ~exist('initMeth', 'var') || isempty(initMeth)
        initMeth = 1;
    end
    switch initMeth
        case 0
             %%disp('random initialization ...') ;
            W = rand(D,K);
        case 1 % NNDSVD
             %%disp('NNDSVD initialization ...') ;
            [W,~] = NNDSVD(X,K,0) ;
        case 2 % NNDSVDa
             %%disp('NNDSVDa initialization ...') ;
            [W,~] = NNDSVD(X,K,1) ;
        case 3 % NNDSVDar
             %%disp('NNDSVDar initialization ...') ;
            [W,~] = NNDSVD(X,K,2) ;
        case 4 % NNDSVD using randomized SVD
             %%disp('NNDSVD initialization using random SVD calculation for efficiency ...') ;
            [W,~] = NNDSVD(X,K,3) ;
        otherwise
             %%disp('NNDSVD initialization ...') ;
            [W,~] = NNDSVD(X,K,0) ;
    end
     %%disp('done') ;
else
    W = w0;
    clear w0 ;
end
% check variables
if ~exist('max_iter', 'var') || isempty(max_iter)
    max_iter = 50000;
end
if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-5;
end
if (~exist('iter0','var') || isempty(iter0))
    iter0 = 1;
end
if (~exist('save_step','var') || isempty(save_step) )
    save_step = floor(max_iter/10) ;
end
if (exist('outputdir','var') && ~isempty(outputdir))    
    % then create the directory
    if(~strcmp(outputdir(end),'/'))
        outputdir=[outputdir '/'];
    end
    if(~exist(outputdir,'dir'))
        success = mkdir(outputdir);
        if(~success)
            error('opnmf:BadDir',['Output directory ' outputdir ' can not be created']);
        end
    end
end
% start optimization
XX = X * X';
for iter=iter0:max_iter
    W_old = W;
    if mod(iter,check_step)==0
         %%fprintf('iter=% 5d ', iter);
    end
     % multiplicative update rule   
    W = W .* (XX*W) ./ (W*(W'*XX*W));
    W = W ./ norm(W);
    
    diffW = norm(W_old-W, 'fro') / norm(W_old, 'fro');
    if diffW<tol
         %%fprintf('converged after %d steps.\n', iter);
        break;
    end
    
    if mod(iter,check_step)==0
         %%fprintf('diff=%.10f, ', diffW);
         %%fprintf('obj=%.10f', norm(X-W*(W'*X), 'fro'));
         %%fprintf('\n');
    end
    
    % save intermediate results
    if ( exist('outputdir','var') && ~isempty(outputdir) )
        if ( mod(iter,save_step) == 0 )
             %%fprintf('Saving intermediate results ...');
            save([outputdir 'IntermResultsExtractBases.mat'],'iter','D','K','W','-v7.3') ;
             %%fprintf('done\n') ;
        end
    end
end
% Reordering the output - putting them in standard form, loosely following
% what is done in the nnmf matlab function
H = W'*X ;
hlen = sqrt(sum(H.^2,2));
if any(hlen==0)
    warning(message('opnmf:LowRank', K - sum( hlen==0 ), K));
    hlen(hlen==0) = 1;
end
clear H
Wh = bsxfun(@times,W,hlen');
% Then order by W
[~,idx] = sort(sum(Wh.^2,1),'descend'); clear Wh
W = W(:,idx);
H = W'*X ;
% put results to original dimension
WW = zeros(Dinit,K) ;
WW(mean_im>0,:) = W ; clear W;
W = WW ; clear WW
end

function [W,H] = NNDSVD(A,k,flag)
%% NNDSVD    
% This function implements the NNDSVD algorithm described in [1] for
% initialization of Nonnegative Matrix Factorization Algorithms.
%
% [W,H] = nndsvd(A,k,flag);
%  
% INPUT
% ------------
%
% A    : the input nonnegative m x n matrix A
% k    : the rank of the computed factors W,H
% flag : indicates the variant of the NNDSVD Algorithm
%        flag = 0 --> NNDSVD
%        flag = 1 --> NNDSVDa
%        flag = 2 --> NNDSVDar
%        flag = 3 --> NNDSVD using random SVD calculation
%
% OUTPUT
% -------------
%   
% W   : nonnegative m x k matrix
% H   : nonnegative k x n matrix
%
% 
% References:
% 
% [1] C. Boutsidis and E. Gallopoulos, SVD-based initialization: A head
%     start for nonnegative matrix factorization, Pattern Recognition,
%     Elsevier
%----------------------check the input matrix------------------------------
if numel(find(A<0)) > 0
    error('The input matrix contains negative elements !')
end
%--------------------------------------------------------------------------
%size of the input matrix
[m,n] = size(A);
%the matrices of the factorization
W = zeros(m,k);
H = zeros(k,n);
% 1st SVD --> partial SVD rank-k to the input matrix A.
if ( flag == 3)
    % use random svd for efficient computation
    l = max(3*k,20) ;
    [U,S,V]=randpca(A,k,true,8,l);
else
    % use standard matlab svn implementation
    [U,S,V] = svds(A,k);
end
%choose the first singular triplet to be nonnegative
W(:,1)     =  sqrt(S(1,1)) * abs(U(:,1) );         
H(1,:)     =  sqrt(S(1,1)) * abs(V(:,1)'); 
% 2nd SVD for the other factors (see table 1 in our paper)
for i=2:k
    uu = U(:,i); vv = V(:,i);
    uup = pos(uu); uun = neg(uu) ;
    vvp = pos(vv); vvn = neg(vv);
    n_uup = norm(uup);
    n_vvp = norm(vvp) ;
    n_uun = norm(uun) ;
    n_vvn = norm(vvn) ;
    termp = n_uup*n_vvp; termn = n_uun*n_vvn;
    if (termp >= termn)
        W(:,i) = sqrt(S(i,i)*termp)*uup/n_uup; 
        H(i,:) = sqrt(S(i,i)*termp)*vvp'/n_vvp;
    else
        W(:,i) = sqrt(S(i,i)*termn)*uun/n_uun; 
        H(i,:) = sqrt(S(i,i)*termn)*vvn'/n_vvn;
    end
end
%------------------------------------------------------------
%actually these numbers are zeros
W(find(W<0.0000000001))=0.1;
H(find(H<0.0000000001))=0.1;
if(exist('flag','var'))
    % NNDSVDa: fill in the zero elements with the average
    if flag==1
        ind1      =  find(W==0) ;
        ind2      =  find(H==0) ;
        average   =  mean(A(:)) ;
        W( ind1 ) =  average    ;
        H( ind2 ) =  average    ;
        % NNDSVDar: fill in the zero elements with random values in the space [0:average/100]
    elseif flag==2
        ind1      =  find(W==0) ;
        ind2      =  find(H==0) ;
        n1        =  numel(ind1);
        n2        =  numel(ind2);        
        average   =  mean(A(:))       ;
        W( ind1 ) =  (average*rand(n1,1)./100)  ;
        H( ind2 ) =  (average*rand(n2,1)./100)  ;
    end
end
end

%%
function p = get_permutation_p(tvalue,rvalue,direction)

if nargin < 3 || isempty(direction)
    direction = 'two';
end

v = [tvalue;rvalue(:)];

if strcmp(direction,'left') || (strcmp(direction,'two') & tvalue <= mean(rvalue))
    [a,b] = sort(v); f = find(b==1);
    p = f/numel(v);
elseif strcmp(direction,'right') || (strcmp(direction,'two') & tvalue > mean(rvalue))
    [a,b] = sort(v,'descend'); f = find(b==1);
    p = f/numel(v);
end

if strcmp(direction,'two')
    p = p*2;
end

end

function [val, idx] = yael_kmax(v, k)
%% 
[val, idx] = sort (v, 'descend');
val = val (1:k, :);
idx = idx (1:k, :);
end

function [Ap] = pos(A)
%%
Ap = (A>=0).*A;
end

function [Am] = neg(A)
%%
Am = (A<0).*(-A);
end

function [dataset,group_info] = classdata(data,class,group,DM)
%%
if strcmp(class,'cat')
   g = unique(group);
   for i = 1:length(g)
       f = find(group==g(i));
       dataset{i} = data(f,:);
   end
   group_info = [];
elseif strcmp(class,'bin') 
   [B,E] = discretize(group,DM.binB);
   group_info.Group = B;
   group_info.Boundary = E;
   g = unique(B);
   for i = 1:length(g)
       f = find(B==g(i));
       dataset{i} = data(f,:);
   end
elseif strcmp(class,'window') 
   g = unique(group);
   g = sort(g);
   P = sliding_window(g,DM.windowL,DM.windowW);
   for i = 1:size(P,1)
       f = P(i,:);
       dataset{i} = data(f,:);
   end  
  group_info = [];
end

end

function P = sliding_window(A,L,W)
%%
% A: vector of data
% L: length of windows 
% W: step width 

if mod(length(A) - L,W) > 0
    N = fix((length(A) - L)/W) + 1;
else
    N = (length(A) - L)/W;
end
     
clear P; P(1,:) = A(1:L);
ep = L;

for i = 1:N   
    if ep+W <= length(A)
        c = P(i,:); c(1:W) = []; 
        c = [c,A([ep+1:ep+W])];
    else
        c = A([length(A)-L+1:length(A)]);
    end
    P(i+1,:) = c;
    ep = ep + W;
end    
end

%% ----------------------------------------------------------------------------------------------------------------------------------------
function [W,H] = factor_analysis(X,K,OS)
%%
if nargin < 3 || isempty(OS)
    OS.fa_scores = 'wls';
    OS.fa_start = 'Rsquared';
    OS.fa_rotate = 'varimax';
    OS.fa_delta = 0.005;
end

if ~isfield(OS,'fa_scores') || isempty(OS.fa_scores)
    OS.fa_scores = 'wls';
end
if ~isfield(OS,'fa_start') || isempty(OS.fa_start)
    OS.fa_start = 'Rsquared';
end
if ~isfield(OS,'fa_rotate') || isempty(OS.fa_rotate)
    OS.fa_rotate = 'varimax';
end
if ~isfield(OS,'fa_delta') || isempty(OS.fa_delta)
    OS.fa_delta = 0.005;
end

X = X';

[numSamples, numFeatures] = size(X);
if numSamples <= numFeatures
    error(['Factor analysis requires the number of samples to be ' ...
           'greater than the number of features. Possible reasons:\n' ...
           '1. Your data in one specific fold has too few samples relative to the number of features.\n' ...
           '2. Consider collecting more data or reducing the dimensionality ' ...
           'of your data using PCA or expanding the size of the validation set, such as using split-half rather than five-folds.']);
end

% computes the maximum likelihood estimate (MLE) of the factor loadings matrix ¦« in the factor analysis model
[W,psi,T,stats,F] = factoran(X,K,'Scores',OS.fa_scores,'Start',OS.fa_start,'Rotate',OS.fa_rotate,'maxit',1000,'Delta',OS.fa_delta) ;
H = [];
end

function [W,H] = principal_component_analysis(X,K,OS)
%%
if nargin < 3 || isempty(OS)
    OS.pca_algorithm = 'svd';
    OS.pca_centered = false;
    OS.pca_options = [];
    OS.pca_VariableWeights = [];
    OS.fa_rotate = []
end

if ~isfield(OS,'pca_algorithm') || isempty(OS.pca_algorithm)
    OS.pca_algorithm = 'svd';
end
if ~isfield(OS,'pca_centered') || isempty(OS.pca_centered)
    OS.pca_centered = false;
end
if ~isfield(OS,'pca_options')  
    OS.pca_options = [];
end
if ~isfield(OS,'pca_VariableWeights')  
    OS.pca_VariableWeights = [];
end
if ~isfield(OS,'fa_rotate')  
    OS.fa_rotate = 'varimax';
end

if strcmp(OS.pca_VariableWeights,'variance')
    W = pca(X','NumComponents',K,'Algorithm',OS.pca_algorithm,'Centered',OS.pca_centered,'VariableWeights','variance','Options',OS.pca_options);
else
    W = pca(X','NumComponents',K,'Algorithm',OS.pca_algorithm,'Centered',OS.pca_centered,'Options',OS.pca_options); 
end
H = [];

if ~isempty(OS.fa_rotate)
    W = rotatefactors(W, 'Method', OS.fa_rotate, 'Maxit', 1000, 'reltol', 1e-8);
end

end

function [W,H] = principal_axis_factoring(X, n_factors,OS)
%%
% Principal Axis Factoring (PAF) implementation in MATLAB
% X: input data matrix
% n_factors: number of factors to extract

if nargin < 3 || isempty(OS)
    OS.fa_rotate = 'varimax';
end
if ~isfield(OS,'fa_rotate') || isempty(OS.fa_rotate)
    OS.fa_rotate = 'varimax';
end

n_iterations = 100;

X = zscore(X');

R = corr(X); 

[~, D] = eig(R);
d = diag(D);
[d, ind] = sort(d, 'descend');
h2 = sum(d(1:n_factors))/size(R,1); 
h2 = repmat(h2,size(R,1),1);

for iter = 1:n_iterations
    R_reduced = R;
    for i = 1:size(R, 1)
        R_reduced(i, i) = h2(i);
    end
    
    [V, D] = eig(R_reduced);
    [d, ind] = sort(diag(D), 'descend');
    V = V(:, ind); 
    
    loadings = V(:, 1:n_factors) * diag(sqrt(d(1:n_factors)));
    h2_new = sum(loadings.^2, 2);
    
    if max(abs(h2_new - h2)) < 1e-6
        break
    end
    h2 = h2_new;
end

loadings = V(:, 1:n_factors) * diag(sqrt(d(1:n_factors)));

try
    loadings = rotatefactors(loadings, 'Method', OS.fa_rotate, 'Maxit', 1000, 'reltol', 1e-8);
catch ME
    warning('Varimax rotation failed to converge. Returning unrotated loadings.');
    loadings = V(:, 1:n_factors) * diag(sqrt(d(1:n_factors)));
end

% eigenvalues = sum(loadings.^2, 1); 
% unrotated_eigenvalues = d(1:n_factors);

W = loadings; H = []; 
end