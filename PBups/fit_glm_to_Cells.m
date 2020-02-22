function stats = fit_glm_to_Cells(Cells,varargin)
    % fits a GLM for spiking data recorded during PBups and contained within
    % an Cells data structure.
    % This function is essentially a wrapper for the neuroGLM package forked from the
    % Pillow lab.
    %% parse and validate inputs
    p=inputParser;
    p.KeepUnmatched=true;    
    p.addParameter('cellno',[]);
    p.addParameter('kfold',[]);
    p.addParameter('save',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('use_parallel',false,@(x)validateattributes(x,{'logical'},{'scalar'}));    % parfor operates over cells by default. But alter the code to operate over cross-validation folds (as indicated with comments in the code) if you are using cross-validation (i.e. kfold>1)
    p.addParameter('maxIter',25,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));
    p.addParameter('minResponsiveFrac',0.5,@(x)validateattributes(x,{'numeric'},{'scalar','positive','<',1})); % fraction of trials on which the cell fired at least one spike
    p.addParameter('minSpkParamRatio',10,@(x)validateattributes(x,{'numeric'},{'scalar','positive'})); % minimum number of spikes per model parameter (i.e. design matrix column)
    p.addParameter('separate_clicks_by_side',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('distribution','poisson',@(x)validatestring(x,{'poisson','normal'}));
    p.addParameter('save_path','');
    p.addParameter('useGPU',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('bin_size_ms',1);
    p.parse(varargin{:});
    params=p.Results;
    %% load Cells file if a filename is passed
    if ischar(Cells) && exist(Cells,'file')
        mat_file_name = Cells;
        [a,b,c] = fileparts(mat_file_name);
        fprintf('Loading %s ... ',[b,c]);tic;        
        Cells = load(Cells);       
        fprintf('took %s.\n',timestr(toc));
        Cells.mat_file_name = mat_file_name;
    end
    %% make rawData and expt structures (i.e. put event times and spike times into neuroGLM format)
    rawData = make_glm_trials_from_Cells(Cells,varargin{:});   
    expt=build_expt_for_pbups(rawData,params.bin_size_ms); 
    nTrials = rawData.nTrials;
    if isempty(params.cellno)
        params.cellno=1:rawData.param.ncells;
        if isempty(params.cellno)
            warning('No cells in Cells!');
            return
        end
    end
    fprintf('\n');
    %% testing save
    if params.save
        fprintf('Testing save now before spending a long time fitting... \n');
        if isempty(params.save_path)
            [a,b,c] = fileparts(Cells.mat_file_name);
        else
            [a,b,c] = fileparts(params.save_path);            
        end
        mat_file_name = fullfile(a,[b,'_glmfits_save_test.mat']);
        test=[];
        save(mat_file_name,'test','-v7');
        delete(mat_file_name);   
        [a,b,c] = fileparts(mat_file_name);
        fprintf(' Success! Saved and deleted %s\nin directory %s.\n',[b,c],a);
    end    
    %% initialize dspec elements that are common to all cells
    % dspec states what the regressors are in your model and how they are parameterized    
    covariates = rawData.timings;
    if ismember('spoke',covariates)
        covariates = union(covariates,{'spoke_left_hit','spoke_right_hit','spoke_left_miss','spoke_right_miss'});
        covariates = setdiff(covariates,{'spoke'});
    end
    if ismember('cpoke_out',covariates)
        covariates = union(covariates,{'cpoke_out_left','cpoke_out_right'});
        covariates = setdiff(covariates,{'cpoke_out'});
    end    
    if params.separate_clicks_by_side
        covariates(strmatch('all_click',covariates))='';
    else
        covariates(strmatch('left_click',covariates))='';
        covariates(strmatch('right_click',covariates))='';        
    end
    dspec_base = buildGLM.initDesignSpec(expt);    
    dspec_base = build_dspec_for_pbups(dspec_base,covariates,[]);
    ncells=rawData.param.ncells;
    %% loop over cells
    % do a non-parallel loop to compute things needed for each cell without
    % passing large necessary structures to workers
    responsive_enough=true(1,ncells);
    if params.useGPU
        bias_column_fun = @(x)[gpuArray.ones(size(x,1),1),x];
    else
        bias_column_fun = @(x)[ones(size(x,1),1),x];
    end     
    n_params=sum([dspec_base.covar.edim])+1+6; % 6 is for spike history filter but this may be an overestimate. roughly good enough for these purposes.
    % non-parallel loop across cells just to compute things that otherwise
    % would require passing the entire list of spike times to all parallel
    % workers
    [X,Y,D] = deal(cell(1,length(params.cellno)));
    fprintf('\nBuilding design matrices for responsive cells ...\n');
    for c=1:length(params.cellno)
        responsiveFrac = sum(arrayfun(@(x)numel(x.(['sptrain',num2str(params.cellno(c))])),expt.trial)>0)./nTrials;         
        totalSpikes = length(cat(1,expt.trial.(['sptrain',num2str(params.cellno(c))])));
        %% determine if cell is responsive enough to fit
        if params.minResponsiveFrac>0
            if responsiveFrac < params.minResponsiveFrac || totalSpikes<n_params*params.minSpkParamRatio        
                responsive_enough(c)=false;
                continue
            end
        end        
        fprintf('%g, ',params.cellno(c));        
        dspec = build_dspec_for_pbups(dspec_base,'spike_history',params.cellno(c)); 
        % build design matrix
        dm = buildGLM.compileSparseDesignMatrix(dspec, 1:nTrials);  
        dm = buildGLM.removeConstantCols(dm);       
        Y{c} = full(buildGLM.getBinnedSpikeTrain(expt, ['sptrain',num2str(params.cellno(c))], dm.trialIndices)); 
        if params.useGPU
            X{c} = gpuArray(dm.X);
            Y{c} = gpuArray(Y{c});
        else
            X{c} = dm.X;
        end
        dm = rmfield(dm,'X');
        dm = rmfield(dm,'dspec');
        D{c}=dm;
    end
    fprintf('\n%g of %g cells are responsive enough to be fit.\n',sum(responsive_enough),ncells);
    responsive_cells = find(responsive_enough);
    X = X(responsive_enough);
    Y = Y(responsive_enough);
    %% trial contains the list of ALL spike times, which is no longer needed    
    trial_fields= fieldnames(dspec.expt.trial);
    is_spk_field = strncmp(trial_fields,'sp',2);
    dspec.expt.trial = rmfield(dspec.expt.trial,trial_fields(is_spk_field)); 
    %% loop over cells
    parfor c=1:sum(responsive_enough) % PARFOR ON THIS LINE WILL FIT CELLS IN PARALLEL
        S=struct();
        S.dm=D{c};
        S.cellno = params.cellno(responsive_cells(c));        
        fprintf('Cell id %g (%g of %g to fit):\n',S.cellno,c,sum(responsive_enough));        
        fields =fieldnames(S);
        for f=1:length(fields)
            stats(c).(fields{f}) = S.(fields{f});
        end          
        %% Fitting UN cross-validated model
        tic;fprintf('   Fitting UN cross-validated model ... ');   
        options = statset('MaxIter',params.maxIter);        
        %% below is code in dev for pass-glm fitting
%         %%
%         Xtrain = bias_column_fun(dm.X);
%         spstrain = Y{c};
%         dtSp=0.001;
%         nparams = size(Xtrain,2);
%         Cinv = 1*eye(nparams); Cinv(1,1)=0;        
%         xlim = [0,3]; % interval of polynomial approximation
%         dx = 0.01; % parameter for computing polynomial coefficients. generally found this to be sufficient in all applications.
%         what_cheby = compute_chebyshev(@(x) exp(x),xlim,dx); % computing the polynomial coefficients
%         a = what_cheby(1); b = what_cheby(2); c = what_cheby(3); % get the coefficients
%         w_hat_exp = (2.0*c*dtSp*Xtrain'*Xtrain + Cinv)\(Xtrain'*spstrain-b*dtSp*sum(Xtrain,1)'); % compute weights using approximation        
%         %%
%         options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'MaxIterations',1000);
%         tic;
%         w = fminunc(@(w) negloglik_grad(w,spstrain,Xtrain,dtSp),zeros(size(Xtrain,2),1),options);        toc;
        %%
        S.init_beta = gather(regress(bias_column_fun(X{c}),Y{c}));
        %% below is code in dev for using non-canonical link functions       
%         link.Inverse = @(x)max(x/50,x);
%         link.Link = @(x)min(x,50*x);
%         link.Derivative = @(x)delta(x);        
%          link.Inverse = @(x)log(1+exp(x));
%          link.Link = @(x)log(exp(x)-1);
%          link.Derivative = @(x) ( 1./(1-exp(-x)));
% %         profile off;tic;mdl = fitglm(X{c},Y{c},'Intercept',true,'Distribution','poisson','DispersionFlag',false,'Options',statset('UseParallel',true));        toc;
%         link.Link=@(x)x;
%         link.Inverse = @(x)x;
%         link.Derivative = @(x)1;
%         link.Link=@(x)log(x);
%         link.Inverse=@(x)exp(x);
%         link.Derivative = @(x)(1./x);
        [~, S.dev, stat_temp] = glmfit(X{c}, Y{c}, params.distribution,'options',options);
        fields_to_copy = fieldnames(stat_temp);
        fields_to_copy = fields_to_copy(~ismember(fields_to_copy,{'residp','residd','resida'}));
        nf=length(fields_to_copy);
        for f=1:nf
            S.(fields_to_copy{f}) = gather(stat_temp.(fields_to_copy{f}));
        end     
        S.dev=gather(S.dev);
        % compute model predicted firing rates
        switch params.distribution
            case 'normal'
                link = 'identity';
            case 'poisson'
                link = 'log';
        end
        S.Yhat=glmval(S.beta,gather(X{c}),link);
        S.Yhat_init_beta = gather([ones(size(X{c},1),1),X{c}]*S.init_beta);
        fprintf('took %s.\n',timestr(toc));
        % reconstruct fitted kernels by weighted combination of basis functions
        S.dm.biasCol=1;
        [S.ws,S.wvars] = buildGLM.combineWeights(S.dm, dspec, S.beta,S.covb,false );
        % determine if least-squared weights are badly scaled. If so, not much point
        % doing cross-validation.
        if any(sqrt(S.wts)~=0 & sqrt(S.wts)<(max(sqrt(S.wts))*eps('double')^(2/3)))
            S.badly_scaled=true;
        else
            S.badly_scaled=false;
        end
        % Fit cross-validated model (if requested and if uncross-validated fit was not badly scaled)
        if ~isempty(params.kfold)
            if S.badly_scaled
                fprintf('Skipping cross-validation since fit to all data was badly scaled.\n');
            else
                S.cvp = cvpartition(nTrials,'KFold',params.kfold);
                combineWeightFun = @(raw_weights,covariances)buildGLM.combineWeights(S.dm,dspec, raw_weights , covariances,false);
                getSpkIdxFun = @(trial_idx)buildGLM.getSpikeIndicesforTrial(dspec.expt,trial_idx);        
                fprintf('   Fitting under %g-fold cross-validation ... ',params.kfold);    
                tic;
                [train_idx,test_idx,init_beta] = deal(cell(1,params.kfold));                
                for i=1:params.kfold
                    train_idx{i} = getSpkIdxFun(S.cvp.training(i));
                    test_idx{i} = getSpkIdxFun(S.cvp.test(i));
                end
                [cv_stats,dev] = zeros(1,params.kfold);
                if params.use_parallel
                    x=X{c};
                    for i=1:params.kfold % PARFOR ON THIS LINE WILL FIT CROSS-VALIDATION FOLDS IN PARALLEL
                        [~, dev(i), cv_stats(i)] = glmfit(x(train_idx{i},:),Y{c}(train_idx{i}), 'poisson','options',options);  
                        init_beta{i} = gather(regress(bias_column_fun(x(train_idx{i},:)),Y{c}(train_idx{i})));                                                     
                    end    
                else
                    for i=params.kfold:-1:1
                        [~, dev(i), cv_stats(i)] = glmfit(X{c}(train_idx{i},:),Y{c}(train_idx{i}), 'poisson','options',options);  
                        init_beta{i} = gather(regress(bias_column_fun(X{c}(train_idx{i},:)),Y{c}(train_idx{i})));                                                     
                    end                  
                end
                cv_stats = rmfield(cv_stats,{'resid','residp','residd','resida'}); % these take up A TON of space, and could always be generated if needed
                for i=1:params.kfold
                    cv_stats(i).dev = dev(i);
                    fields = fieldnames(cv_stats);
                    for f=1:length(fields)
                        cv_stats(i).(fields{f}) = gather(cv_stats(i).(fields{f}));
                    end                    
                    cv_stats(i).Yhat=glmval(cv_stats(i).beta,gather(X{c}(test_idx{i},:)),'log',cv_stats(i));
                    cv_stats(i).init_beta = init_beta{i};
                    cv_stats(i).Yhat_init_beta = gather([ones(size(X{c}(test_idx{i},:),1),1),X{c}(test_idx{i},:)]*cv_stats(i).init_beta);                    
                    [cv_stats(i).ws,cv_stats(i).wvars] =combineWeightFun(cv_stats(i).beta,cv_stats(i).covb);
                end
                S.cv_stats=cv_stats;
                fprintf('Took %s.\n',timestr(toc));            
            end
        end
        %S.covariate_stats = get_covariate_stats(S);        
        fields =fieldnames(S);
        for f=1:length(fields)
            stats(c).(fields{f}) = S.(fields{f});
        end
    end
    params.dspec = dspec;
    if params.save && isfield(stats,'dev')
        if isempty(params.save_path)
            mat_file_name = strrep(mat_file_name,'glmfits_save_test','glmfits');
            save(mat_file_name,'stats','params','-v7');
            fprintf('Saved fit stats successfully to %s.\n',mat_file_name);
        else
            mat_file_name = params.save_path;
            save(mat_file_name,'stats','params','-v7');
            fprintf('Saved fit stats successfully to %s.\n',mat_file_name);            
        end
    end
end

function beta = regress(x,y)
    [Q,R] = qr(x,0);
    beta=R\(Q'*y);
end

function y = delta(x)
    if x<0
        y=50;
    else
        y=1;
    end
end