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
    p.addParameter('use_parallel',false,@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.addParameter('maxIter',25,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));
    p.addParameter('minResponsiveFrac',0,@(x)validateattributes(x,{'numeric'},{'scalar','positive','<',1}));
    p.addParameter('minSpkParamRatio',0,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    p.addParameter('separate_clicks_by_side',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('distribution','poisson',@(x)validatestring(x,{'poisson','normal'}));
    p.addParameter('save_path','');
    p.addParameter('useGPU',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('bin_size_ms',4);
    p.parse(varargin{:});
    params=p.Results;
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
        fprintf('Testing save now before spending a long time fitting... ');
        if isempty(params.save_path)
            [a,b,c] = fileparts(Cells.mat_file_name);
        else
            [a,b,c] = fileparts(params.save_path);            
        end
        mat_file_name = fullfile(a,[b,'_glmfits_save_test.mat']);
        test=[];
        save(mat_file_name,'test','-v7.3');
        delete(mat_file_name);   
        fprintf(' Success! Saved and deleted %s.\n',mat_file_name);
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
    dspec = buildGLM.initDesignSpec(expt);    
    dspec = build_dspec_for_pbups(dspec,covariates,[]);
    %% loop over cells
    for c=1:length(params.cellno)
        S=struct();
        S.cellno = params.cellno(c);          
        fprintf('Cell id %g (%g of %g to fit):\n',S.cellno,c,rawData.param.ncells);
        %% determine if cell is responsive enough to fit
        S.responsiveFrac = sum(arrayfun(@(x)numel(x.(['sptrain',num2str(S.cellno)])),expt.trial)>0)./nTrials;        
        if params.minResponsiveFrac>0
            if S.responsiveFrac < params.minResponsiveFrac
                fprintf('Cell %g only fired a spike on %0.2g%% of trials. Moving on without fitting.\n',S.cellno,S.responsiveFrac*100);
                continue
            end
        end
        %% add spike history to dspec for current cell
        S.dspec = build_dspec_for_pbups(dspec,'spike_history',S.cellno); 
        dm = buildGLM.compileSparseDesignMatrix(S.dspec, 1:nTrials);  
        dm = buildGLM.removeConstantCols(dm);
        Y = full(buildGLM.getBinnedSpikeTrain(expt, ['sptrain',num2str(S.cellno)], dm.trialIndices)); 
        S.Y=Y;
        S.dm=dm;
        %% determine if spike/parameter ratio is acceptable
        S.totalSpks = sum(Y);
        S.spkParamRatio = S.totalSpks ./ (size(dm.X,2)+1);   
        fields =fieldnames(S);
        for f=1:length(fields)
            stats(c).(fields{f}) = S.(fields{f});
        end        
        if params.minSpkParamRatio>0
            if S.spkParamRatio < params.minSpkParamRatio
                fprintf('Cell %g only has %g spikes to %g params to be fit. Moving on without fitting.\n',S.cellno,S.totalSpks,size(dm.X,2)+1);
                continue
            end
        end    
        if rank(dm.X) < size(dm.X,2)
            warning('Design matrix is not full rank for cell %g. Skipping cell.',S.cellno);
            continue;
        end        
        %% Fitting UN cross-validated model
        tic;fprintf('   Fitting UN cross-validated model ... ');   
        options = statset('MaxIter',params.maxIter);        
        if params.useGPU
            dm.X = gpuArray(dm.X);
            Y = gpuArray(Y);
            bias_column_fun = @(x)[gpuArray.ones(size(x,1),1),x];
        else
            bias_column_fun = @(x)[ones(size(x,1),1),x];
        end
%         %%
%         Xtrain = bias_column_fun(dm.X);
%         spstrain = Y;
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
        options = statset('MaxIter',params.maxIter);                
        tic;S.init_beta = gather(regress(bias_column_fun(dm.X),Y));toc
%         link.Inverse = @(x)max(x/50,x);
%         link.Link = @(x)min(x,50*x);
%         link.Derivative = @(x)delta(x);        
%          link.Inverse = @(x)log(1+exp(x));
%          link.Link = @(x)log(exp(x)-1);
%          link.Derivative = @(x) ( 1./(1-exp(-x)));
% %         profile off;tic;mdl = fitglm(dm.X,Y,'Intercept',true,'Distribution','poisson','DispersionFlag',false,'Options',statset('UseParallel',true));        toc;
%         link.Link=@(x)x;
%         link.Inverse = @(x)x;
%         link.Derivative = @(x)1;
%         link.Link=@(x)log(x);
%         link.Inverse=@(x)exp(x);
%         link.Derivative = @(x)(1./x);
        params.distribution='poisson';
        tic;[~, S.dev, stat_temp] = glmfit(dm.X, Y, params.distribution,'options',options);toc
        fields_to_copy = fieldnames(stat_temp);
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
        S.Yhat=glmval(S.beta,gather(dm.X),link);
        S.Yhat_init_beta = gather([ones(size(dm.X,1),1),dm.X]*S.init_beta);
        fprintf('took %s.\n',timestr(toc));
        % reconstruct fitted kernels by weighted combination of basis functions
        [S.ws,S.wvars] = buildGLM.combineWeights(buildGLM.addBiasColumn(dm), S.beta,S.covb );
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
                combineWeightFun = @(raw_weights,covariances)buildGLM.combineWeights(buildGLM.addBiasColumn(dm), raw_weights , covariances);
                getSpkIdxFun = @(trial_idx)buildGLM.getSpikeIndicesforTrial(expt,trial_idx);        
                fprintf('   Fitting under %g-fold cross-validation ... ',params.kfold);    
                tic;
                for i=params.kfold:-1:1
                    train_idx{i} = getSpkIdxFun(S.cvp.training(i));
                    test_idx{i} = getSpkIdxFun(S.cvp.test(i));
                end
                clear dev cv_stats init_beta
                if params.use_parallel
                    X=dm.X;
                    parfor i=1:params.kfold                       
                        [~, dev(i), cv_stats(i)] = glmfit(X(train_idx{i},:),Y(train_idx{i}), 'poisson','options',options);  
                        init_beta{i} = gather(regress(bias_column_fun(X(train_idx{i},:)),Y(train_idx{i})));                                                     
                    end    
                else
                    for i=params.kfold:-1:1
                        [~, dev(i), cv_stats(i)] = glmfit(dm.X(train_idx{i},:),Y(train_idx{i}), 'poisson','options',options);  
                        init_beta{i} = gather(regress(bias_column_fun(dm.X(train_idx{i},:)),Y(train_idx{i})));                                                     
                    end                  
                end
                cv_stats = rmfield(cv_stats,{'resid','residp','residd','resida'}); % these take up A TON of space, and could always be generated if needed
                for i=1:params.kfold
                    cv_stats(i).dev = dev(i);
                    fields = fieldnames(cv_stats);
                    for f=1:length(fields)
                        cv_stats(i).(fields{f}) = gather(cv_stats(i).(fields{f}));
                    end                    
                    cv_stats(i).Yhat=glmval(cv_stats(i).beta,gather(dm.X(test_idx{i},:)),'log',cv_stats(i));
                    cv_stats(i).init_beta = init_beta{i};
                    cv_stats(i).Yhat_init_beta = gather([ones(size(dm.X(test_idx{i},:),1),1),dm.X(test_idx{i},:)]*cv_stats(i).init_beta);                    
                    [cv_stats(i).ws,cv_stats(i).wvars] =combineWeightFun(cv_stats(i).beta,cv_stats(i).covb);
                end
                S.cv_stats=cv_stats;
                fprintf('Took %s.\n',timestr(toc));            
            end
        end
        S.params = params;        
        %S.covariate_stats = get_covariate_stats(S);        
        fields =fieldnames(S);
        for f=1:length(fields)
            stats(c).(fields{f}) = S.(fields{f});
        end
    end
    if params.save && isfield(stats,'dev')
        if isempty(params.save_path)
            mat_file_name = strrep(mat_file_name,'glmfits_save_test','glmfits');
            save(mat_file_name,'stats','-v7.3');
            fprintf('Saved fit stats successfully to %s.\n',mat_file_name);
        else
            mat_file_name = params.save_path;
            save(mat_file_name,'stats','-v7.3');
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