function stats = fit_glm_to_Cells(Cells,varargin)
    % fits a GLM for spiking data recorded during PBups and contained within
    % an Cells data structure.
    % This function is essentially a wrapper for the neuroGLM package forked from the
    % Pillow lab.
    %% parse and validate inputs
    p=inputParser;
    p.KeepUnmatched=true;    
    p.addParameter('cellno',[]);
    p.addParameter('kfold',1,@(x)validateattributes(x,{'numeric'},{'scalar','integer','>',0}));
    p.addParameter('save',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('create_pool',false,@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.addParameter('n_workers',1,@(x)validateattributes(x,{'numeric'},{'scalar','integer','>',0}));    % parallelization operates over cells unless cross-validation is used (i.e. kfold>1) in which case it operates over cross-validation folds
    p.addParameter('maxIter',25,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));
    p.addParameter('bin_size_s',0.001,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));  % resolution of the model. predictions have this bin size.  
    p.addParameter('minResponsiveFrac',0.5,@(x)validateattributes(x,{'numeric'},{'scalar','positive','<',1}));
    p.addParameter('minSpkParamRatio',10,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    p.addParameter('distribution','poisson',@(x)validateattributes(x,{'char'},{}));   
    p.addParameter('link','canonical',@(x)validatestring(x,{'log','identity','softplus'}));
    p.addParameter('save_path','',@(x)validateattributes(x,{'char'},{})); % directory where files will be saved
    p.addParameter('useGPU',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('phi',1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'})); % value of C (adaptation state) at time lag 0
    p.addParameter('tau_phi',0,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'})); % recovery time constant of C
    p.addParameter('within_stream',false,@(x)validateattributes(x,{'logical'},{'scalar'})); % within_stream adaptation flag. false by default, given empirical findings in Brunton et al. 2015    
    p.addParameter('fit_adaptation',true,@(x)validateattributes(x,{'logical'},{'scalar'})); % whether or not to fit phi and tau_phi for each neuron using gradient descent
    p.parse(varargin{:});
    params=p.Results;
    validatestring(params.distribution,{'poisson','normal'},mfilename,'distribution');
    switch params.link
        % deal with special cases
        case 'canonical'
            if strcmp(params.distribution,'poisson')
                params.link='log';
            else
                params.link='identity';                
            end
        case 'softplus'
            params.link.Link = @(mu)log(1+exp(mu));
            params.link.Derivative = @(mu)1./(1+exp(-mu));
            params.link.Inverse = @(eta)log(exp(eta)-1);
    end        
    if ischar(params.link)
        [link.Link,link.Derivative,link.Inverse] = stattestlink(params.link,'double');        
        params.link=link;
    end
    %% load Cells file if a filename is passed
    if ischar(Cells) && exist(Cells,'file')
        mat_file_name = Cells;
        [~,b,c] = fileparts(mat_file_name);
        fprintf('Loading %s ... ',[b,c]);tic;        
        Cells = load(Cells);       
        fprintf('took %s.\n',timestr(toc));
        Cells.mat_file_name = mat_file_name;
    end
    %% make rawData and expt structures (i.e. put event times and spike times into neuroGLM format)
    rawData = make_glm_trials_from_Cells(Cells,varargin{:}); 
    expt=build_expt_for_pbups(rawData,'bin_size_s',params.bin_size_s); 
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
            [a,b,~] = fileparts(Cells.mat_file_name);
            save_subdir = [datestr(now,'YYYY_mm_dd_HH_MM_SS'),'_glmfits'];
            params.save_path = fullfile(a,b,save_subdir);
        end
        mat_file_name = fullfile(params.full_save_path,'_glmfits_save_test.mat');
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
    dspec_base = buildGLM.initDesignSpec(expt);    
    dspec_base = build_dspec_for_pbups(dspec_base,covariates,[]);
    dm = buildGLM.compileSparseDesignMatrix(dspec_base, 1:nTrials, params);
    X_base = dm.X;
    ncells=rawData.param.ncells;
    %% loop over cells
    % non-parallel loop to do some preparatory computation that will
    % reduce I/O to the workers for the main fitting stage
    responsive_enough=true(1,ncells);  
    n_params=sum([dspec_base.covar.edim])+1+6; % 6 is for spike history filter but this may be an overestimate. roughly good enough for these purposes.
    %% first (non-parallel) loop across cells
    % non-parallel loop across cells just to compute the design matrices
    % (to avoid passing a structure with all the spike times to the
    % workers)
    % this is a bit of a hack -- you could just extract the spike times
    % for each cell and pass it as a cell array. but this would require
    % significant rewriting of code.
    X = cell(1,length(params.cellno));
    tic;fprintf('\nBuilding design matrices for responsive cells ...');
    spikes=struct();    
    for c=length(params.cellno):-1:1
        responsiveFrac(c) = sum(arrayfun(@(x)numel(x.(['sptrain',num2str(params.cellno(c))])),expt.trial)>0)./nTrials;         
        totalSpikes(c) = length(cat(1,expt.trial.(['sptrain',num2str(params.cellno(c))])));
        %% determine if cell is responsive enough to fit
        if params.minResponsiveFrac>0
            if responsiveFrac(c) < params.minResponsiveFrac || totalSpikes(c)<n_params*params.minSpkParamRatio        
                responsive_enough(c)=false;
                continue
            end
        end        
        dspec = build_dspec_for_pbups(dspec_base,'spike_history',params.cellno(c)); 
        %% build design matrices        
        dm = buildGLM.compileSparseDesignMatrix(dspec, 1:nTrials, params , numel(dspec.covar));   % only remake the columns associated with the spike history term since the rest is common across cells
        dm.X(:,1:(dspec.edim-dspec.covar(end).edim))=X_base;
        dm = buildGLM.removeConstantCols(dm);       
        spikes(c).Y = full(buildGLM.getBinnedSpikeTrain(expt, ['sptrain',num2str(params.cellno(c))], dm.trialIndices)); 
        spikes(c).cellno=params.cellno(c);
        X{c} = dm.X;
    end
    if ~any(responsive_enough)
        stats=struct();
        fprintf(' took %s.\n%g of %g cells are sufficiently responsive. Returning.\n',timestr(toc),sum(responsive_enough),ncells);        
    else
        fprintf(' took %s.\n%g of %g cells are sufficiently responsive. Fitting now:\n',timestr(toc),sum(responsive_enough),ncells);
    end
    responsive_cells = find(responsive_enough);
    X = X(responsive_enough);
    spikes = spikes(responsive_enough);
    %% trial contains the list of ALL spike times, which is no longer needed  
    dm = rmfield(dm,'X');    
    trial_fields= fieldnames(dm.dspec.expt.trial);
    is_spk_field = strncmp(trial_fields,'sptrain',7);
    params.dm=dm;        
    params.dm.dspec.expt.trial = rmfield(dm.dspec.expt.trial,trial_fields(is_spk_field)); 
    %% loop over cells
    if params.n_workers>1 
        if params.create_pool
            delete(gcp('nocreate'));
            parpool(params.n_workers);
        end
        parfor c=1:sum(responsive_enough)
            cellno = params.cellno(responsive_cells(c));
            fprintf('Cell id %g (%g of %g to fit):\n',cellno,c,sum(responsive_enough));           
            [stats(c),spikes(c).Yhat,spikes(c).Yhat_cv] = mainLoop(X{c},spikes(c).Y,params);
            stats(c).cellno=cellno;               
        end    
    else
        for c=sum(responsive_enough):-1:1
            cellno = params.cellno(responsive_cells(c));            
            fprintf('Cell id %g (%g of %g to fit):\n',cellno,c,sum(responsive_enough));           
            [stats(c),spikes(c).Yhat,spikes(c).Yhat_cv] = mainLoop(X{c},spikes(c).Y,params);
            stats(c).cellno=cellno;   
        end
    end  
    %% save
    params.rat = Cells.rat;
    params.sess_date = Cells.sess_date;
    params.sessid = Cells.sessid;
    params.responsive_enough=responsive_enough;
    params.responsiveFrac=responsiveFrac;
    params.totalSpikes=totalSpikes;
    if params.save
        save(fullfile(params.save_path,'glmfit_stats.mat'),'params','stats','-v7');
        fprintf('Saved %s successfully.\n',fullfile(params.save_path,'glmfit_stats.mat'));
        save(fullfile(params.save_path,'glmfit_spikes.mat'),'params','spikes','-v7');
        fprintf('Saved %s successfully.\n',fullfile(params.save_path,'glmfit_spikes.mat'));        
    end
end

function [stats,Yhat,Yhat_cv] = mainLoop(X,Y,params)   
    % all the main fitting function needs is the binned observations (Y),
    % the design matrix (X, no bias column), the dm structure (which is the same
    % for all cells) and the parent function params
    % if z-scoring was performed on the design matrix, you'd need to have a
    % separate dm structure for each cell storing this
    tic;fprintf('   Fitting UN cross-validated model ... ');drawnow; 
    nTrials = numel(params.dm.trialIndices);
    Yhat_cv=[];
    [stats,Yhat] = fit(X, Y, params, 1:nTrials);   
    % compute model predicted firing rates
    fprintf('took %s.\n',timestr(toc));
    % reconstruct fitted kernels by weighted combination of basis functions
    params.dm.biasCol=1;
    [stats.ws,stats.wvars] = buildGLM.combineWeights(params.dm, params.dm.dspec, stats.beta,stats.covb,true );
    % determine if least-squared weights are badly scaled. If so, not much point
    % doing cross-validation.
    if any(sqrt(stats.wts)~=0 & sqrt(stats.wts)<(max(sqrt(stats.wts))*eps('double')^(2/3)))
        stats.badly_scaled=true;
    else
        stats.badly_scaled=false;
    end
    stats=rmfield(stats,'wts');
    % Fit cross-validated model (if requested and if uncross-validated fit was not badly scaled)
    if ~isempty(params.kfold) && params.kfold>1
        if stats.badly_scaled
            fprintf('Skipping cross-validation since fit to all data was badly scaled.\n');
        else
            stats.cvp = cvpartition(nTrials,'KFold',params.kfold);
            combineWeightFun = @(raw_weights,covariances)buildGLM.combineWeights(params.dm,params.dm.dspec, raw_weights , covariances,false);
            getSpkIdxFun = @(trial_idx)buildGLM.getSpikeIndicesforTrial(params.dm.dspec.expt,trial_idx);        
            tic;fprintf('   Fitting under %g-fold cross-validation ... ',params.kfold);    
            [train_idx,test_idx] = deal(cell(1,params.kfold));                
            for i=1:params.kfold
                train_idx{i} = getSpkIdxFun(stats.cvp.training(i));
                test_idx{i} = getSpkIdxFun(stats.cvp.test(i));
            end
            if params.n_workers>1 && isempty(getCurrentTask()) %% only parallelize here if not already in a parallel loop over cells
                if params.create_pool
                    delete(gcp('nocreate'));    
                    parpool(params.n_workers);
                end
                for i=params.kfold:-1:1
                    Xs{i}=X(train_idx{i},:);
                    idx{i} = find(stats.cvp.training(i));
                    Ys{i}=Y(train_idx{i});
                end
                parfor i=1:params.kfold 
                    cv_stats(i) = fit(Xs{i},Ys{i}, params, idx{i});  
                end    
            else
                for i=params.kfold:-1:1
                    cv_stats(i) = fit(X(train_idx{i},:),Y(train_idx{i}), params, find(stats.cvp.training(i)));  
                end                  
            end
            Yhat_cv=zeros(size(Yhat));
            for i=1:params.kfold
                [cv_stats(i).ws,cv_stats(i).wvars] =combineWeightFun(cv_stats(i).beta,cv_stats(i).covb);
                covar_idx=find(ismember({params.dm.dspec.covar.label},{'left_clicks','right_clicks'}));                         
                X_test=buildGLM.updateSparseDesignMatrix_covar(params.dm.dspec, find(stats.cvp.test(i)), cv_stats(i).adaptation_stats, covar_idx, X(test_idx{i},:));
                Yhat_cv(test_idx{i})=params.link.Inverse(cv_stats(i).beta(1)+X_test*cv_stats(i).beta(2:end));                
            end
            cv_stats = rmfield(cv_stats,'wts');            
            stats.cv_stats=cv_stats;
            fprintf('Took %s.\n',timestr(toc));            
        end
    end
    stats.covariate_stats = get_covariate_stats(stats,params);
    fields = fieldnames(stats.wvars);
    for f=1:length(fields) % you can remove covariance structure now that summary has been computed
        stats.wvars.(fields{f}) = rmfield(stats.wvars.(fields{f}),'cov');
    end
end

function [stats,Yhat] = fit(X,Y,params,trials)
    if params.fit_adaptation
        [transform,itransform] = deal(@(x)x);        
        %transform=@(x)log(1+x); % identity near zero, but logarithmic as x->Inf
        %itransform=@(x)max(0,exp(x)-1);
        % makes fits converge faster when tau_phi is large
        % I tried tanh but this needs to be scaled (by ~2e4) to allow large tau_phi
        % and this leads to slow convergence when tau_phi is small.
        params.transform=transform;
        params.itransform=itransform;
        if params.useGPU
            Y=gpuArray(Y);
        end   
        % init
        params.phi=0.5;
        params.tau_phi=0.3;
        covar_idx=find(ismember({params.dm.dspec.covar.label},{'left_clicks','right_clicks'}));        
        X_update_fun = @(phi,tau_phi)buildGLM.updateSparseDesignMatrix_covar(params.dm.dspec, trials, struct('phi',phi,'tau_phi',tau_phi,'within_stream',params.within_stream), covar_idx, X);
        time_at_start=tic;              
        osf = @(x,optimValues,state)optim_status_fun(x,optimValues,state,time_at_start,itransform);
        options=optimoptions('fmincon','UseParallel',false,'OutputFcn',osf,'Algorithm','interior-point','FunctionTolerance',eps,'StepTolerance',1e-8);  % stop if you are taking tiny steps but not if the change in LL is small -- sometimes the gradient is really small far from the optimum and you should keep going until you get near the basin              
        optim_fun = @(x)NLL_fun(X_update_fun,Y,itransform(x(1)),itransform(x(2)),rmfield(params,'dm'));                        
        [adaptation_stats.beta,adaptation_stats.NLL,adaptation_stats.exitflag,adaptation_stats.output,adaptation_stats.lambda,...
            adaptation_stats.grad,adaptation_stats.hessian] = fmincon(optim_fun,transform([params.phi;params.tau_phi]),[],[],[],[],transform([1e-10 2e-3]),transform([1e1 1e6]),[],options);                        
        [params.phi,adaptation_stats.phi]=deal(itransform(adaptation_stats.beta(1)));
        [params.tau_phi,adaptation_stats.tau_phi] = deal(itransform(adaptation_stats.beta(2)));            
        adaptation_stats.se=sqrt(diag(inv(adaptation_stats.hessian)));
        adaptation_stats.phi_range = itransform(adaptation_stats.beta(1) +[-1 1]*adaptation_stats.se(1));
        adaptation_stats.tau_phi_range = itransform(adaptation_stats.beta(2) +[-1 1]*adaptation_stats.se(2));   
        adaptation_stats.within_stream=params.within_stream;
        X = X_update_fun(params.phi,params.tau_phi);          
    end
    options = statset('MaxIter',params.maxIter);   
    if params.useGPU
        X=gpuArray(X);
    end
    [~,dev,stats] = glmfit(X,Y,params.distribution,'options',options);            
    stats.dev=dev;
    if nargout>1
        Yhat=gather(params.link.Inverse(stats.beta(1)+X*stats.beta(2:end)));
    end
    if params.fit_adaptation
        stats.adaptation_stats=adaptation_stats;
    else
        switch params.distribution
            case 'poisson'
                NLL = -sum(log(poisspdf(Y,Yhat)));
            case 'normal'
                NLL = -sum(log(normpdf(Y,Yhat,1)));
        end 
        stats.NLL=NLL;
    end
    fields = fieldnames(stats);
    nf=length(fields);
    for f=1:nf
        stats.(fields{f}) = gather(stats.(fields{f}));
    end    
end

function NLL = NLL_fun(X_update_fun,Y,phi,tau_phi,params)
    X = X_update_fun(phi,tau_phi);  
    options = statset('MaxIter',params.maxIter);                  
    if params.useGPU
        X=gpuArray(X);
    end        
    beta = glmfit(X,Y,params.distribution,'options',options);       
    pred = params.link.Inverse(beta(1)+X*beta(2:end));    
    switch params.distribution
        case 'poisson'
            NLL = -sum(log(poisspdf(Y,pred)));
        case 'normal'
            NLL = -sum(log(normpdf(Y,pred,1)));
    end 
    NLL=gather(NLL);
end

function stop = optim_status_fun(x,optimValues,state,time_at_start,itransform)
    stop=false;
    switch state
        case 'iter'
            if isempty(optimValues.stepsize)
                if isempty(optimValues.firstorderopt)
                    fprintf('   %2d           %3d         %10.10e                                             %3.3e      %4.1f         %6.1f\n',...
                        optimValues.iteration,optimValues.funccount,optimValues.fval,...
                        itransform(x(1)),1000*itransform(x(2)),toc(time_at_start));
                else
                    fprintf('   %2d           %3d         %10.10e                         %3.3e           %3.3e      %4.1f         %6.1f\n',...
                        optimValues.iteration,optimValues.funccount,optimValues.fval,optimValues.firstorderopt,...
                       itransform(x(1)),1000*itransform(x(2)),toc(time_at_start));                        
                end
            else
                fprintf('   %2d           %3d         %10.10e     %3.3e           %3.3e           %3.3e      %4.1f         %6.1f\n',...
                  optimValues.iteration,optimValues.funccount,optimValues.fval,optimValues.stepsize,optimValues.firstorderopt,...
                  itransform(x(1)),1000*itransform(x(2)),toc(time_at_start));                
            end
        case 'interrupt'
              % Probably no action here. Check conditions to see  
              % whether optimization should quit.
        case 'init'
              fprintf('\nIteration    Func Count    NLL per timepoint     Step Size      1st Order Optimality       Phi         Tau (ms)     Time Elapsed (s)\n');
        case 'done'
              % Cleanup of plots, guis, or final plot
    otherwise
    end
end