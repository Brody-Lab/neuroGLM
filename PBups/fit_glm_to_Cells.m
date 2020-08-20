function [stats,params] = fit_glm_to_Cells(Cells,varargin)
    % fits a GLM for spiking data recorded during PBups and contained within
    % an Cells data structure.
    % This function is essentially a wrapper for the neuroGLM package forked from the
    % Pillow lab.
    %% parse and validate inputs
    p=inputParser;
    p.addParameter('cellno',[]);
    p.addParameter('kfold',1,@(x)validateattributes(x,{'numeric'},{'scalar','integer','>',0}));
    p.addParameter('save',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('create_pool',false,@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.addParameter('parallelize_by_cells',true,@(x)validateattributes(x,{'logical'},{'scalar'}));        
    p.addParameter('n_workers',1,@(x)validateattributes(x,{'numeric'},{'scalar','integer','>',0}));    % parallelization operates over cells unless cross-validation is used (i.e. kfold>1) in which case it operates over cross-validation folds
    p.addParameter('maxIter',100,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));
    p.addParameter('lambda',0.05,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'})); % amplitude of ridge penalty applied in the GLM fitting
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
    p.addParameter('choice_time_back_s',0.75,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'})); % choice kernels extend backwards acausally in time before stimulus end by this many seconds
    p.addParameter('include_mono_clicks',true,@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.addParameter('save_by_cell',true,@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.addParameter('fit',true,@(x)validateattributes(x,{'logical'},{'scalar'})); % setting to false is useful when you just want the params structure
    p.addParameter('save_params',true,@(x)validateattributes(x,{'logical'},{'scalar'})); % setting to false useful when you are simultaneously running jobs across cells. you don't want params to overwrite each other.
    p.addParameter('test_cv',false,@(x)validateattributes(x,{'logical'},{'scalar'})); % go straight to cross-validation for testing purposes 
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
        fields=fieldnames(Cells);
        if length(fields)==1
            Cells=Cells.(fields{1}); % if not saved with -struct flag
        end
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
            save_subdir = [datestr(now,'YYYY_mm_dd_HH_MM_SS'),'_glmfits'];
            params.save_path = fullfile(fileparts(Cells.mat_file_name),save_subdir);
        end
        if ~isdir(params.save_path)
            mkdir(params.save_path);
        end
        if params.save_by_cell
            spikes_dir=fullfile(params.save_path,'spikes');
            stats_dir=fullfile(params.save_path,'stats');
            if ~isdir(spikes_dir)
                mkdir(spikes_dir);
            end
            if ~isdir(stats_dir)
                mkdir(stats_dir);            
            end
        end        
        if params.save_by_cell
            mat_file_name = fullfile(params.save_path,'spikes','glmfits_save_test.mat');
        else
            mat_file_name = fullfile(params.save_path,'glmfits_save_test.mat');            
        end
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
    dspec_base = build_dspec_for_pbups(dspec_base,covariates,[],'choice_time_back_s',params.choice_time_back_s,'include_mono_clicks',params.include_mono_clicks);
    dm = buildGLM.compileSparseDesignMatrix(dspec_base, 1:nTrials, params);
    X_base = dm.X;
    ncells=rawData.param.ncells;
    responsive_enough=true(1,ncells);  
    n_params=sum([dspec_base.covar.edim])+1+6; % 6 is for spike history filter but this may be an overestimate. roughly good enough for these purposes.
    %% determine if cells are responsive enough to fit    
    for c=length(params.cellno):-1:1
        responsiveFrac(c) = sum(arrayfun(@(x)numel(x.(['sptrain',num2str(params.cellno(c))])),expt.trial)>0)./nTrials;         
        totalSpikes(c) = length(cat(1,expt.trial.(['sptrain',num2str(params.cellno(c))])));
        if params.minResponsiveFrac>0
            if responsiveFrac(c) < params.minResponsiveFrac || totalSpikes(c)<n_params*params.minSpkParamRatio        
                responsive_enough(c)=false;
            end
        end     
    end    
    %% build design matrices for each cell to fit
    % non-parallel loop across cells just to compute the design matrices
    % (to avoid passing a structure with all the spike times to the
    % workers).
    % this is a bit of a hack -- you could just extract the spike times
    % for each cell and pass it as a cell array. but this would require
    % significant rewriting of code.
    X = cell(1,length(params.cellno));
    tic;fprintf('\nBuilding design matrices for responsive cells ...');
    spikes=struct(); 
    for c=length(params.cellno):-1:1      
        dspec = build_dspec_for_pbups(dspec_base,'spike_history',params.cellno(c)); 
        dm = buildGLM.compileSparseDesignMatrix(dspec, 1:nTrials, params , numel(dspec.covar));   % only remake the columns associated with the spike history term since the rest is common across cells
        dm.X(:,1:(dspec.edim-dspec.covar(end).edim))=X_base;
        dm = buildGLM.removeConstantCols(dm);       
        spikes(c).Y = full(buildGLM.getBinnedSpikeTrain(expt, ['sptrain',num2str(params.cellno(c))], dm.trialIndices)); 
        spikes(c).cellno=params.cellno(c);
        X{c} = dm.X;
        if ~params.fit
           break % let it compute dm for one cell so it can be stored in the params
        end          
    end
    if ~params.fit
        stats=struct([]);
        fprintf(' took %s.\nSkipping fitting.\n',timestr(toc));        
    elseif ~any(responsive_enough)
        stats=struct([]);
        fprintf(' took %s.\n%g of %g cells are sufficiently responsive. Skipping fitting.\n',timestr(toc),sum(responsive_enough),ncells);        
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
    if params.fit
        if params.n_workers>1 && params.parallelize_by_cells
            if params.create_pool
                delete(gcp('nocreate'));
                parpool(params.n_workers);
            end
            parfor c=1:sum(responsive_enough)
                cellno(c)=params.cellno(responsive_cells(c));                 
                fprintf('Cell id %g (%g of %g to fit):\n',cellno(c),c,sum(responsive_enough));           
                [stats(c),spikes(c).Yhat,spikes(c).Yhat_cv] = mainLoop(X{c},spikes(c).Y,params);
            end    
        else
            for c=sum(responsive_enough):-1:1
                cellno(c)=params.cellno(responsive_cells(c));                             
                fprintf('Cell id %g (%g of %g to fit):\n',cellno(c),c,sum(responsive_enough));           
                [stats(c),spikes(c).Yhat,spikes(c).Yhat_cv] = mainLoop(X{c},spikes(c).Y,params);
            end
        end  
        for c=sum(responsive_enough):-1:1
            stats(c).cellno= cellno(c); % add in separate loop to avoid assignment between dissimilar structures           
        end
    end
    %% save
    params.rat = Cells.rat;
    if isfield(Cells,'sess_date')
        params.sess_date = Cells.sess_date;
    end
    params.sessid = Cells.sessid;
    params.responsive_enough=responsive_enough;
    params.responsiveFrac=responsiveFrac;
    params.totalSpikes=totalSpikes;
    [params.git_branch,params.git_commit] = return_git_status(fileparts(which(mfilename)));
    [~,params.hostname] = system('hostname');
    params.hostname=deblank(params.hostname);
    params.save_time=datestr(now,'YYYY_mm_DD_HH_MM_SS');
    if params.save               
        if params.save_by_cell
            if params.save_params
                save(fullfile(params.save_path,'glmfit_params.mat'),'params','-v7');  
                fprintf('Saved %s successfully.\n',fullfile(params.save_path,'glmfit_params.mat'));                            
            end
            for i=1:length(stats)
                these_stats=stats(i);
                these_spikes=spikes(i);
                save(fullfile(params.save_path,'stats',sprintf('glmfit_stats_cell%g.mat',stats(i).cellno)),'-struct','these_stats');
                fprintf('Saved %s successfully.\n',fullfile(params.save_path,'stats',sprintf('glmfit_stats_cell%g.mat',stats(i).cellno)));                
                save(fullfile(params.save_path,'spikes',sprintf('glmfit_spikes_cell%g.mat',stats(i).cellno)),'-struct','these_spikes');      
                fprintf('Saved %s successfully.\n',fullfile(params.save_path,'spikes',sprintf('glmfit_spikes_cell%g.mat',stats(i).cellno)));                                
            end
        else
            save(fullfile(params.save_path,'glmfit_stats.mat'),'params','stats','-v7');
            fprintf('Saved %s successfully.\n',fullfile(params.save_path,'glmfit_stats.mat'));
            save(fullfile(params.save_path,'glmfit_spikes.mat'),'params','spikes','-v7');
            fprintf('Saved %s successfully.\n',fullfile(params.save_path,'glmfit_spikes.mat'));        
        end
    end
end

function [stats,Yhat,Yhat_cv] = mainLoop(X,Y,params)   
    % all the main fitting function needs is the binned observations (Y),
    % the design matrix (X, no bias column), the dm structure (which is the same
    % for all cells) and the parent function params
    % if z-scoring was performed on the design matrix, you'd need to have a
    % separate dm structure for each cell storing this
    nTrials = numel(params.dm.trialIndices);
    params.dm.biasCol=1;  
    if params.test_cv
        stats.badly_scaled=false;
        Yhat=[];
    else
        tic;fprintf('   Fitting UN cross-validated model ... ');drawnow;         
        [stats,Yhat] = fit(X, Y, params, 1:nTrials);          
    end
    Yhat_cv = zeros(size(Y));
    fprintf('took %s.\n',timestr(toc));
    % Fit cross-validated model (if requested and if uncross-validated fit was not badly scaled)
    if ~isempty(params.kfold) && params.kfold>1
        if stats.badly_scaled
            fprintf('Skipping cross-validation since fit to all data was badly scaled.\n');
            stats.cvp=[];
            stats.cv_stats=[];
        else
            stats.cvp = cvpartition(nTrials,'KFold',params.kfold);
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
                    idx{i} = stats.cvp.training(i);
                    Ys{i}=Y(train_idx{i});
                end
                parfor i=1:params.kfold 
                    [cv_stats(i),~,Yhat_cv_cell{i}] = fit(Xs{i},Ys{i}, params, idx{i});  
                end    
                for i=1:params.kfold
                    Yhat_cv(test_idx{i}) = Yhat_cv_cell{i};
                end
            else
                for i=params.kfold:-1:1
                    [cv_stats(i),~,Yhat_cv(test_idx{i})] = fit(X(train_idx{i},:),Y(train_idx{i}), params, stats.cvp.training(i));  
                end                  
            end       
            stats.cv_stats=cv_stats;
            fprintf('Took %s.\n',timestr(toc));            
        end
    end
    %stats.covariate_stats = get_covariate_stats(stats,params);
%     fields = fieldnames(stats.wvars);
%     for f=1:length(fields) % you can remove covariance structure now that summary has been computed
%         stats.wvars.(fields{f}) = rmfield(stats.wvars.(fields{f}),'cov');
%     end
end

function [stats,Yhat,Yhat_test] = fit(X,Y,params,trials)
    if params.useGPU
        Y=gpuArray(Y);
    end  
    if ~all(trials)
        cv=true;
        getSpkIdxFun = @(trial_idx)buildGLM.getSpikeIndicesforTrial(params.dm.dspec.expt,trial_idx);                
        train_idx = find(getSpkIdxFun(trials));
        test_idx = find(getSpkIdxFun(~trials));
        Y_test = Y(test_idx);
        X_test = X(test_idx,:);
        X = X(train_idx,:);
        Y = Y(train_idx);
    else
        cv=false;
    end
    glm_options = statset('MaxIter',params.maxIter);       
    if params.fit_adaptation
        [params.transform,params.itransform] = deal(@(x)x);        
        %transform=@(x)log(1+x); % identity near zero, but logarithmic as x->Inf
        %itransform=@(x)max(0,exp(x)-1);
        % makes fits converge faster when tau_phi is large
        % I tried tanh but this needs to be scaled (by ~2e4) to allow large tau_phi
        % and this leads to slow convergence when tau_phi is small.
        % init
        %params.phi=0.5;
        %params.tau_phi=0.3; % good starting points but let user decide
        covar_idx=find(ismember({params.dm.dspec.covar.label},{'left_clicks','right_clicks'}));        
        X_update_fun = @(phi,tau_phi)buildGLM.updateSparseDesignMatrix_covar(params.dm.dspec, trials, struct('phi',phi,'tau_phi',tau_phi,'within_stream',params.within_stream), covar_idx, X);
        time_at_start=tic;      
        if cv
            X_test_update_fun = @(phi,tau_phi)buildGLM.updateSparseDesignMatrix_covar(params.dm.dspec, ~trials, struct('phi',phi,'tau_phi',tau_phi,'within_stream',params.within_stream), covar_idx, X_test);                    
            return_cov=false;
        else      
            X_test_update_fun=[];
            Y_test=[];
            return_cov=true;
        end
        options=optimoptions('fmincon','UseParallel',false,'OutputFcn',@optim_status_fun,'Algorithm','interior-point','FunctionTolerance',eps,'StepTolerance',1e-8);  % stop if you are taking tiny steps but not if the change in LL is small -- sometimes the gradient is really small far from the optimum and you should keep going until you get near the basin              
        optim_fun = @(x)NLL_fun(params.itransform(x(1)),params.itransform(x(2)));                        
        [adaptation_stats.beta,adaptation_stats.NLL,adaptation_stats.exitflag,adaptation_stats.output,adaptation_stats.lambda,...
            adaptation_stats.grad,adaptation_stats.hessian] = fmincon(optim_fun,params.transform([params.phi;params.tau_phi]),[],[],[],[],params.transform([1e-10 2e-3]),params.transform([1e1 5e1]),[],options);                        
        [params.phi,adaptation_stats.phi]=deal(params.itransform(adaptation_stats.beta(1)));
        [params.tau_phi,adaptation_stats.tau_phi] = deal(params.itransform(adaptation_stats.beta(2)));            
        adaptation_stats.se=sqrt(diag(inv(adaptation_stats.hessian)));
        adaptation_stats.phi_range = params.itransform(adaptation_stats.beta(1) +[-1 1]*adaptation_stats.se(1));
        adaptation_stats.tau_phi_range = params.itransform(adaptation_stats.beta(2) +[-1 1]*adaptation_stats.se(2));   
        adaptation_stats.within_stream=params.within_stream;
        X = X_update_fun(params.phi,params.tau_phi);          
    end
    if params.useGPU
        X=gpuArray(X);
    end
    [~,dev,stats] = glmfit(X,Y,params.distribution,'options',glm_options,'lambda',params.lambda);            
    stats.dev=dev;
    Yhat=gather(params.link.Inverse(stats.beta(1)+X*stats.beta(2:end)));       
    if params.fit_adaptation
        stats.adaptation_stats=adaptation_stats;
        stats.NLL = compute_NLL(Y,Yhat,params.distribution,true);
        adaptation_stats = rmfield(adaptation_stats,'NLL');
    else
        stats.NLL = compute_NLL(Y,Yhat,params.distribution,true);
    end
    if cv
        if params.fit_adaptation
            X_test = X_test_update_fun(params.phi,params.tau_phi);          
        end     
        Yhat_test=gather(params.link.Inverse(stats.beta(1)+X_test*stats.beta(2:end)));   
        stats.NLL_test = compute_NLL(Y_test,Yhat_test,params.distribution,true);
        stats.NLL_train = stats.NLL;
        stats=rmfield(stats,'NLL');        
    end
    fields = fieldnames(stats);
    nf=length(fields);
    for f=1:nf
        stats.(fields{f}) = gather(stats.(fields{f}));
    end  
    [stats.ws,stats.wvars]=buildGLM.combineWeights(params.dm,params.dm.dspec, stats.beta , stats.covb,return_cov);

    %% this is the function being minimized
    function NLL = NLL_fun(phi,tau_phi)
        global beta
        thisX = X_update_fun(params.itransform(phi),params.itransform(tau_phi));  
        if params.useGPU
            thisX=gpuArray(thisX);
        end        
        [beta,~,these_stats] = glmfit(thisX,Y,params.distribution,'options',glm_options,'lambda',params.lambda);     
        if 1/these_stats.cond<eps || these_stats.badly_scaled % disallow adaptation params which cause numerical instability in GLM fitting
            NLL=Inf;
            return
        end
        pred = params.link.Inverse(beta(1)+thisX*beta(2:end));    
        NLL = compute_NLL(Y,pred,params.distribution,false); % has to be false. if you do by timepoint here, the LL surface is too flat for the algorithm (i.e. scaling matters!)
    end
    
    %% this nested function is the fmincon output function which prints and stores information about each iteration
    function stop = optim_status_fun(x,optimValues,state)
        global beta
        optimValues.phi=x(1);
        optimValues.tau_phi=x(2);
        stop=false;
        switch state
            case {'iter','done'}
                if ~isempty(Y_test) && any(Y_test)
                    X_test = X_test_update_fun(params.itransform(optimValues.phi),params.itransform(optimValues.tau_phi));  
                    pred_test = params.link.Inverse(beta(1)+X_test*beta(2:end));
                    optimValues.NLL_test_per_timepoint = compute_NLL(Y_test,pred_test,params.distribution,true);
                else
                    optimValues.NLL_test_per_timepoint=NaN;
                end
                optimValues.NLL_train_per_timepoint = optimValues.fval ./ size(Y,1);
                optimValues.state=state;
                adaptation_stats.iter_info(optimValues.iteration+1) = optimValues;
                if isempty(optimValues.stepsize)
                    if isempty(optimValues.firstorderopt)
                        fprintf('   %2d           %3d            %10.10e             %10.10e         %3.3e      %4.1f         %6.1f\n',...
                            optimValues.iteration,optimValues.funccount,optimValues.NLL_train_per_timepoint,optimValues.NLL_test_per_timepoint,...
                            params.itransform(optimValues.phi),1000*params.itransform(optimValues.tau_phi),toc(time_at_start));
                    else
                        fprintf('   %2d           %3d            %10.10e             %10.10e         %3.3e           %3.3e      %4.1f         %6.1f\n',...
                            optimValues.iteration,optimValues.funccount,optimValues.NLL_train_per_timepoint,optimValues.NLL_test_per_timepoint,optimValues.firstorderopt,...
                           params.itransform(optimValues.phi),1000*params.itransform(optimValues.tau_phi),toc(time_at_start));                        
                    end
                else
                        fprintf('   %2d           %3d            %10.10e             %10.10e         %3.3e           %3.3e           %3.3e      %4.1f         %6.1f\n',...
                      optimValues.iteration,optimValues.funccount,optimValues.NLL_train_per_timepoint,optimValues.NLL_test_per_timepoint,optimValues.stepsize,optimValues.firstorderopt,...
                      params.itransform(optimValues.phi),1000*params.itransform(optimValues.tau_phi),toc(time_at_start));                
                end
                if state=="done"
                    fprintf('...done.\n');
                end
            case 'interrupt'
                  % Probably no action here. Check conditions to see  
                  % whether optimization should quit.
            case 'init'
                  fprintf('\nIteration    Func Count    NLL per timepoint (train)    NLL per timepoint (test)     Step Size      1st Order Optimality       Phi         Tau (ms)     Time Elapsed (s)\n');
        end
    end

end

function NLL = compute_NLL(Y,Y_hat,distribution,per_timepoint)
    if per_timepoint
        func=@mean;
    else
        func=@sum;
    end
    switch distribution
        case 'poisson'
            NLL = gather(-func(log(poisspdf(Y,Y_hat))));
        case 'normal'
            NLL = gather(-func(log(normpdf(Y,Y_hat,1))));
        otherwise
            error('Unrecognized distribution.');
    end
end