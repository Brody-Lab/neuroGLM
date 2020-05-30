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
    p.addParameter('use_parallel',false,@(x)validateattributes(x,{'logical'},{'scalar'}));    % parfor operates over cells unless cross-validation is used (i.e. kfold>1) in which case it operates over cross-validation folds
    p.addParameter('maxIter',25,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));
    p.addParameter('bin_size_s',0.001,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));  % resolution of the model. predictions have this bin size.  
    p.addParameter('minResponsiveFrac',0.5,@(x)validateattributes(x,{'numeric'},{'scalar','positive','<',1}));
    p.addParameter('minSpkParamRatio',10,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    p.addParameter('distribution','poisson',@(x)validateattributes(x,{'char'},{}));   
    p.addParameter('link','canonical',@(x)validatestring(x,{'log','identity','softplus'}));
    p.addParameter('save_path','');
    p.addParameter('useGPU',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('phi',0.1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'})); % adaptation param
    p.addParameter('tau_phi',0.1,@(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'})); % adaptation param
    p.addParameter('within_stream',false,@(x)validateattributes(x,{'logical'},{'scalar'})); % within_stream adaptation flag    
    p.addParameter('fit_adaptation',true,@(x)validateattributes(x,{'logical'},{'scalar'})); % whether or not to fit phi and tau_phi for each neuron using gradient descent
    p.parse(varargin{:});
    params=p.Results;
    validatestring(params.distribution,{'poisson','normal'},mfilename,'distribution');
    switch params.link
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
        else
            [a,b,~] = fileparts(params.save_path);            
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
    dspec_base = buildGLM.initDesignSpec(expt);    
    dspec_base = build_dspec_for_pbups(dspec_base,covariates,[]);
    dm = buildGLM.compileSparseDesignMatrix(dspec_base, 1:nTrials, params);  
    X_base = dm.X;
    ncells=rawData.param.ncells;
    %% loop over cells
    % do a non-parallel loop to compute things needed for each cell without
    % passing large necessary structures to workers
    responsive_enough=true(1,ncells);  
    n_params=sum([dspec_base.covar.edim])+1+6; % 6 is for spike history filter but this may be an overestimate. roughly good enough for these purposes.
    %% first (non-parallel) loop across cells
    % non-parallel loop across cells just to compute the design matrices
    % (to avoid passing a structure with all the spike times to the
    % workers)
    % this is a hack -- better rewriting could just extract the spike times
    % for each cell and pass it as a cell array. but this would require some coding.
    [X,Y] = deal(cell(1,length(params.cellno)));
    tic;fprintf('\nBuilding design matrices for responsive cells ...');
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
        dspec = build_dspec_for_pbups(dspec_base,'spike_history',params.cellno(c)); 
        %% build design matrices        
        dm = buildGLM.compileSparseDesignMatrix(dspec, 1:nTrials, params , X_base, numel(dspec.covar));   % only remake the columns associated with the spike history term since the rest is common across cells
        dm = buildGLM.removeConstantCols(dm);       
        Y{c} = full(buildGLM.getBinnedSpikeTrain(expt, ['sptrain',num2str(params.cellno(c))], dm.trialIndices)); 
        X{c} = dm.X;
    end
    fprintf(' took %s.\n%g of %g cells are sufficiently responsive. Fitting now:\n',timestr(toc),sum(responsive_enough),ncells);
    responsive_cells = find(responsive_enough);
    X = X(responsive_enough);
    Y = Y(responsive_enough);
    %% trial contains the list of ALL spike times, which is no longer needed  
    dm = rmfield(dm,'X');    
    trial_fields= fieldnames(dm.dspec.expt.trial);
    is_spk_field = strncmp(trial_fields,'sp',2);
    params.dm=dm;        
    params.dm.dspec.expt.trial = rmfield(dm.dspec.expt.trial,trial_fields(is_spk_field)); 
    %% loop over cells
    if params.use_parallel && (isempty(params.kfold) || params.kfold<2)
        parfor c=1:sum(responsive_enough)
            fprintf('Cell id %g (%g of %g to fit):\n',params.cellno(responsive_cells(c)),c,sum(responsive_enough));           
            stats(c) = mainLoop(X{c},Y{c},params);
        end    
    else
        for c=sum(responsive_enough):-1:1
            fprintf('Cell id %g (%g of %g to fit):\n',params.cellno(responsive_cells(c)),c,sum(responsive_enough));           
            stats(c) = mainLoop(X{c},Y{c},params);
        end
    end
    for c=1:sum(responsive_enough)
        stats(c).cellno = params.cellno(responsive_cells(c));
        stats(c).covariate_stats.cellno= stats(c).cellno;
        covariate_stats(c) = stats(c).covariate_stats;
    end    
    %% save
    params.dm=dm; % get back list of all spike times
    params.rat = Cells.rat;
    params.sess_date = Cells.sess_date;
    params.sessid = Cells.sessid;
    if params.save && isfield(stats,'dev')
        if isempty(params.save_path)
            mat_file_name = strrep(mat_file_name,'glmfits_save_test','glmfits');
            save(mat_file_name,'stats','params','-v7.3'); % v7.3 needed because often larger than 2GB
            fprintf('Saved fit stats successfully to %s.\n',mat_file_name);
            mat_file_name = strrep(mat_file_name,'glmfits','glmfits_summary');
            save(mat_file_name,'covariate_stats','params','-v7'); % v7.3 needed because often larger than 2GB
            fprintf('Saved fit summary successfully to %s.\n',mat_file_name);            
        else
            mat_file_name = params.save_path;
            save(mat_file_name,'stats','params','-v7.3');
            fprintf('Saved fit stats successfully to %s.\n',mat_file_name);            
            [a,b,c] = fileparts(mat_file_name);
            mat_file_name = [a,filesep,b,'_summary',c];
            save(mat_file_name,'covariate_stats','params','-v7'); % v7.3 needed because often larger than 2GB
            fprintf('Saved fit summary successfully to %s.\n',mat_file_name);             
        end
    end
end

function stats = mainLoop(X,Y,params)   
    % all the main fitting function needs is the binned observations (Y),
    % the design matrix (X, no bias column), the dm structure (which is the same
    % for all cells) and the parent function params
    % if z-scoring was performed on the design matrix, you'd need to have a
    % separate dm structure for each cell storing this
    tic;fprintf('   Fitting UN cross-validated model ... ');drawnow; 
    nTrials = numel(params.dm.trialIndices);
    stats = fit(X, Y, params, 1:nTrials);
    fields = fieldnames(stats);
    rm_fields = {'residp','residd','resida'};
    nf=length(fields);
    for f=1:nf
        if ismember(fields{f},rm_fields)
            stats=rmfield(stats,fields{f});
        else
            stats.(fields{f}) = gather(stats.(fields{f}));
        end
    end     
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
            [train_idx,test_idx,init_beta] = deal(cell(1,params.kfold));                
            for i=1:params.kfold
                train_idx{i} = getSpkIdxFun(stats.cvp.training(i));
                test_idx{i} = getSpkIdxFun(stats.cvp.test(i));
            end
            [cv_stats,dev] = zeros(1,params.kfold);
            if params.use_parallel
                parfor i=1:params.kfold 
                    [~, dev(i), cv_stats(i)] = fit(X(train_idx{i},:),Y(train_idx{i}), params, stats.cvp.training(i));  
                end    
            else
                for i=params.kfold:-1:1
                    [~, dev(i), cv_stats(i)] = fit(X(train_idx{i},:),Y(train_idx{i}), params, stats.cvp.training(i));  
                end                  
            end
            cv_stats = rmfield(cv_stats,{'resid','residp','residd','resida'}); % these take up A TON of space, and could always be generated if needed
            for i=1:params.kfold
                cv_stats(i).dev = dev(i);
                fields = fieldnames(cv_stats);
                for f=1:length(fields)
                    cv_stats(i).(fields{f}) = gather(cv_stats(i).(fields{f}));
                end                    
                cv_stats(i).Yhat=glmval(cv_stats(i).beta,gather(X(test_idx{i},:)),'log',cv_stats(i));
                cv_stats(i).init_beta = init_beta{i};
                cv_stats(i).Yhat_init_beta = gather([ones(size(X(test_idx{i},:),1),1),X(test_idx{i},:)]*cv_stats(i).init_beta);                    
                [cv_stats(i).ws,cv_stats(i).wvars] =combineWeightFun(cv_stats(i).beta,cv_stats(i).covb);
            end
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

function stats = fit(X,Y,params,trials)
    if params.fit_adaptation
        global time_at_start itransform
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
        covar_idx=find(ismember({params.dm.dspec.covar.label},{'left_clicks','right_clicks'}));  
        X_update_fun = @(phi,tau_phi)buildGLM.compileSparseDesignMatrix(params.dm.dspec, trials,...
            struct('phi',phi,'tau_phi',tau_phi,'within_stream',params.within_stream) ,X, covar_idx);  % only remake the columns associated with the click terms since the rest is independent of the adaptation parameters
        time_at_start=tic;              
        options=optimoptions('fmincon','UseParallel',false,'OutputFcn',@optim_status_fun,'Algorithm','active-set');  % stop if you are taking tiny steps but not if the change in LL is small -- sometimes the gradient is really small far from the optimum and you should keep going until you get near the basin              
        optim_fun = @(x)NLL_fun(X_update_fun,Y,itransform(x(1)),itransform(x(2)),rmfield(params,'dm'));                        
        [adaptation_stats.beta,adaptation_stats.NLL,adaptation_stats.exitflag,adaptation_stats.output,adaptation_stats.lambda,...
            adaptation_stats.grad,adaptation_stats.hessian] = fmincon(optim_fun,transform([params.phi;params.tau_phi]),[],[],[],[],transform([1e-10 2e-3]),transform([1e1 1e6]),[],options);                        
        [params.phi,adaptation_stats.phi]=deal(adaptation_stats.beta(1));
        [params.tau_phi,adaptation_stats.tau_phi] = deal(adaptation_stats.beta(2));            
        adaptation_stats.se=sqrt(diag(inv(adaptation_stats.hessian)));
        adaptation_stats.phi_range = (adaptation_stats.beta(1) +[-1 1]*adaptation_stats.se(1));
        adaptation_stats.tau_phi_range = (adaptation_stats.beta(2) +[-1 1]*adaptation_stats.se(2));     
        dm = X_update_fun(params.phi,params.tau_phi);          
        X=dm.X;
    end
    options = statset('MaxIter',params.maxIter);   
    if params.useGPU
        X=gpuArray(X);
    end
    [~,dev,stats] = glmfit(X,Y,params.distribution,'options',options);            
    stats.dev=dev;
    stats.Yhat=glmval(stats.beta,gather(X),params.link);    
    stats.adaptation_stats=adaptation_stats;
end

function NLL = NLL_fun(X_update_fun,Y,phi,tau_phi,params)
    dm = X_update_fun(phi,tau_phi);  
    options = statset('MaxIter',params.maxIter);                  
    if params.useGPU
        dm.X=gpuArray(dm.X);
    end        
    beta = glmfit(dm.X,Y,params.distribution,'options',options);       
    pred = params.link.Inverse(beta(1)+dm.X*beta(2:end));    
    switch params.distribution
        case 'poisson'
            NLL = -sum(log(poisspdf(Y,pred)));
        case 'normal'
            NLL = -sum(log(normpdf(Y,pred,1)));
    end 
    NLL=gather(NLL);
end

function stop = optim_status_fun(x,optimValues,state)
    stop=false;
    global time_at_start itransform
    switch state
        case 'iter'
            if isempty(optimValues.stepsize)
                if isempty(optimValues.firstorderopt)
                    fprintf('   %2d            %3d        %10.10e                                             %3.3e      %4.1f         %6.1f\n',...
                        optimValues.iteration,optimValues.funccount,optimValues.fval,...
                        itransform(x(1)),1000*itransform(x(2)),toc(time_at_start));
                else
                    fprintf('   %2d            %3d        %10.10e                         %3.3e           %3.3e      %4.1f         %6.1f\n',...
                        optimValues.iteration,optimValues.funccount,optimValues.fval,optimValues.firstorderopt,...
                       itransform(x(1)),1000*itransform(x(2)),toc(time_at_start));                        
                end
            else
                fprintf('   %2d            %3d        %10.10e     %3.3e           %3.3e           %3.3e      %4.1f         %6.1f\n',...
                  optimValues.iteration,optimValues.funccount,optimValues.fval,optimValues.stepsize,optimValues.firstorderopt,...
                  itransform(x(1)),1000*itransform(x(2)),toc(time_at_start));                
            end
        case 'interrupt'
              % Probably no action here. Check conditions to see  
              % whether optimization should quit.
        case 'init'
              fprintf('\nIteration      Func Count         NLL            Step Size      1st Order Optimality       Phi         Tau (ms)     Time Elapsed (s)\n');
        case 'done'
              % Cleanup of plots, guis, or final plot
    otherwise
    end
end
