function stats = fit_glm_to_Cells(Cells,varargin)
    % fits a GLM for spiking data recorded during PBups and contained within
    % an Cells data structure.
    % This function is essentially a wrapper for the neuroGLM package forked from the
    % Pillow lab.
    %% parse and validate inputs
    p=inputParser;
    p.addParameter('cellno',[]);
    p.addParameter('kfold',[]);
    p.KeepUnmatched=true;
    p.addParameter('save',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('use_parallel',false,@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.addParameter('maxIter',25,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));
    p.addParameter('minResponsiveFrac',0,@(x)validateattributes(x,{'numeric'},{'scalar','positive','<',1}));
    p.addParameter('minSpkParamRatio',0,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    p.addParameter('separate_clicks_by_side',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('distribution','poisson',@(x)validatestring(x,{'poisson','normal'}));
    p.parse(varargin{:});
    params=p.Results;
    %% make rawData and expt structures (i.e. put event times and spike times into neuroGLM format)
    rawData = make_glm_trials_from_Cells(Cells,varargin{:}); 
    expt=build_expt_for_pbups(rawData); 
    nTrials = rawData.nTrials;
    if isempty(params.cellno)
        params.cellno=1:length(Cells.tt);
        if isempty(params.cellno)
            warning('No cells in Cells!');
            return
        end
    end
    fprintf('\n');
    %% testing save
    if params.save
        fprintf('Testing save now before spending a long time fitting... ');
        mat_file_name = strrep(Cells.mat_file_name,'Cells','glmfits_save_test');
        test=[];
        save(mat_file_name,'test','-v7.3');
        delete(mat_file_name);   
        fprintf(' Success!\n');
    end    
    cell_count=0;    
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
    n_cells = length(params.cellno);
    for c=1:n_cells
        S=struct();
        S.cellno = params.cellno(c);          
        fprintf('Cell id %g (%g of %g to fit):\n',S.cellno,c,length(params.cellno));
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
        %% determine if spike/parameter ratio is acceptable
        S.totalSpks = sum(Y);
        S.spkParamRatio = S.totalSpks ./ (size(dm.X,2)+1);   
        if params.minSpkParamRatio>0
            if S.spkParamRatio < params.minSpkParamRatio
                fprintf('Cell %g only has %g spikes to %g params to be fit. Moving on without fitting.\n',S.cellno,S.totalSpks,size(dm.X,2)+1);
                fields =fieldnames(S);
                for f=1:length(fields)
                    stats(c).(fields{f}) = S.(fields{f});
                end
                continue
            end
        end        
        %% Fitting UN cross-validated model
        tic;fprintf('   Fitting UN cross-validated model ... ');                    
        [~, S.dev, stat_temp] = glmfit(dm.X, Y, params.distribution);
        %tic;B=pinv(X'*X) * (X'*Y);toc;
        fields_to_copy = fieldnames(stat_temp);
        nf=length(fields_to_copy);
        for f=1:nf
            S.(fields_to_copy{f}) = stat_temp.(fields_to_copy{f});
        end       
        % compute model predicted firing rates
        switch params.distribution
            case 'normal'
                link = 'identity';
            case 'poisson'
                link = 'log';
        end
        S.Yhat=glmval(S.beta,dm.X,link,S.beta);
        fprintf('took %s.\n',timestr(toc));
        % reconstruct fitted kernels by weighted combination of basis functions
        [S.ws,S.wvars] = buildGLM.combineWeights(buildGLM.addBiasColumn(dm), S.beta , S.covb);
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
                [Xtrain,Ytrain,Xtest] = deal(cell(1,params.kfold));
                for i=1:params.kfold
                    train_idx = getSpkIdxFun(S.cvp.training(i));
                    test_idx = getSpkIdxFun(S.cvp.test(i));
                    Xtrain{i} = dm.X(train_idx,:);
                    Ytrain{i} = Y(train_idx);
                    Xtest{i} = dm.X(test_idx,:);
                end
                options = statset('MaxIter',params.maxIter);
                [cv_stats,dev]=deal([]);
                if params.use_parallel
                    for i=1:params.kfold
                        [~, dev(i), cv_stats(i)] = glmfit(Xtrain{i}, Ytrain{i}, 'poisson','options',options);  
                    end    
                else
                    for i=1:params.kfold
                        [~, dev(i), cv_stats(i)] = glmfit(Xtrain{i}, Ytrain{i}, 'poisson','options',options);  
                    end                  
                end
                cv_stats = rmfield(cv_stats,{'resid','residp','residd','resida'}); % these take up A TON of space, and could always be generated if needed
                for i=1:params.kfold
                    cv_stats(i).dev=dev(i);
                    cv_stats(i).Yhat=glmval(cv_stats(i).beta,Xtest{i},'log',cv_stats(i));
                    [cv_stats(i).ws,cv_stats(i).wvars] =combineWeightFun(cv_stats(i).beta,cv_stats(i).covb);
                end
                S.cv_stats=cv_stats;
                fprintf('Took %s.\n',timestr(toc));            
            end
        end
        fields =fieldnames(S);
        for f=1:length(fields)
            stats(c).(fields{f}) = S.(fields{f});
        end
    end
    if params.save && isfield(stats,'dev')
        mat_file_name = strrep(Cells.mat_file_name,'Cells','glmfits');
        save(mat_file_name,'stats','-v7.3');
        fprintf('Saved fit stats successfully to %s.\n',mat_file_name);
    end
end