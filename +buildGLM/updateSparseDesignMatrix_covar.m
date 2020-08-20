function X = updateSparseDesignMatrix_covar(dspec, trialIndices, adaptation_params, covar_idx, X)
    % Compile information from experiment according to given DesignSpec
    % AGB: Now you can choose to only compute the bits of the design matrix for a
    % subset of the covariates by supplying an index of the covariates you
    % want to compute (covar_idx). Useful when fitting the adaptation
    % params where the columns corresponding to the clicks need to be
    % recomputed repeatedly during fitting.
    
    expt = dspec.expt;
    expt.trial=expt.trial(trialIndices);
    subIdxs = buildGLM.getGroupIndicesFromDesignSpec(dspec);
    trialT = expt.binfun([expt.trial.duration]);
    if size(X,1)~=sum(trialT)
        error('cannot reconcile trial indices supplied with design matrix supplied.');
    end
    left_click_idx={dspec.covar.label}=="left_clicks";
    right_click_idx={dspec.covar.label}=="right_clicks";
    if nargin>2 &&  abs(adaptation_params.phi-1)>eps && abs(adaptation_params.tau_phi)>eps && any(ismember(covar_idx,find(left_click_idx|right_click_idx)))
        adapt=true;
    else
        adapt=false;
    end
    for k = 1:length(expt.trial)
        ndx = (sum(trialT(1:k))-(trialT(k)-1)):sum(trialT(1:k));
        if adapt
            % if adapting clicks, perform the adaptation at the sampling resolution
            % of the original data, not the resolution of the model.
            left_clicks = expt.trial(k).left_clicks ./ expt.param.samplingFreq;
            right_clicks = expt.trial(k).right_clicks ./ expt.param.samplingFreq;
            stereo_click = expt.trial(k).stereo_click ./ expt.param.samplingFreq;
            [left_adapted,right_adapted] = adapt_clicks(left_clicks,right_clicks,stereo_click,...
                adaptation_params.phi,adaptation_params.tau_phi,adaptation_params.within_stream);
            left_clicks = full(basisFactory.deltaStim(expt.binfun(left_clicks*expt.param.samplingFreq),trialT(k),left_adapted));
            right_clicks = full(basisFactory.deltaStim(expt.binfun(right_clicks*expt.param.samplingFreq),trialT(k),right_adapted));
        end
        for kCov = covar_idx % for each covariate
            covar = dspec.covar(kCov);
            sidx = subIdxs{kCov};
            if isfield(covar, 'cond') && ~isempty(covar.cond) && ~covar.cond(expt.trial(k))
                continue;
            end
            if adapt
                switch covar.label
                    case 'left_clicks'
                        stim = left_clicks;
                    case 'right_clicks'
                        stim = right_clicks;
                    otherwise
                        stim = covar.stim(expt.trial(k)); % either dense or sparse
                        stim = full(stim);
                end
            else
                stim = covar.stim(expt.trial(k)); % either dense or sparse
                stim = full(stim);                
            end
            if isfield(covar, 'basis') && ~isempty(covar.basis)
                if covar.offset==0
                    offset=0;
                else
                    offset=expt.binfun(expt.param.samplingFreq*covar.offset);
                end
                switch covar.basis.type
                    case 'makeNonlinearRaisedCos'
                        X(ndx, sidx) = basisFactory.convBasis(stim, covar.basis, offset); % offset should be in the base units of dspec.expt
                    case 'raised cosine@makeSmoothTemporalBasis'
                        X(ndx, sidx) = basisFactory.convBasis(stim, covar.basis, offset); % offset should be in the base units of dspec.expt
                end
            else
                X(ndx, sidx) = stim; % adaptation not applied here, since without a basis set these covariates don't extend in time and therefore are not subject to adaptation
            end
        end
    end
    
    %% Check sanity of the design
    if any(~isfinite(X(:)))
        warning('Design matrix contains NaN or Inf...this is not good!');
    end
end
