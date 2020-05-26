function dm = compileSparseDesignMatrix(dspec, trialIndices, adaptation_params, X, covar_idx)
    % Compile information from experiment according to given DesignSpec
    expt = dspec.expt;
    subIdxs = buildGLM.getGroupIndicesFromDesignSpec(dspec);
    trialT = expt.binfun([expt.trial(trialIndices).duration]);
    totalT = sum(trialT);
    if nargin<4
        dm.X = zeros(totalT, sum([dspec.covar.edim]));
        if nargin<5
            covar_idx = 1:numel(dspec.covar);
        end
    else
        dm.X = X;
    end
    trialIndices = trialIndices(:)';
    left_click_idx={dspec.covar.label}=="left_clicks";
    right_click_idx={dspec.covar.label}=="right_clicks";
    stereo_click_idx={dspec.covar.label}=="stereo_click";
    if nargin>2 &&  abs(adaptation_params.phi-1)>eps && abs(adaptation_params.tau_phi)>eps
        adapt=true;
    else
        adapt=false;
    end
    for k = trialIndices
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
        else
             right_clicks = full(dspec.covar(right_click_idx).stim(expt.trial(k)));
             left_clicks = full(dspec.covar(left_click_idx).stim(expt.trial(k)));
        end
        for kCov = covar_idx % for each covariate
            covar = dspec.covar(kCov);
            sidx = subIdxs{kCov};
            if isfield(covar, 'cond') && ~isempty(covar.cond) && ~covar.cond(expt.trial(k))
                continue;
            end
            switch covar.label
                case 'left_clicks'
                    stim = left_clicks;
                case 'right_clicks'
                    stim = right_clicks;
                otherwise
                    stim = covar.stim(expt.trial(k)); % either dense or sparse
                    stim = full(stim);
            end
            if isfield(covar, 'basis') && ~isempty(covar.basis)
                switch covar.basis.type
                    case 'makeNonlinearRaisedCos'
                        dm.X(ndx, sidx) = basisFactory.convBasis(stim, covar.basis, covar.offset);
                    case 'raised cosine@makeSmoothTemporalBasis'
                        dm.X(ndx, sidx) = basisFactory.convBasis(stim, covar.basis, covar.basis.param.binfun(covar.offset));
                end
            else
                dm.X(ndx, sidx) = stim; % adaptation not applied here, since without a basis set these covariates don't extend in time and therefore are not subject to adaptation
            end
        end
    end

    dm.trialIndices = trialIndices;
    dm.dspec = dspec;

    %% Check sanity of the design
    if any(~isfinite(dm.X(:)))
        warning('Design matrix contains NaN or Inf...this is not good!');
    end
end
