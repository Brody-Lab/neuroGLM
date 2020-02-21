function dm = compileSparseDesignMatrix(dspec, trialIndices)
    % Compile information from experiment according to given DesignSpec
    expt = dspec.expt;
    subIdxs = buildGLM.getGroupIndicesFromDesignSpec(dspec);
    trialT = expt.binfun([expt.trial(trialIndices).duration]);
    totalT = sum(trialT);
    dm.X = zeros(totalT, sum([dspec.covar.edim]));
    trialIndices = trialIndices(:)';
    for k = trialIndices    
        ndx = (sum(trialT(1:k))-(trialT(k)-1)):sum(trialT(1:k));    
        for kCov = 1:numel(dspec.covar) % for each covariate
            covar = dspec.covar(kCov);
            sidx = subIdxs{kCov}; 
            if isfield(covar, 'cond') && ~isempty(covar.cond) && ~covar.cond(expt.trial(k))
                continue;
            end       
            stim = covar.stim(expt.trial(k)); % either dense or sparse
            stim = full(stim);
            if isfield(covar, 'basis') && ~isempty(covar.basis)
                switch covar.basis.type
                    case 'makeNonlinearRaisedCos'
                        dm.X(ndx, sidx) = basisFactory.convBasis(stim, covar.basis, covar.offset);
                    case 'raised cosine@makeSmoothTemporalBasis'
                        dm.X(ndx, sidx) = basisFactory.convBasis(stim, covar.basis, covar.basis.param.binfun(covar.offset));                        
                end
            else
                dm.X(ndx, sidx) = stim;
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
