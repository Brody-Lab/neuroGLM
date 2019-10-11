function w_samples = sampleWeights(dm, w , wcov, varargin)
    % Sample the weights per covariate given the weights per design matrix
    % column and their full covariance
    
    % Because the full covariance has been taken into account, you can do
    % meaningful non-parametric statistics on the samples.
    
    % If you provide a cell array of covariate names as a variable input
    % argument, it only provides samples for those

    % AGB 2019
    %
    % Input
    %   dm: design matrix structure
    %   w: weight on the basis functions
    %   wcov: weight covariance matrix
    %
    % Output
    %   wsamples.(label) = combined samples

    p=inputParser;
    p.addRequired('wcov',@(x)validateattributes(x,{'numeric'},{'size',[1 1].*numel(w)}));
    p.addParameter('covariates',{dm.dspec.covar.label},@(x)validateattributes(x,{'cell'},{}));
    p.addParameter('nsamples',1e4,@(x)validateattributes(x,{'numeric'},{'>',1,'finite','nonnan','scalar'}));
    p.parse(wcov,varargin{:});
    params=p.Results;
    w=w(:)';
    dspec = dm.dspec;
    try
        chol(params.wcov);
    catch
        error('weight covariance matrix must be positive definite');
    end
    covar_idx = ismember({dm.dspec.covar.label},params.covariates);
    dm.dspec.covar = dm.dspec.covar(covar_idx);
    if nargout>1 && ~covSupplied
        error('You must supply a weight covariance matrix as the third argument.');
    end

    if isfield(dm, 'biasCol') % undo z-score operation
        w(dm.biasCol) = [];
        if covSupplied
            params.wcov(dm.biasCol,:)=[];
            params.wcov(:,dm.biasCol)=[];
        end
    end

    if isfield(dm, 'zscore') % undo z-score operation
        w = (w .* dm.zscore.sigma(:)) + dm.zscore.mu(:);
        if covSupplied
            params.wcov = (params.wcov .* dm.zscore.sigma(:)) + dm.zscore.mu(:);
        end
    end

    if isfield(dm, 'constCols') % put back the constant columns
        w2 = zeros(dm.dspec.edim, 1);
        w2(~dm.constCols) = w; % first term is bias
        w = w2;
        if covSupplied
            wcov2 = zeros(dm.dspec.edim);
            wcov2(~dm.constCols,~dm.constCols)=params.wcov;
            params.wcov=wcov2;
        end
    end

    if numel(w) ~= dm.dspec.edim
        error('Expecting w to be %d dimension but it''s [%d]', ...
        dspec.edim, numel(w));
    end

    if numel(params.wcov) ~= dm.dspec.edim^2
        error('Expecting w to have %d^2 elements but it''s [%d]', ...
        dspec.edim, sqrt(numel(params.wcov)));
    end

    startIdx = [1 (cumsum([dspec.covar(:).edim]) + 1)];

    if covSupplied
        w=mvnrnd(w,params.wcov,params.nsamples);
    end
    for kCov = 1:numel(dspec.covar)
        covar = dspec.covar(kCov);
        basis = covar.basis;
        assert(isstruct(basis), 'Basis structure is not a structure?');
        sdim = covar.edim / basis.edim;
        for sIdx = 1:sdim
            w_sub = w(:,startIdx(kCov) + (1:basis.edim)-1 + basis.edim * (sIdx - 1));
            w_samples.(covar.label) = w_sub*basis.B';    
        end
    end
end
