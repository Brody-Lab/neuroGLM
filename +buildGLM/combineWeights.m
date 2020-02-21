function [wout,wvarout] = combineWeights(dm, w , wcov)
    % Combine the weights per column in the design matrix per covariate 

    % wvarout is a structure organized the same as wout but with the estimated
    % variances of the recombined fitted weights. This requires supplying a
    % weight covariance matrix (wcov) used to resample weights from a
    % corresponding Gaussian.

    % AGB 2019
    %
    % Input
    %   dm: design matrix structure
    %   w: weight on the basis functions
    %
    % Output
    %   wout.(label).data = combined weights
    %   wout.(label).tr = time axis

    nsamples=1e4;
    w=w(:)';
    dspec = dm.dspec;
    binSize = dspec.expt.binSize;
    if nargin>2
        covSupplied=true;
        assert(size(wcov,1)==numel(w) && size(wcov,2)==numel(w)); %
    else
        covSupplied=false;
    end

    if nargout>1 && ~covSupplied
        error('You must supply a weight covariance matrix as the third argument.');
    end

    if isfield(dm, 'biasCol') % undo z-score operation
        wout.bias = w(dm.biasCol);
        w(dm.biasCol) = [];
        if covSupplied
            wcov(dm.biasCol,:)=[];
            wcov(:,dm.biasCol)=[];
        end
    end

    if isfield(dm, 'zscore') % undo z-score operation
        w = (w .* dm.zscore.sigma(:)) + dm.zscore.mu(:);
        if covSupplied
            wcov = (wcov .* dm.zscore.sigma(:)) + dm.zscore.mu(:);
        end
    end

    if isfield(dm, 'constCols') % put back the constant columns
        w2 = zeros(dm.dspec.edim, 1);
        w2(~dm.constCols) = w; % first term is bias
        w = w2;
        if covSupplied
            wcov2 = zeros(dm.dspec.edim);
            wcov2(~dm.constCols,~dm.constCols)=wcov;
            wcov=wcov2;
        end
    end

    if numel(w) ~= dm.dspec.edim
        error('Expecting w to be %d dimension but it''s [%d]', ...
        dspec.edim, numel(w));
    end

    if numel(wcov) ~= dm.dspec.edim^2
        error('Expecting w to have %d^2 elements but it''s [%d]', ...
        dspec.edim, sqrt(numel(wcov)));
    end

    startIdx = [1 (cumsum([dspec.covar(:).edim]) + 1)];
    wout = struct();

    if covSupplied
        w=mvnrnd(w,wcov,nsamples);
    end
    for kCov = 1:numel(dspec.covar)
        covar = dspec.covar(kCov);
        basis = covar.basis;
        if isempty(basis)
            w_sub = w(:,startIdx(kCov) + (1:covar.edim) - 1);
            wout.(covar.label).tr = ((1:size(w_sub, 2))-1 + covar.offset/binSize) * binSize;
            wout.(covar.label).data = mean(w_sub*basis.B',1);
        else
            assert(isstruct(basis), 'Basis structure is not a structure?');
            sdim = covar.edim / basis.edim;
            wout.(covar.label).data = zeros(size(basis.B, 1), sdim);
            for sIdx = 1:sdim
                w_sub = w(:,startIdx(kCov) + (1:basis.edim)-1 + basis.edim * (sIdx - 1));
                wout.(covar.label).data = mean(w_sub*basis.B',1);    
                wout.(covar.label).tr = (basis.tr(:, 1) + covar.offset/binSize) * binSize * ones(1, sdim);    
                if covSupplied
                    wvarout.(covar.label).cov = cov(w_sub*basis.B',1);                
                    wvarout.(covar.label).data = diag(wvarout.(covar.label).cov);
                    wvarout.(covar.label).tr=(basis.tr(:, 1) + covar.offset/binSize) * binSize * ones(1, sdim);
                end
            end
        end
    end
end
