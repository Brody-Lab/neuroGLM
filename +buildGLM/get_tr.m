function t = get_tr(tr,idx)
    % reconstruct indices of a covariate
    if nargin==1
        t = tr.lims(1):tr.binsize:tr.lims(2);
    else
        t  = tr.lims(1) + tr.binsize*(idx-1);
    end
end