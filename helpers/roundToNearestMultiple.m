function rounded = roundToNearestMultiple(x,step)
    % round to nearest multiple of an arbitrary value
    if ~isnumeric(step) || ~isnumeric(x)
        error('All inputs must be numeric.')
    elseif isempty(x) || all(isnan(x(:)))
        rounded=x;
    elseif all(isinf(x(:)))
        error('nearest multiple undefined for infinite x');
    elseif any(~isfinite(step(:)))
        error('step must be finite and non-nan.');
    else
        nans=~isfinite(x);
        rounded=round(x/step)*step;
        rounded(nans)=NaN;
    end
end