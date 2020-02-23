function rounded = nearestEven(x)
    % round to nearest even integer
    if ~isnumeric(x) || ~isreal(x) || any(isinf(x(:)))
        error('X must be real and numeric and finite.')
    elseif isempty(x) || all(isnan(x(:))) 
        rounded=x;
    else
        rounded=roundToNearestMultiple(x,2);
    end
end