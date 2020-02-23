function p=empirical_p(x,distribution,tail)
    validateattributes(x,{'numeric'},{},'','x',1);
    validateattributes(distribution,{'numeric'},{},'','distribution',2);     
    if all(isnan(x))
        p=NaN(size(x));
        warning('No non-NaN elements of distribution');
        return
    end        
    xSz = size(x);
    %% repmat distribution if needed
    distribution=squeeze(distribution);
    distSz = size(distribution);
    if length(distSz)>length(xSz)
        if length(distSz)-length(xSz)==1
            xSz(end+1)=1;
        else
            error('distribution cannot have more than one dimension more than sample.');
        end
    end
    if any(xSz~=distSz & xSz~=1) % not fixable by bsxfun, need to repmat distSz
        if sum(distSz~=1)>1
            error('There does not appear to be a vector for each element of x.');
        end
        origSz=xSz;
        if sum(xSz>1)==1
            x=x(:);
            xSz=size(x);                        
            firstSingleton=2;
        else
            x=squeeze(x);
            xSz=size(x);            
            firstSingleton=length(xSz)+1;
            xSz(end+1)=1;
        end
        shape = ones(1,firstSingleton);
        shape(end) = numel(distribution);
        distribution=reshape(distribution,shape);
        repetitions = xSz;
        repetitions(firstSingleton) = 1;
        distribution = repmat(distribution,repetitions);
        distSz=size(distribution);          
    else
        origSz=xSz;
    end
    %%
    mismatch = xSz~=distSz;
    dim=find(mismatch);
    if any(xSz>1 & mismatch)
        error('Non-singleton dimensions of input arguments mustch match.');
    elseif numel(dim)>1
        error('Cannot match up a vector in distribution for each element of x.');
    end
    if ~any(mismatch)
        warning('Distribution is a scalar for each x. P-values will be arbitrarily set to NaN.');
        p = NaN(size(x));
        return   
    end
    nl = distSz(mismatch);
    if nl<100
        warning('Distribution vectors contain less than 100 elements. P-value may be noisy.');
    end
    if nargin<3
        tail='both';
    end
    nans = isnan(distribution);    
    switch tail
        case 'low'     
            if isscalar(x)
                GE = distribution>=x;
            else
                GE = bsxfun(@ge,distribution,x);
            end              
            p = (1+sum(GE,dim)) ./ sum(~nans,dim);
        case 'high'
            if isscalar(x)
                LE = distribution<=x;
            else
                LE = bsxfun(@le,distribution,x);
            end              
            p = (1+sum(LE,dim)) ./ sum(~nans,dim);
        case 'both'
            if isscalar(x)
                LE = distribution<=x;
                GE = distribution>=x;
            else
                LE = bsxfun(@le,distribution,x);
                GE = bsxfun(@ge,distribution,x);
            end              
            p = 2 * (1 + min ( sum(GE,dim) , sum(LE,dim) )) ./ sum(~nans,dim) ;
    end
    p( ~(p<1) | isnan(x) )=NaN;
    p=reshape(p,origSz);
end