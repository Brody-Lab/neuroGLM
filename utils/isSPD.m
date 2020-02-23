function [isSPD,distance] = isSPD(A,varargin)
    distance=0;
    try
        chol(A);
        isSPD=true;
    catch
        isSPD=false;      
        distance=spdDistance(A,varargin{:});
    end
end