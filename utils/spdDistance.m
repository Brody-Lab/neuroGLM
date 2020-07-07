function distance = spdDistance(A,varargin)
    if issparse(A)
        A=full(A);
    end
    if ~any(strcmp(varargin,'symmetric'))
        % force symmetry
        A=(A+A')/2;
    end   
    [~,S,Q] = svd(A);  % Economy size. % faster if you assign a first output variable!
    H = Q*S*Q';   % H isn't necessarily symmetric due to floating point trash  
    distance= norm( (A-H) ,'fro')/2;       
end