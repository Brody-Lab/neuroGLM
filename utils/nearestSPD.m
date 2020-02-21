function Ahat = nearestSPD(A,varargin)
    % A must be a real,square matrix
    % for speed, if A is symmetric, use flag 'symmetric'
    % if you don't care if floating point trash makes Ahat fail the chol
    % test, use flag 'fast'
    
    % Adrian Bondy, 2015
    if ~any(strcmp(varargin,'symmetric'))
        % force symmetry
        A=(A+A')/2;
    end   
    %% fast polar decomposition
    [~,S,Q] = svd(A);  % Economy size. % faster if you assign a first output variable!
    H = Q*S*Q';   % H isn't necessarily symmetric due to floating point trash  
    Ahat=(A+H)./2;
    %% make sure it passes the chol test 
    if ~any(strcmp(varargin,'fast'))
        Ahat=(Ahat+Ahat')/2;        
        if ~isSPD(Ahat)
            Ahat = fixSPD(Ahat);
        end
    end
end    

function Ahat = fixSPD(A)
    Ahat = (A + A')/2;
    %% test that A is in fact PD. if it is not so, then tweak it just a bit.
    p = 1;
    k = 0;
    while p ~= 0
      [~,p] = chol(Ahat);
      k = k + 1;
      if k>1000 % If you've been tweaking with a ridge with no success for a long time, add a tiny random matrix, lives on ([-eps eps])
          mssg(0,'NearestSPD is having issues unrelated to floating point trash.  Adding tiny random matrix.');
          A=A+eps(size(A,1)).*rand(size(A,1));
          Ahat=nearestSPD(A);
          return
      end
      if p ~= 0
        % A failed the chol test. It must have been just a hair off,
        % due to floating point trash, so it is simplest now just to
        % tweak by adding a tiny multiple of an identity matrix.
        mineig = min(eig(A));
        if mineig==0
            Ahat = A + (-mineig*k.^2 + eps*min(diag(A)))*eye(size(A));
        else
            Ahat = A + (-mineig*k.^2 + eps(mineig))*eye(size(A));            
        end
      end
    end             
end