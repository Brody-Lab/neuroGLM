function dspec = addCovariateSpiketrain(dspec, covLabel, stimLabel, desc, basisStruct, varargin)

if nargin < 4 || isempty(desc); desc = covLabel; end

bin_size_ms =  1e3 * dspec.expt.binSize / dspec.expt.param.samplingFreq;
max_lag_ms=50;


% you should check these filter banks make sense given a different lag
% (pretty robust to bin_size)


if nargin < 5
    if max_lag_ms==100
        if bin_size_ms==1
            basisStruct = basisFactory.makeNonlinearRaisedCos(10, bin_size_ms, [0 max_lag_ms], 1);        
        elseif bin_size_ms==2
            basisStruct = basisFactory.makeNonlinearRaisedCos(9, bin_size_ms, [0 max_lag_ms], 1);
        elseif bin_size_ms<=4
            basisStruct = basisFactory.makeNonlinearRaisedCos(7, bin_size_ms, [0 max_lag_ms], 1);        
        elseif bin_size_ms==5
            basisStruct = basisFactory.makeNonlinearRaisedCos(6, bin_size_ms, [0 max_lag_ms], 1);   
        elseif bin_size_ms<=7
            basisStruct = basisFactory.makeNonlinearRaisedCos(5, bin_size_ms, [0 max_lag_ms], 1);        
        elseif bin_size_ms<=10
            basisStruct = basisFactory.makeNonlinearRaisedCos(4, bin_size_ms, [0 max_lag_ms], 1);  
        elseif bin_size_ms<=20
            basisStruct = basisFactory.makeNonlinearRaisedCos(3, bin_size_ms, [0 max_lag_ms], 1);          
        elseif bin_size_ms<=max_lag_ms
            basisStruct = basisFactory.makeNonlinearRaisedCos(2, bin_size_ms, [0 max_lag_ms], 1);        
        else
            error('Case of a single spike history basis not implemented yet. Are you sure you want to be doing this?');
        end
    elseif max_lag_ms==50
        if bin_size_ms<=2
            basisStruct = basisFactory.makeNonlinearRaisedCos(6, bin_size_ms, [0 max_lag_ms], 1);        
        elseif bin_size_ms==3
            basisStruct = basisFactory.makeNonlinearRaisedCos(5, bin_size_ms, [0 max_lag_ms], 1);
        elseif bin_size_ms<=6
            basisStruct = basisFactory.makeNonlinearRaisedCos(4, bin_size_ms, [0 max_lag_ms], 1);        
        elseif bin_size_ms<=12
            basisStruct = basisFactory.makeNonlinearRaisedCos(3, bin_size_ms, [0 max_lag_ms], 1);                
        elseif bin_size_ms<=max_lag_ms
            basisStruct = basisFactory.makeNonlinearRaisedCos(2, bin_size_ms, [0 max_lag_ms], 1);  
        else
            error('Case of a single spike history basis not implemented yet. Are you sure you want to be doing this?');            
        end
    end      
end

assert(ischar(desc), 'Description must be a string');

offset = 1; % Make sure to be causal. No instantaneous interaction allowed.

binfun = dspec.expt.binfun;
stimHandle = @(trial, expt) basisFactory.deltaStim(binfun(trial.(stimLabel)), binfun(trial.duration));

dspec = buildGLM.addCovariate(dspec, covLabel, desc, stimHandle, basisStruct, offset, varargin{:});