function dspec = addCovariateSpiketrain(dspec, covLabel, stimLabel, desc, basisStruct, varargin)

if nargin < 4 || isempty(desc); desc = covLabel; end

bin_size_ms =  1e3 * dspec.expt.binSize / dspec.expt.param.samplingFreq;
max_lag_s=0.05;


% you should check these filter banks make sense given a different lag
% (pretty robust to bin_size)


if nargin < 5
    if max_lag_s==0.1
        if bin_size_ms==1
            basisStruct = basisFactory.makeNonlinearRaisedCos(10,  [0 max_lag_s], 0.001, dspec.expt);        
        elseif bin_size_ms==2
            basisStruct = basisFactory.makeNonlinearRaisedCos(9,  [0 max_lag_s], 0.001, dspec.expt);
        elseif bin_size_ms<=4
            basisStruct = basisFactory.makeNonlinearRaisedCos(7,  [0 max_lag_s], 0.001, dspec.expt);        
        elseif bin_size_ms==5
            basisStruct = basisFactory.makeNonlinearRaisedCos(6,  [0 max_lag_s], 0.001, dspec.expt);   
        elseif bin_size_ms<=7
            basisStruct = basisFactory.makeNonlinearRaisedCos(5,  [0 max_lag_s], 0.001, dspec.expt);        
        elseif bin_size_ms<=10
            basisStruct = basisFactory.makeNonlinearRaisedCos(4,  [0 max_lag_s], 0.001, dspec.expt);  
        elseif bin_size_ms<=20
            basisStruct = basisFactory.makeNonlinearRaisedCos(3,  [0 max_lag_s], 0.001, dspec.expt);          
        elseif bin_size_ms<=max_lag_s
            basisStruct = basisFactory.makeNonlinearRaisedCos(2,  [0 max_lag_s], 0.001, dspec.expt);        
        else
            error('Case of a single spike history basis not implemented yet. Are you sure you want to be doing this?');
        end
    elseif max_lag_s==0.05
        if bin_size_ms<=2
            basisStruct = basisFactory.makeNonlinearRaisedCos(6,  [0 max_lag_s], 0.001, dspec.expt);        
        elseif bin_size_ms==3
            basisStruct = basisFactory.makeNonlinearRaisedCos(5,  [0 max_lag_s], 0.001, dspec.expt);
        elseif bin_size_ms<=6
            basisStruct = basisFactory.makeNonlinearRaisedCos(4,  [0 max_lag_s], 0.001, dspec.expt);        
        elseif bin_size_ms<=12
            basisStruct = basisFactory.makeNonlinearRaisedCos(3,  [0 max_lag_s], 0.001, dspec.expt);                
        elseif bin_size_ms<=max_lag_s
            basisStruct = basisFactory.makeNonlinearRaisedCos(2,  [0 max_lag_s], 0.001, dspec.expt);  
        else
            error('Case of a single spike history basis not implemented yet. Are you sure you want to be doing this?');            
        end
    end      
end

assert(ischar(desc), 'Description must be a string');

offset = basisStruct.tr(2); % Make sure to be causal. No instantaneous interaction allowed.

binfun = dspec.expt.binfun;
stimHandle = @(trial, expt) basisFactory.deltaStim(binfun(trial.(stimLabel)), binfun(trial.duration));

dspec = buildGLM.addCovariate(dspec, covLabel, desc, stimHandle, basisStruct, offset, varargin{:});