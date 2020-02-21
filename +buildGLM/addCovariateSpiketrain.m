function dspec = addCovariateSpiketrain(dspec, covLabel, stimLabel, desc, basisStruct, varargin)

if nargin < 4 || isempty(desc); desc = covLabel; end

if nargin < 5
    if dspec.expt.binSize<=2
        basisStruct = basisFactory.makeNonlinearRaisedCos(6, dspec.expt.binSize, [0 50], 1);        
    elseif dspec.expt.binSize==3
        basisStruct = basisFactory.makeNonlinearRaisedCos(5, dspec.expt.binSize, [0 50], 1);
    elseif dspec.expt.binSize<=6
        basisStruct = basisFactory.makeNonlinearRaisedCos(4, dspec.expt.binSize, [0 50], 1);        
    elseif dspec.expt.binSize<=12
        basisStruct = basisFactory.makeNonlinearRaisedCos(3, dspec.expt.binSize, [0 50], 1);                
    elseif dspec.expt.binSize<=50
        basisStruct = basisFactory.makeNonlinearRaisedCos(2, dspec.expt.binSize, [0 50], 1);        
    else
        error('Case of a single spike history basis not implemented yet. Are you sure you want to be doing this?');
    end
        
end

assert(ischar(desc), 'Description must be a string');

offset = 1; % Make sure to be causal. No instantaneous interaction allowed.

binfun = dspec.expt.binfun;
stimHandle = @(trial, expt) basisFactory.deltaStim(binfun(trial.(stimLabel)), binfun(trial.duration));

dspec = buildGLM.addCovariate(dspec, covLabel, desc, stimHandle, basisStruct, offset, varargin{:});