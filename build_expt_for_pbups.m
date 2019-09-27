function expt = build_expt_for_pbups(rawData,varargin)
    % rawData must be in units of ms
    expt = buildGLM.initExperiment('ms', 1, [rawData.param.rat,', ',rawData.param.sess_date],rawData.param);
    for f=1:length(rawData.timings)
        if isfield(rawData.trial,rawData.timings{f})
            expt = buildGLM.registerTiming(expt,rawData.timings{f},rawData.timings{f});
        end
    end
    for c=1:rawData.param.ncells
        expt = buildGLM.registerSpikeTrain(expt,['sptrain',num2str(c)],[rawData.param.rat,', ',rawData.param.sess_date,', cell ',num2str(c)]);
    end
    expt = buildGLM.registerValue(expt,'pokedR','chose right');
    expt = buildGLM.registerValue(expt,'ishit','correct trial');    
    expt.trial = rawData.trial;
end