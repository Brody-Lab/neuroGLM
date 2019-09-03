function expt = build_expt_for_pbups(rawData,varargin)
    % rawData must be in units of ms
    timingFields = {'cpoke_in','cpoke_out','clicks_on','left_clicks','right_clicks','spoke','right_reward','left_reward'};
    expt = buildGLM.initExperiment('ms', 1, [rawData.param.rat,', ',rawData.param.sess_date],rawData.param);
    for f=1:length(timingFields)
        if isfield(rawData.trial,timingFields{f})
            expt = buildGLM.registerTiming(expt,timingFields{f},timingFields{f});
        end
    end
    for c=1:rawData.param.ncells
        expt = buildGLM.registerSpikeTrain(expt,['sptrain',num2str(c)],[rawData.param.rat,', ',rawData.param.sess_date,', cell ',num2str(c)]);
    end
    expt = buildGLM.registerValue(expt,'pokedR','chose right');
    expt.trial = rawData.trial;
end