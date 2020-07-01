function expt = build_expt_for_pbups(rawData,varargin)
    p=inputParser;
    p.addParameter('bin_size_s',0.001,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    p.parse(varargin{:});
    params=p.Results;
    if isfield(rawData.param,'sess_date')   
        id= [char(rawData.param.rat),', ',char(rawData.param.sess_date)];        
    elseif isfield(rawData.param,'sessid')
        id=rawData.param.sessid;
    elseif isfield(rawData.param,'mat_file_name')
        id=rawData.param.mat_file_name;
    end 
    expt = buildGLM.initExperiment(rawData.param.samplingFreq*params.bin_size_s, id,rawData.param);
    for f=1:length(rawData.timings)
        if isfield(rawData.trial,rawData.timings{f})
            expt = buildGLM.registerTiming(expt,rawData.timings{f},rawData.timings{f});
        end
    end
    for c=1:rawData.param.ncells
        expt = buildGLM.registerSpikeTrain(expt,['sptrain',num2str(c)],[num2str(id),', cell ',num2str(c)]);
    end
    expt = buildGLM.registerValue(expt,'pokedR','chose right');
    expt = buildGLM.registerValue(expt,'ishit','correct trial');
    expt.trial = rawData.trial;
end
