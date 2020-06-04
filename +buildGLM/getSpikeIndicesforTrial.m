function ndx = getSpikeIndicesforTrial(expt,trialIndex)
    trialT = cumsum(expt.binfun([expt.trial.duration]));
    ndx=[];
    if islogical(trialIndex)
        trialIndex = find(trialIndex(:)');
    end
    for t=trialIndex(:)'
        if t==1
            ndx = [ndx 1:trialT(t)];    
        else
            ndx = [ndx trialT(t-1)+1:trialT(t)];                
        end
    end
end