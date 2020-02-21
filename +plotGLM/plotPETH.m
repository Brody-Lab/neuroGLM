function [h,filtered_spikes] = plotPETH(expt,variable,trialIndices,varargin)
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('align_to','');
    p.addParameter('desired_delay',500);    
    p.addParameter('bin_size',10);
    p.addParameter('color',[0 0 0]);
    p.addParameter('plot',true);
    p.addParameter('sd',50);
    p.parse(varargin{:});
    params=p.Results;
    sd=params.sd;
    count=0;
    for k=trialIndices
        count=count+1;
        if ischar(variable)
            spikes = buildGLM.getBinnedSpikeTrain(expt, variable, k)*1000;
        else
            spikes = variable(buildGLM.getSpikeIndicesforTrial(expt,k))*1000;
        end
        filtered_spikes(count,:) = filterArray(full(spikes),my_gauss_kernel(sd*5,sd));
    end
    aligned_to = expt.param.aligned_to;
    event_time = expt.trial(1).(aligned_to); 
    if ~isempty(params.align_to)
        filtered_spikes = realign_spikes(expt,filtered_spikes,trialIndices,params.align_to,params.desired_delay);
        aligned_to = params.align_to;
        event_time = params.desired_delay;
    end
    min_trials =10;
    times = ((1:expt.trial(1).duration)-event_time)/1000;
    idx = sum(~isnan(filtered_spikes))>min_trials;
    if numel(trialIndices)>1
        if params.plot
            n_bins=100;
            [a,b] = discretize(times,n_bins);
            bin_centers=mean([b(1:end-1);b(2:end)]);
            for i=1:n_bins
                filtered_binned_spikes(:,i) = nanmean(filtered_spikes(:,find(a==i)),2);
            end
            h=shadedErrorBar(bin_centers,filtered_binned_spikes,{@nanmean,@SE});
            h.mainLine.Color = params.color;
            if ~isempty(h.patch)
                h.patch.FaceColor = params.color;
                h.patch.FaceAlpha = 0.1;
            end
            xlabel(['Time (s) relative to ',strrep(aligned_to,'_',' ')]);
            ylabel('Firing Rate (sp/s)');        
        end
    elseif params.plot
        h=plot(times,filtered_spikes,'color',params.color);
        xlabel(['Time (s) relative to ',strrep(aligned_to,'_',' ')]);
        ylabel('Firing Rate (sp/s)');
    end
    if ~params.plot
        h=[];
    end
end

function realigned_spikes = realign_spikes(expt,spikes,trialIndices,align_to,desired_delay)
    realigned_spikes = NaN(size(spikes));
    for i=1:length(trialIndices)
        event_time = expt.trial(trialIndices(i)).(align_to);
        start = event_time - desired_delay;
        start = max(1,start);
        finish = expt.trial(trialIndices(i)).duration;%max(max(expt.trial(trialIndices(i)).left_clicks),max(expt.trial(trialIndices(i)).right_clicks));
        idx = round(start):finish;
        realigned_spikes(i,1:length(idx)) = spikes(i,idx);
    end
end

