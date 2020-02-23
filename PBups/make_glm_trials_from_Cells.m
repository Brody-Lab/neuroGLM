function rawData = make_glm_trials_from_Cells(Cells,varargin)
    % Makes the rawData trial structure used by neuroGLM from a Cells
    % file.
    % The final structure has fields giving the time of various events on
    % each trial and the spike times for each unit (those fields are called
    % "sptrain1",'sptrain2", etc
    %units in s
    %% parse and validate inputs
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('ref_event','cpoke_in',@(x)validateattributes(x,{'char'},{'nonempty'}));
    p.addParameter('samplingFreq',1e3,@(x)validateattributes(x,{'numeric'},{'scalar'})); %Hz
    p.addParameter('removeViolations',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('removeStimTrials',true,@(x)validateattributes(x,{'logical'},{'scalar'})); 
    p.addParameter('nClickBins',3,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
    p.addParameter('separate_stereo_click',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('separate_clicks_by','latency',@(x)validatestring(x,{'latency','time'}));    
    p.parse(varargin{:});
    % time relative to reference event over which you include spikes (make sure the window 
    % over which you include spiking data in your input structure is the same or smaller)    
    if isfield(Cells,'kSpikeWindowS')
        p.addParameter('spikeWindowS',Cells.kSpikeWindowS.(p.Results.ref_event),@(x)validateattributes(x,{'numeric'},{'numel',2})); 
    else
        p.addParameter('spikeWindowS',[-1.245 3.25],@(x)validateattributes(x,{'numeric'},{'numel',2}));         
    end
    p.parse(varargin{:});
    params = p.Results;
    fields_to_copy = {'rat','sess_date','sess_time','mat_file_name','sessid'};
    for f=1:length(fields_to_copy)
        if isfield(Cells,fields_to_copy{f})
            rawData.param.(fields_to_copy{f}) = Cells.(fields_to_copy{f});
        else
            warning('Could not copy field %s because it does not exist.',fields_to_copy{f});
        end
    end
    events = {'cpoke_in','cpoke_out','spoke'};
    rawData.timings = events;
    duration = round(diff(params.spikeWindowS)*params.samplingFreq);
    if params.removeViolations
        exclude_trials = Cells.Trials.violated;
    end
    if params.removeStimTrials
        exclude_trials = exclude_trials | Cells.Trials.laser.isOn;
    end   
    % remove side LED and free choice trials by default
    exclude_trials = exclude_trials | abs(Cells.Trials.gamma) > 90;
    %% preallocate structure in memory
    rawData.trial = struct();
    rawData.trial(sum(~exclude_trials)).duration = duration; % preallocate
    rawData.nTrials = sum(~exclude_trials);
    switch params.separate_clicks_by
        case 'time'
            right = linspace(0,params.samplingFreq*max([Cells.Trials.stim_dur_s_theoretical]),params.nClickBins+1);      
        case 'latency'
            total_rate=40; % this information isn't in the cells file now so make a good assumption
            binEdges = [0 expinv((1:params.nClickBins-1)/params.nClickBins,params.samplingFreq/total_rate) Inf];
    end
    good_trials = find(~exclude_trials);
    for t = 1:length(good_trials)
        original_t = good_trials(t);
        trial_start_time = Cells.Trials.stateTimes.(params.ref_event)(original_t) + params.spikeWindowS(1);
        rawData.trial(t).duration = duration;
        rawData.trial(t).ishit = Cells.Trials.is_hit(original_t);
        for e=1:length(events)
            rawData.trial(t).(events{e}) = Cells.Trials.stateTimes.(events{e})(original_t) - trial_start_time;                    
            if isnan(rawData.trial(t).(events{e}))
                rawData.trial(t).(events{e}) = [];
            else
                rawData.trial(t).(events{e}) = round(rawData.trial(t).(events{e})*params.samplingFreq);
            end
        end
        rawData.trial(t).gamma = Cells.Trials.gamma(original_t);
        if isfield(Cells.Trials,'stim_dur_s_theoretical')
            rawData.trial(t).stim_dur = Cells.Trials.stim_dur_s_theoretical(original_t);
        else
            rawData.trial(t).stim_dur = Cells.Trials.stim_dur_s(original_t);            
        end
        rawData.trial(t).pokedR = Cells.Trials.pokedR(original_t);
        for c = 1:length(Cells.spike_time_s.(params.ref_event))
            rawData.trial(t).(['sptrain',num2str(c)]) = round( (Cells.spike_time_s.(params.ref_event){c}{original_t} - params.spikeWindowS(1) ) * params.samplingFreq);
        end
        %% clicks (stereo click gets its own covariate, and depending on the value of nclickbins
        this_trial_left_clicks = params.samplingFreq*(Cells.Trials.stateTimes.left_clicks(Cells.Trials.stateTimes.left_click_trial==original_t) - trial_start_time); 
        this_trial_right_clicks = params.samplingFreq*(Cells.Trials.stateTimes.right_clicks(Cells.Trials.stateTimes.right_click_trial==original_t) - trial_start_time);   
        if strcmp(params.separate_clicks_by,'latency')
            all_clicks = [this_trial_left_clicks this_trial_right_clicks];
            click_idx = [zeros(1,numel(this_trial_left_clicks)),ones(1,numel(this_trial_right_clicks))];
            [all_clicks,sort_idx] = sort(all_clicks);
            click_idx = click_idx(sort_idx);
            itis = [Inf,diff(all_clicks)];
            if itis(2) == 0
                itis(2)=Inf;
            else
                warning('No stereo click at start of this trial!');
            end                
            this_trial_left_latencies = itis(click_idx==0);
            this_trial_right_latencies = itis(click_idx==1);
        end
        if params.separate_stereo_click
            if this_trial_left_clicks(1)==this_trial_right_clicks(1)
                rawData.trial(t).stereo_click = this_trial_left_clicks(1);
                this_trial_left_clicks(1)=[];
                this_trial_right_clicks(1)=[];
                if strcmp(params.separate_clicks_by,'latency')
                    this_trial_left_latencies(1)=[];
                    this_trial_right_latencies(1)=[];
                end                    
            else
                warning('No stereo click at start of this trial!');
            end
            if t==1
                rawData.timings = union(rawData.timings,{'stereo_click'});
            end            
        end
        for i=1:params.nClickBins
            if params.nClickBins==1
                idx_str = '';
            else
                idx_str = num2str(i);
            end
            if i==1
                if strcmp(params.separate_clicks_by,'time')
                    if params.separate_stereo_click
                        these_bin_edges = binEdges + rawData.trial(t).stereo_click;
                    else
                        these_bin_edges = binEdges + min([this_trial_left_clicks(:);this_trial_right_clicks(:)]);                    
                    end
                end
            end
            if t==1
                rawData.timings = union(rawData.timings,{['left_clicks',idx_str] ['right_clicks',idx_str] ['all_clicks',idx_str]});
            end
            switch params.separate_clicks_by
                case 'time'
                    rawData.trial(t).(['left_clicks',idx_str])  = this_trial_left_clicks(this_trial_left_clicks>=these_bin_edges(i) & this_trial_left_clicks<these_bin_edges(i+1)) ;       
                    rawData.trial(t).(['right_clicks',idx_str])  = this_trial_right_clicks(this_trial_right_clicks>=these_bin_edges(i) & this_trial_right_clicks<these_bin_edges(i+1)) ; 
                    rawData.trial(t).(['all_clicks',idx_str])  = union(rawData.trial(t).(['left_clicks',idx_str]),rawData.trial(t).(['right_clicks',idx_str])) ;        
                case 'latency'
                    rawData.trial(t).(['left_clicks',idx_str])  = this_trial_left_clicks(this_trial_left_latencies>=binEdges(i) & this_trial_left_latencies<binEdges(i+1)) ;       
                    rawData.trial(t).(['right_clicks',idx_str])  = this_trial_right_clicks(this_trial_right_latencies>=binEdges(i) & this_trial_right_latencies<binEdges(i+1)) ; 
                    rawData.trial(t).(['all_clicks',idx_str])  = union(rawData.trial(t).(['left_clicks',idx_str]),rawData.trial(t).(['right_clicks',idx_str])) ;                       
            end
        end
    end
    rawData.param.samplingFreq = params.samplingFreq; % Hz
    rawData.param.aligned_to = params.ref_event;
    rawData.param.ncells = length(Cells.spike_time_s.(params.ref_event));
end