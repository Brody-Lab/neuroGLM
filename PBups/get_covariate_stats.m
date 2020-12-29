function [covariate_stats,stats] = get_covariate_stats(stats,params,varargin)
    p = inputParser;
    p.addParameter('save',false);
    p.addParameter('save_path','');
    p.parse(varargin{:});
    if nargin==0 % if no args, prompt user to locate glmfits file
        [file,dir] = uigetfile('*glmfits.mat','Choose a *glmfits.mat file');
        glmfits_file = fullfile(dir,file); 
        load(glmfits_file);
        p.Results.save_path = strrep(glmfits_file,'glmfits','glmfits_summary');
    elseif ischar(stats) % user specified a filename
        p.Results.save_path = strrep(stats,'glmfits','glmfits_summary');   
        fprintf('Loading %s.\n',stats);              
        load(stats);
    end
    for i=1:length(stats)
       covariate_stats(i) = get_covariate_stats_internal(stats(i),params.dm,params.dm.dspec,varargin{:});
    end
    for i=1:length(stats)
        if isfield(stats,'cellno')
            covariate_stats(i).cellno = stats(i).cellno;
        end
        covariate_stats(i).badly_scaled = stats(i).badly_scaled;
    end
    if p.Results.save
       if isempty(p.Results.save_path)
           warning('No save path specified.');
           return
       else
           save(p.Results.save_path,'covariate_stats');
           fprintf('Saved covariate_stats to %s.\n',p.Results.save_path);      
       end
    end
end




function [covariate_stats,stats]=get_covariate_stats_internal(stats,dm,dspec,varargin)
    p=inputParser;
    p.addParameter('ilink',@(x)exp(x),@(x)validateattributes(x,{'function_handle'},{''}));
    p.parse(varargin{:});
    params=p.Results;
    if ismember('ilink',p.UsingDefaults) && isfield(stats,'params') && isfield(stats.params,'distribution') && strcmp(stats.params.distribution,'poisson')
        params.ilink = @(x)exp(x);
    end  
    if isempty(stats.ws)
        covariate_stats=[];
        warning('Empty structure.');
        return
    end
    MIfun = @(x,y)(x-y)./(x+y);        
    DC = params.ilink(0);
    covars = {dspec.covar.label};
    if isfield(stats,'cellno')
        fprintf('get_covariate_stats: cell %g\n',stats.cellno);
    end
    if ~isfield(dspec.expt.param,'nClickBins')
        dspec.expt.param.nClickBins = sum(contains({dspec.covar.label},'left_click'));    
    end
    %% check that cov is in struct
    dm.biasCol=1;
    if ~isfield(stats.wvars.(covars{1}),'cov')
        fprintf('get_covariate_stats: adding covariance matrix at solution grouped by covariate ... ');a=tic;     
        [stats.ws,stats.wvars] = buildGLM.combineWeights(dm, dspec,stats.beta , stats.covb,true);
        fprintf('took %s.\n',timestr(toc(a)));
    end     
    if ~isfield(stats,'wsamples')
        stats.wsamples = buildGLM.sampleWeights(dm,stats.beta,stats.covb,'nsamples',1e3);                    
    end
    %% get some basic stats from all covariates
    for i=1:length(covars)
        data = stats.ws.(covars{i}).data;
        tr = buildGLM.get_tr(stats.ws.(covars{i}).tr);
        samples = stats.wsamples.(covars{i});
        covariate_stats.(covars{i}).average = mean(params.ilink(data));
        covariate_stats.(covars{i}).average_pval = empirical_p(DC,mean(params.ilink(samples),2)); 
        [~,peak_idx] = max(abs(params.ilink(data)-DC));
        covariate_stats.(covars{i}).max_deviation_time = buildGLM.get_tr(stats.ws.(covars{i}).tr,peak_idx);
        covariate_stats.(covars{i}).max_deviation = params.ilink(data(peak_idx));  
        covariate_stats.(covars{i}).max_deviation_pval = empirical_p(0,samples(:,peak_idx));
        ops = {@lt,@gt};
        for k=1:2
            idx = ops{k}(tr,0);
            data = stats.ws.(covars{i}).data(idx);
            covariate_stats.(covars{i}).prepost_idx{k} = find(idx);            
            if any(idx)
                covariate_stats.(covars{i}).average_prepost(k) = mean(params.ilink(data));
                covariate_stats.(covars{i}).average_pval_prepost(k) = empirical_p(DC,mean(params.ilink(samples(:,idx)),2));                     
            else
                covariate_stats.(covars{i}).average_prepost(k) =NaN;
                covariate_stats.(covars{i}).average_pval_prepost(k)=NaN;                
            end
        end
    end             
    %% left and right clicks    
    if isfield(stats.wsamples,'left_clicks3')
        sides = {'left_clicks3','right_clicks3'};
        left_allclicks = params.ilink((stats.wsamples.left_clicks1 + stats.wsamples.left_clicks2 + stats.wsamples.left_clicks3)/3);
        right_allclicks = params.ilink((stats.wsamples.right_clicks1 + stats.wsamples.right_clicks2 + stats.wsamples.right_clicks3)/3);        
        %% average across all clicks
        tmp = mean(mean(left_allclicks,2)>mean(right_allclicks,2));
        covariate_stats.clicks_allclicks_average_LR_MI = MIfun(mean(left_allclicks(:)),mean(right_allclicks(:)));     
        covariate_stats.clicks_allclicks_average_LR_MI_pval = 2*min(tmp,1-tmp);        
        if tmp>0.5
            covariate_stats.pref_click_by_allclicks_average = 'left';
        else
            covariate_stats.pref_click_by_allclicks_average = 'right';
        end        
        %% get max deviation statistics of average across all click kernels
        tr = buildGLM.get_tr(stats.ws.(sides{1}).tr);
        [~,peak_idx] = max(abs(params.ilink(mean(left_allclicks,1))-DC));
        covariate_stats.max_deviation_time_left_allclicks = tr(peak_idx);
        covariate_stats.max_deviation_left_allclicks = params.ilink(mean(left_allclicks(:,peak_idx)));
        covariate_stats.max_deviation_left_allclicks_pval = empirical_p(0,left_allclicks(:,peak_idx));
        [~,peak_idx] = max(abs(params.ilink(mean(right_allclicks,1))-DC));
        covariate_stats.max_deviation_time_right_allclicks = tr(peak_idx);
        covariate_stats.max_deviation_right_allclicks = params.ilink(mean(right_allclicks(:,peak_idx)));        
        covariate_stats.max_deviation_right_allclicks_pval = empirical_p(0,right_allclicks(:,peak_idx));        
        if covariate_stats.max_deviation_left_allclicks>covariate_stats.max_deviation_right_allclicks
            covariate_stats.pref_click_by_allclick_deviation = 'left';
            idx = find(tr==covariate_stats.max_deviation_time_left_allclicks);  
            
            
            covariate_stats.pref_click_max_deviation_allclicks = covariate_stats.max_deviation_left_allclicks;
            covariate_stats.pref_click_max_deviation_time = covariate_stats.max_deviation_time_left_allclicks;                
            covariate_stats.pref_click_max_deviation_pval = covariate_stats.max_deviation_left_allclicks_pval;       
            covariate_stats.nonpref_click_max_deviation_allclicks = covariate_stats.max_deviation_right_allclicks;
            covariate_stats.nonpref_click_max_deviation_time = covariate_stats.max_deviation_time_right_allclicks;                
            covariate_stats.nonpref_click_max_deviation_pval = covariate_stats.max_deviation_right_allclicks_pval;      
            
            covariate_stats.clicks_max_deviation_allclicks_LR_MI = MIfun(covariate_stats.max_deviation_left_allclicks ,params.ilink(mean(right_allclicks(:,idx)))) ;                       
            covariate_stats.clicks_max_deviation_allclicks_LR_MI_pval =    2* mean(left_allclicks(:,idx) < right_allclicks(:,idx));            
        else
            covariate_stats.pref_click_by_allclick_deviation = 'right';
            idx = find(tr==covariate_stats.max_deviation_time_right_allclicks);   
            covariate_stats.nonpref_click_max_deviation_allclicks = covariate_stats.max_deviation_left_allclicks;
            covariate_stats.nonpref_click_max_deviation_time = covariate_stats.max_deviation_time_left_allclicks;                
            covariate_stats.nonpref_click_max_deviation_pval = covariate_stats.max_deviation_left_allclicks_pval;       
            covariate_stats.pref_click_max_deviation_allclicks = covariate_stats.max_deviation_right_allclicks;
            covariate_stats.pref_click_max_deviation_time = covariate_stats.max_deviation_time_right_allclicks;                
            covariate_stats.pref_click_max_deviation_pval = covariate_stats.max_deviation_right_allclicks_pval;               
            covariate_stats.clicks_max_deviation_allclicks_LR_MI = MIfun(params.ilink(mean(left_allclicks(:,idx))) ,covariate_stats.max_deviation_right_allclicks ) ;                       
            covariate_stats.clicks_max_deviation_allclicks_LR_MI_pval =    2* mean(left_allclicks(:,idx) > right_allclicks(:,idx));           
        end  
        %% get max deviation of side selectivity across all clicks
        [~,peak_idx] = max(abs(mean(left_allclicks) - mean(right_allclicks))./(mean(left_allclicks)+mean(right_allclicks)));
        covariate_stats.max_deviation_time_LR_difference_allclicks = tr(peak_idx);
        covariate_stats.max_deviation_LR_difference_max_MI_allclicks = MIfun(mean(left_allclicks(:,peak_idx)),mean(right_allclicks(:,peak_idx)));
        tmp = mean(left_allclicks(:,peak_idx)<right_allclicks(:,peak_idx));
        covariate_stats.max_deviation_LR_difference_max_MI_allclicks_pval = 2*min(tmp,1-tmp);
            
    else
        sides = {'left_clicks','right_clicks'};
    end
    tmp = mean(mean(stats.wsamples.(sides{1}),2)>mean(stats.wsamples.(sides{2}),2));
    covariate_stats.clicks_average_LR_MI = MIfun(params.ilink(mean(stats.wsamples.(sides{1})(:))),params.ilink(mean(stats.wsamples.(sides{2})(:))));
    covariate_stats.clicks_average_LR_MI_pval = 2*min(tmp,1-tmp);    
    if tmp>0.5
        covariate_stats.pref_click_by_average = 'left';
    else
        covariate_stats.pref_click_by_average = 'right';
    end
    tr = buildGLM.get_tr(stats.ws.(sides{1}).tr);
    if covariate_stats.(sides{1}).max_deviation>covariate_stats.(sides{2}).max_deviation
        covariate_stats.pref_click_by_deviation = 'left';
        idx = find(tr==covariate_stats.(sides{1}).max_deviation_time);  
        covariate_stats.pref_click_max_deviation = covariate_stats.(sides{1}).max_deviation;
        covariate_stats.pref_click_max_deviation_time = covariate_stats.(sides{1}).max_deviation_time;                
        covariate_stats.pref_click_max_deviation_pval = covariate_stats.(sides{1}).max_deviation_pval;
        covariate_stats.nonpref_click_max_deviation = covariate_stats.(sides{2}).max_deviation;
        covariate_stats.nonpref_click_max_deviation_time = covariate_stats.(sides{2}).max_deviation_time;                
        covariate_stats.nonpref_click_max_deviation_pval = covariate_stats.(sides{2}).max_deviation_pval;        
        covariate_stats.clicks_max_deviation_LR_MI = MIfun(covariate_stats.(sides{1}).max_deviation ,params.ilink(stats.ws.(sides{2}).data(idx))) ;                       
        covariate_stats.clicks_max_deviation_LR_MI_pval =    2* mean(stats.wsamples.(sides{1})(:,idx) < stats.wsamples.(sides{2})(:,idx));            
    else
        covariate_stats.pref_click_by_deviation = 'right';       
        idx = find(tr==covariate_stats.(sides{2}).max_deviation_time);     
        covariate_stats.pref_click_max_deviation = covariate_stats.(sides{2}).max_deviation;
        covariate_stats.pref_click_max_deviation_time = covariate_stats.(sides{2}).max_deviation_time;                
        covariate_stats.pref_click_max_deviation_pval = covariate_stats.(sides{2}).max_deviation_pval;       
        covariate_stats.nonpref_click_max_deviation = covariate_stats.(sides{1}).max_deviation;
        covariate_stats.nonpref_click_max_deviation_time = covariate_stats.(sides{1}).max_deviation_time;                
        covariate_stats.nonpref_click_max_deviation_pval = covariate_stats.(sides{1}).max_deviation_pval;        
        covariate_stats.clicks_max_deviation_LR_MI = MIfun(params.ilink(stats.ws.(sides{1}).data(idx)),covariate_stats.(sides{2}).max_deviation ) ;         
        covariate_stats.clicks_max_deviation_LR_MI_pval =    2* mean(stats.wsamples.(sides{1})(:,idx) > stats.wsamples.(sides{2})(:,idx));            
    end
    
    %% get max deviation of side selectivity
    [~,peak_idx] = max(abs(mean(stats.wsamples.(sides{1})) - mean(stats.wsamples.(sides{2})))./(mean(stats.wsamples.(sides{1}))+mean(stats.wsamples.(sides{2}))));
    covariate_stats.max_deviation_time_LR_difference = tr(peak_idx);
    covariate_stats.max_deviation_LR_difference_max_MI = MIfun(mean(stats.wsamples.(sides{1})(:,peak_idx)),mean(stats.wsamples.(sides{2})(:,peak_idx)));
    tmp = mean(stats.wsamples.(sides{1})(:,peak_idx)<stats.wsamples.(sides{2})(:,peak_idx));
    covariate_stats.max_deviation_LR_difference_max_MI_pval = 2*min(tmp,1-tmp);    

    %% choice modulation
    covariate_stats.choice_MI_cpoke_out_prepost = MIfun(covariate_stats.cpoke_out_left.average_prepost,covariate_stats.cpoke_out_right.average_prepost);    
    for i=1:2
        tmp = mean(mean(stats.wsamples.cpoke_out_left(:,covariate_stats.cpoke_out_left.prepost_idx{i}),2)>mean(stats.wsamples.cpoke_out_right(:,covariate_stats.cpoke_out_right.prepost_idx{i}),2));
        covariate_stats.choice_MI_pval(i) = 2*min(tmp,1-tmp);
    end

    
    %% reward modulation at left side poke
    covariate_stats.reward_MI_left = MIfun(covariate_stats.spoke_left_hit.average_prepost(2),covariate_stats.spoke_left_miss.average_prepost(2));
    tmp = mean(mean(stats.wsamples.spoke_left_hit(:,covariate_stats.spoke_left_hit.prepost_idx{2}),2)>mean(stats.wsamples.spoke_left_miss(:,covariate_stats.spoke_left_miss.prepost_idx{2}),2));
    covariate_stats.reward_MI_left_pval = 2*min(tmp,1-tmp);    
    %% reward modulation at right side poke
    covariate_stats.reward_MI_right = MIfun(covariate_stats.spoke_right_hit.average_prepost(2),covariate_stats.spoke_right_miss.average_prepost(2));
    tmp = mean(mean(stats.wsamples.spoke_right_hit(:,covariate_stats.spoke_right_hit.prepost_idx{2}),2)>mean(stats.wsamples.spoke_right_miss(:,covariate_stats.spoke_right_miss.prepost_idx{2}),2));
    covariate_stats.reward_MI_right_pval = 2*min(tmp,1-tmp);  
    
end