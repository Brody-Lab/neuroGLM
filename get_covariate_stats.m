function [covariate_stats,stats]=get_covariate_stats(stats,varargin)
    p=inputParser;
    p.addParameter('ilink',@(x)exp(x),@(x)validateattributes(x,{'function_handle'},{''}));
    p.parse(varargin{:});
    params=p.Results;
    if ismember('ilink',p.UsingDefaults) && isfield(stats.params,'distribution') && strcmp(stats.params.distribution,'poisson')
        params.ilink = @(x)exp(x);
    end  
    DC = params.ilink(0);
    covars = {stats.dspec.covar.label};
    fprintf('get_covariate_stats: cell %g\n',stats.cellno);
    if ~isfield(stats.dspec.expt.param,'nClickBins')
        stats.dspec.expt.param.nClickBins = sum(contains({stats.dspec.covar.label},'left_click'));    
    end
    %% check that cov is in struct
    if ~isfield(stats.wvars.(covars{1}),'cov')
        fprintf('get_covariate_stats: adding covariance matrix at solution grouped by covariate ... ');a=tic;
        dm = buildGLM.compileSparseDesignMatrix(stats.dspec, 1:length(stats.dspec.expt.trial));                
        dm = buildGLM.removeConstantCols(dm);      
        [stats.ws,stats.wvars] = buildGLM.combineWeights(buildGLM.addBiasColumn(dm), stats.beta , stats.covb);
        fprintf('took %s.\n',timestr(toc(a)));
    end    
    %% get some basic stats from all covariates
    for i=1:length(covars)
        if contains(covars{i},'left_click') || contains(covars{i},'right_click')
            continue
        end
        data = stats.ws.(covars{i}).data;
        tr = stats.ws.(covars{i}).tr;
        cov = stats.wvars.(covars{i}).cov;  
        vars = stats.wvars.(covars{i}).data;        
        covariate_stats.(covars{i}).average(i) = mean(params.ilink(data));
        covariate_stats.(covars{i}).average_pval(i) = empirical_p(DC,mean(params.ilink(mvnrnd(data,cov,1e3)),2)); 
        if any(tr<0)
            ops = {@gt,@lt};
            for k=1:2
                idx = ops{k}(tr,0);
                data = stats.ws.(covars{i}).data(idx);
                covariate_stats.(covars{i}).average_prepost(k) = mean(params.ilink(data));
                covariate_stats.(covars{i}).average_pval_prepost(k) = empirical_p(DC,mean(params.ilink(mvnrnd(data,cov(idx,idx),1e3)),2));                     
            end
        end
        if strcmp(covars{i},'stereo_click')
            [~,peak_idx] = max(abs(params.ilink(data)-DC));
            covariate_stats.(covars{i}).max_deviation_time = stats.ws.(covars{i}).tr(peak_idx);
            covariate_stats.(covars{i}).max_deviation = params.ilink(data(peak_idx));  
            covariate_stats.(covars{i}).max_deviation_pval = 2*min(normcdf(0,[-1 1].*data(peak_idx),sqrt(vars(peak_idx))));
        end            
    end             

    %% left and right clicks    
    for i=1:stats.dspec.expt.param.nClickBins
       for side_idx=1:2
            sides = {'left_clicks','right_clicks'};           
            label=[sides{side_idx},num2str(i)];
            data = stats.ws.(label).data;
            vars = stats.wvars.(label).data;
            cov = stats.wvars.(label).cov;                   
            [~,peak_idx] = max(abs(params.ilink(data)-DC));
            covariate_stats.(sides{side_idx}).max_deviation_time(i) = stats.ws.(label).tr(peak_idx);
            covariate_stats.(sides{side_idx}).max_deviation(i) = params.ilink(data(peak_idx));  
            covariate_stats.(sides{side_idx}).max_deviation_pval(i) = 2*min(normcdf(0,[-1 1].*data(peak_idx),sqrt(vars(peak_idx))));
       end
    end
    MIfun = @(x,y)(x-y)./(x+y);    
    tr = stats.ws.left_clicks1.tr;
    if covariate_stats.(sides{1}).max_deviation>covariate_stats.(sides{2}).max_deviation
        covariate_stats.pref_click = 'left';
        idx = find(tr==covariate_stats.(sides{1}).max_deviation_time(3));        
        covariate_stats.left_right_max_deviation = MIfun(covariate_stats.(sides{1}).max_deviation ,params.ilink(stats.ws.right_clicks3.data(idx))) ;                       
        covariate_stats.left_right_max_deviation_pval = 2*min(normcdf(stats.ws.right_clicks3.data(idx),[-1 1].*stats.ws.left_clicks3.data(idx),sqrt(stats.wvars.left_clicks3.data(idx))));
    else
        covariate_stats.pref_click = 'right';       
        idx = find(tr==covariate_stats.(sides{2}).max_deviation_time(3));        
        covariate_stats.left_right_max_deviation = MIfun(params.ilink(stats.ws.left_clicks3.data(idx)),covariate_stats.(sides{2}).max_deviation ) ;         
        covariate_stats.left_right_max_deviation_pval = 2*min(normcdf(stats.ws.left_clicks3.data(idx),[-1 1].*stats.ws.right_clicks3.data(idx),sqrt(stats.wvars.right_clicks3.data(idx))));
        
    end

    %% choice modulation
    covariate_stats.choice_MI_cpoke_out_prepost = MIfun(covariate_stats.cpoke_out_left.average_prepost,covariate_stats.cpoke_out_right.average_prepost);
    covariate_stats.choice_MI_pval = min(covariate_stats.cpoke_out_left.average_pval_prepost,covariate_stats.cpoke_out_right.average_pval_prepost);

    
    %% reward modulation at prefered side poke
    if covariate_stats.choice_MI_cpoke_out_prepost(2)>0
        covariate_stats.reward_MI_pref_side = MIfun(covariate_stats.spoke_left_hit.average_prepost(2),covariate_stats.spoke_left_miss.average_prepost(2));
        covariate_stats.reward_pval = min(covariate_stats.spoke_left_hit.average_pval_prepost(2),covariate_stats.spoke_left_miss.average_pval_prepost(2));        
    else
        covariate_stats.reward_MI_pref_side = MIfun(covariate_stats.spoke_right_hit.average_prepost(2),covariate_stats.spoke_right_miss.average_prepost(2));        
        covariate_stats.reward_pval = min(covariate_stats.spoke_right_hit.average_pval_prepost(2),covariate_stats.spoke_right_miss.average_pval_prepost(2));
    end
    
end