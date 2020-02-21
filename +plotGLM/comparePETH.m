function [h,r_square] = comparePETH(stats,varargin)

    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('group_by','gamma');
    p.addParameter('ngroups',1);
    p.addParameter('align_to','stereo_click');
    p.addParameter('onlyCorrect',false);
    p.addParameter('subplots',{[2 1 1],[2 1 2]});
    p.parse(varargin{:});
    params=p.Results;
    subplot = @(m,n,p)subtightplot(m,n,p);

    ntrials = length(stats.dspec.expt.trial);        
    if params.ngroups>1
        if isfield(stats.dspec.expt.trial,params.group_by)
            group_vals = arrayfun(@(x)x.(params.group_by),stats.dspec.expt.trial);
        elseif strcmp(params.group_by,'click_diff')
            fields = fieldnames(stats.dspec.expt.trial(1));
            left_fields = find(strncmp(fields,'left_click',8));
            right_fields = find(strncmp(fields,'right_click',8));
            for i=1:length(stats.dspec.expt.trial)
                n_left(i) = 0;
                for k=1:length(left_fields)
                    n_left(i) = n_left(i) + length(stats.dspec.expt.trial(i).(fields{left_fields(k)}));
                end
                n_right(i) = 0;
                for k=1:length(left_fields)
                    n_right(i) = n_right(i) + length(stats.dspec.expt.trial(i).(fields{right_fields(k)}));
                end            
            end
            group_vals = n_left - n_right;
        end
        group_idx = discretize(group_vals,params.ngroups);
    else
        group_idx = ones(1,ntrials);
    end
    if params.onlyCorrect
        error = ~[stats.dspec.expt.trial.ishit];
        group_idx(error)=0;
    end
    colors = redblue(params.ngroups);
    whites = find(sum(colors,2)==3);
    colors(whites,:) = ones(length(whites),3)*0.9;
    yl=[1 1];
    for i=1:params.ngroups
        % plot spikes
        tmp=num2cell(params.subplots{1});
        subplot(tmp{:});
        [h{i,1},filtered_spikes{i}] = plotGLM.plotPETH(stats.dspec.expt,['sptrain',num2str(stats.cellno)],find(group_idx==i),'align_to',params.align_to,'color',colors(i,:),'sd',50,varargin{:});hold on
        mean_spikes(i,:) = nanmean(filtered_spikes{i},1);
        tmp=num2cell(params.subplots{2});
        drawnow;ylnow=get(gca,'ylim');
        yl(1) = min(yl(1),ylnow(1));
        yl(2) = max(yl(2),ylnow(2));
        set(gca,'ylim',yl);
        subplot(tmp{:});
        [h{i,2},filtered_spikes_pred{i}] = plotGLM.plotPETH(stats.dspec.expt,stats.Yhat,find(group_idx==i),'align_to',params.align_to,'color',colors(i,:),'sd',40,varargin{:});hold on
        mean_spikes_pred(i,:) = nanmean(filtered_spikes_pred{i},1);   
        drawnow;ylnow=get(gca,'ylim');
        yl(1) = min(yl(1),ylnow(1));
        yl(2) = max(yl(2),ylnow(2));
        set(gca,'ylim',yl);
    end
    filtered_spikes_pred = filtered_spikes_pred{:};
    filtered_spikes = filtered_spikes{:};
    r_square(1) = rsquare(filtered_spikes(:),filtered_spikes_pred(:));
    r_square(2) = rsquare(mean_spikes(:),mean_spikes_pred(:));
end