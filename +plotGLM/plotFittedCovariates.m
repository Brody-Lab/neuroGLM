function plotFittedCovariates(stats,varargin)
    p=inputParser;
    p.addParameter('matchylim',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('ilink',@(x)exp(x),@(x)validateattributes(x,{'function_handle'},{}));
    p.addParameter('showOffClicks',false);
    p.parse(varargin{:});
    params=p.Results;
    if ismember('link',p.UsingDefaults) && isfield(stats.params,'distribution') && strcmp(stats.params.distribution,'poisson')
        params.ilink = @(x)exp(x);
    end
    covars = {stats.dspec.covar.label};
    matches={'left','right'};
    groups = group_covars(covars,matches);
    color.left='b';
    color.right='r';
    n_subplot_columns = ceil(sqrt(length(groups)));
    mmx=[];
    group_count=0;
    for kCov = 1:length(groups)
        count=0;
        if contains(groups{kCov}{1},'click') && params.showOffClicks
            continue
        end
        group_count=group_count+1;
        for j=1:length(groups{kCov})
            count=count+1;
            label = groups{kCov}{j};
            if contains(label,'click') && params.showOffClicks
                continue
            end
            ts=stats.ws.(label).tr;          
            subplot(n_subplot_columns,n_subplot_columns,group_count);hold on;
            %subplot(2,3,group_count);
            if length(groups{kCov})==1
                shadedErrorBar(ts, params.ilink(stats.ws.(label).data), params.ilink(stats.ws.(label).data+sqrt(stats.wvars.(label).data(:)'))-params.ilink(stats.ws.(label).data));hold on
            else
                h(count) =  shadedErrorBar(ts, params.ilink(stats.ws.(label).data), params.ilink(stats.ws.(label).data+sqrt(stats.wvars.(label).data(:)'))-params.ilink(stats.ws.(label).data),color.(matches{j}));  hold on
                if count>1
                    h(count).patch.FaceAlpha=0.55555;
                    if strcmp(label,'left_clicks') || strcmp(label,'right_clicks')
                        legend([h.mainLine],{'left clicks','right clicks'})
                    else
                        legend([h.mainLine],{'left choice','right choice'})
                    end
                end
            end
            mmx = max([mmx(:);params.ilink(stats.ws.(label).data(:))]);
            title_str = strrep(groups{kCov}{j},'_',' ');
            title_str = strrep(title_str,matches{j},'');
            title(title_str);
            set(gca,'xlim',[min(ts) max(ts)]);
            xlabel('time (s)');
            ylabel('Gain');
        end
    end
    fig=gcf;
    if params.showOffClicks
        figure(1253);
        clf;
        subplot(1,3,1);
        shadedErrorBar(stats.ws.stereo_click.tr, params.ilink(stats.ws.stereo_click.data), params.ilink(stats.ws.stereo_click.data+sqrt(stats.wvars.stereo_click.data(:)'))-params.ilink(stats.ws.stereo_click.data));
        title('Stereo Click');
        xlabel('Time (s)');
        ylabel('Gain');
        subplot(1,3,2);
        colors=jet(6);
        for i=1:3 % right
            label = ['left_clicks',num2str(i)];
            h(i)=shadedErrorBar(stats.ws.(label).tr, params.ilink(stats.ws.(label).data), params.ilink(stats.ws.(label).data+sqrt(stats.wvars.(label).data(:)'))-params.ilink(stats.ws.(label).data),{'color',colors(3-i+1,:)}); hold on;
            h(i).patch.FaceAlpha=0;
            title('Left Clicks');
            xlabel('Time (s)');
            ylabel('Gain');
        end
        legend([h.mainLine],{'most adapted third','middle adapted third','least adapted third'});
        subplot(1,3,3);
        for i=1:3 % right
            label = ['right_clicks',num2str(i)];
            h(i)=shadedErrorBar(stats.ws.(label).tr, params.ilink(stats.ws.(label).data), params.ilink(stats.ws.(label).data+sqrt(stats.wvars.(label).data(:)'))-params.ilink(stats.ws.(label).data),{'color',colors(i+3,:)});   hold on;
            h(i).patch.FaceAlpha=0;
            title('Right Clicks');
            xlabel('Time (s)');
            ylabel('Gain');
        end
        legend([h.mainLine],{'most adapted third','middle adapted third','least adapted third'});
    end
    matchylim(gcf);
    figure(fig);
    if params.matchylim
        try
            if params.ilink(1)~=1
                matchylim(gcf,'ylim',[0 1]*max(abs(mmx))*1.1);
            else
                matchylim(gcf,'ylim',[-1 1]*max(abs(mmx))*1.1);
            end
        end
    end
    set(gcf,'position',[ 1000.3          221       1101.3         1110]);
end

function groups = group_covars(covars,matches)
    count=0;
    groups={};
    for i=1:length(covars)
        matched=false(1,length(matches));
        if ismember(covars{i},[groups{:}])
            continue
        end
        for k=1:length(matches)
            if contains(covars{i},matches{k})
                matched(k)=true;
                count=count+1;
                groups{count}=covars(i);
                other_match_patterns = setdiff(matches,matches{k});
                for j=1:length(other_match_patterns)
                    other_matches = ismember(covars,strrep(covars{i},matches{k},other_match_patterns{j}));
                    groups{count} = unique([groups{count} covars(other_matches)]);
                end
                break
            end
            if ~any(matched) && k==length(matches)
                count=count+1;
                groups{count}=covars(i);
            end
        end
    end
end
