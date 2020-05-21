function matchylim(figs,varargin)
    % puts all heatmaps in the input list of figure handles into a common colormap,
    % bounded by the most extreme values across the entire set of heatmaps
    j=1;
    method='';
    ylim=[];
    while j<=length(varargin)
        if strncmpi(varargin{j},'usedata',4)
            method='tight';
        elseif strncmpi(varargin{j},'ylim',4)
            j=j+1;
            ylim=varargin{j};
        elseif strncmpi(varargin{j},'tight',5)
            method='tight';
        end
        j=j+1;
    end  
    if isempty(ylim)
        switch method
            case 'tight'
                yl=[];
                for f=1:length(figs)
                   hs=get(figs(f),'children');
                   data = cat(1,get(hs,'YData'));
                   yl=[yl ; data(:)];
                   if iscell(yl)
                    yl=[yl{:}];
                   end
                end               
            otherwise
                yl=[];
                for f=1:length(figs)
                    kids=get(figs(f),'children');
                    for k=1:length(kids)
                        if isprop(kids(k),'ylim')
                            yl = [yl get(kids(k),'ylim')];                
                        end
                    end
                end
        end
        ylim=minmax(yl);
    end
    for f=1:length(figs)
        figure(figs(f));
        kids=get(gcf,'children');
        for k=1:length(kids)
            if isprop(kids(k),'ylim')
                set(kids(k),'ylim',ylim);
            end
        end
    end
end