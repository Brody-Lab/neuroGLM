function [mean_values,pval] = significance_histogram(values,varargin)
    p=inputParser;
    p.KeepUnmatched=true;
    dividingLine = 0;
    p.addRequired('values',@(x)validateattributes(x,{'numeric'},{}));
    p.addRequired('pvals',@(x)validateattributes(x,{'numeric'},{'nonnegative'}));
    p.addParameter('nbins',15,@(x)validateattributes(x,{'numeric'},{'positive','integer','scalar'}));
    p.addParameter('bin_limits',[NaN NaN]);
    p.addParameter('value_label','');
    p.parse(values,varargin{:});
    [values,pvals,value_label]=struct2var(p.Results,{'values','pvals','value_label'});
    hold on
    if ~ismember('pvals',p.UsingDefaults)
        if any(size(pvals)~=size(values)) 
            error('values and pvals must be the same size.');
        end
        sig=histogram(values(pvals<0.05));
        maxOffset=max(abs([min(values) max(values)]-dividingLine));
        if any(isnan(p.Results.bin_limits))
           sig.BinLimits=dividingLine + [-maxOffset maxOffset];
        else
           sig.BinLimits = p.Results.bin_limits; 
        end
        sig.NumBins=nearestEven(p.Results.nbins);
        nonsig=histogram(values(pvals>0.05));
        if any(isnan(p.Results.bin_limits))
           nonsig.BinLimits=dividingLine + [-maxOffset maxOffset];
        else
           nonsig.BinLimits = p.Results.bin_limits; 
        end
        nonsig.NumBins=nearestEven(p.Results.nbins);
        xval=sig.BinEdges([1:end-1])+sig.BinWidth/2;
        sigval=sig.Values;
        nonsigval=nonsig.Values;
        xl=sig.BinLimits;
        cla;
        bar_handle=bar(xval',[sigval;nonsigval]','stacked');
        bar_handle(2).FaceColor=[1 1 1];
        bar_handle(1).FaceColor=[0 0 0];
    else
        h = histogram(values);
        h.NumBins=p.Results.nbins;
        h.BinLimits=minmax(values);
    end
    mean_values=nanmean(values);
    bootmean=bootstrp(1000,@nanmean,values);
    pval=empirical_p(dividingLine,bootmean);
    if pval<0.05
        sigstring='*';
        if pval<0.01
            sigstring='**';
        end
    else
        sigstring='n.s.';
    end
    line([dividingLine dividingLine],get(gca,'ylim'),'LineStyle',':','color',[1 1 1]/2);
    yl=get(gca,'ylim');
    set(gca,'xlim',xl);
    patch([-abs(diff(xl))/20 0 abs(diff(xl))/20]+mean_values,[0 -abs(diff(yl))/20 0]+yl(2)-abs(diff(yl))/20,[1 1 1]/2);
    text(mean_values,yl(2)-abs(diff(yl))/20,sigstring,'HorizontalAlignment','center','FontSize',15,'VerticalAlignment','baseline');
    title(['mean = ',num2str(round(mean_values,3)),', sig only = ',num2str(round(mean(values(pvals<0.05)),3))]);
    box off;
    set(gca,'xgrid','on');
    ylabel('Number');    
    xlabel(value_label);
end