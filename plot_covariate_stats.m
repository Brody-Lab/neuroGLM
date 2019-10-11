function plot_covariate_stats(stats)

pval=0.05;


struct2var([stats.covariate_stats]);
struct2var(stereo_click,'all',1,'_stereo');
struct2var(left_clicks,'all',1,'_left_clicks');
struct2var(right_clicks,'all',1,'_right_clicks');

click_modulated = max_deviation_pval_stereo<pval;
prefers_left = strcmp(pref_click,'left');


%% plot stereo clicks
figure;
count=0;
for i=1:length(stats)
    if ~isempty(stats(i).ws) && stats(i).badly_scaled==0
        bad_fit(i)=false;
        if click_modulated(i)
            count=count+1;
            data(count,:) = stats(i).ws.stereo_click.data;
            vars = stats(i).wvars.stereo_click.data;
            xs = stats(i).ws.stereo_click.tr;
            if prefers_left(i)
                max_deviation(count,:) = max_deviation_left_clicks(count,:);
                max_deviation_pval(count,:) = max_deviation_pval_left_clicks(count,:);
                max_deviation_time(count) = max_deviation_time_left_clicks(count,3);
            else
                max_deviation(count,:) = max_deviation_right_clicks(count,:);
                max_deviation_pval(count,:) = max_deviation_pval_right_clicks(count,:);   
                max_deviation_time(count) = max_deviation_time_right_clicks(count,3);                
            end
            if max_deviation(count,end)<1
                max_deviation_is_negative(count)=true;
            else
                max_deviation_is_negative(count)=false;
            end
        end
    else
        bad_fit(i)=true;
    end
end

idx = ~bad_fit(:) & click_modulated(:);
max_deviation_time_stereo = max_deviation_time_stereo(idx);
data_pos = exp(data(max_deviation_stereo(idx)>1 & max_deviation_time_stereo>40 ,:));
data_pos2 = exp(data(max_deviation_stereo(idx)>1 & max_deviation_time_stereo<40 ,:));

data_neg = exp(data(max_deviation_stereo(idx)<1 ,:));
%data_pos = bsxfun(@times,data_pos,1./max(data_pos,[],2));
%data_neg = bsxfun(@times,data_neg,1./min(data_neg,[],2));

figure;plot(xs,data_pos','r')
hold on;plot(xs,data_pos2','g')
hold on;plot(xs,data_neg','b')
figure;shadedErrorBar(xs,data_pos,{@mean,@SE},'r')
hold on;shadedErrorBar(xs,data_neg,{@mean,@SE},'b')
hold on;shadedErrorBar(xs,data_pos2,{@mean,@SE},'g')
figure;
%% plot adaptation
max_deviation_to_plot = max_deviation;
for i=1:5
    max_deviation_to_plot(abs(zscore(max_deviation_to_plot))>5)=NaN;
end
max_deviation_to_plot(max_deviation_is_negative) = 1./max_deviation_to_plot(max_deviation_is_negative);
max_deviation_stereo_to_plot = max_deviation_stereo(~bad_fit(:) & click_modulated(:));
max_deviation_stereo_to_plot(max_deviation_stereo_to_plot<1) = 1./max_deviation_stereo_to_plot(max_deviation_stereo_to_plot<1);
[bads,~] = find(max_deviation_to_plot>50);
bads=unique(bads);
 max_deviation_to_plot(bads,:) = [];
 max_deviation_stereo_to_plot(bads)=[];
figure;errorbar([1 2 3 4],mean([max_deviation_stereo_to_plot(:) max_deviation_to_plot]),std([max_deviation_stereo_to_plot(:) max_deviation_to_plot])./sqrt(size(max_deviation_to_plot,1)));
%figure;plot([1 2 3 4],([max_deviation_stereo_to_plot(:) max_deviation_to_plot])');
set(gca,'xtick',1:4,'xticklabel',{'Stim Start (stereo click)','most adapted clicks (pref side)','','least adapted clicks (pref side)'});
ylabel('Spike Rate Fold Change at Peak Response');
xlabel('Click Covariate');

%% plot click side modulation (need to fix to hold them to the same tr) did that, but need a better measure of significance. really need that fucking covariance across covariates.
left_right_max_deviation_to_plot = left_right_max_deviation(~bad_fit(:) & click_modulated(:) ,3);
left_right_max_deviation_to_plot(max_deviation_is_negative) = -left_right_max_deviation_to_plot(max_deviation_is_negative);
%left_right_max_deviation_to_plot = left_right_max_deviation_to_plot(max_deviation_pval(: ,3)<1);
pvals = left_right_max_deviation_pval(~bad_fit(:) & click_modulated(:)); 
%figure;h=histogram(left_right_max_deviation_to_plot);h.NumBins=20;
figure;CPhistogram(left_right_max_deviation_to_plot,pvals)
end