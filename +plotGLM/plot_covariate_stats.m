function plot_covariate_stats(stats,stats_cat)

%% get rid of bad fits right off the bat
for i=1:length(stats)
    if stats(i).badly_scaled==0
        bad_fit(i)=false;
    else
        bad_fit(i)=true;
    end
end
stats = stats(~bad_fit);

stats_cat = stats_cat(~bad_fit);

%% define significance threshold
pval=0.05;

%% import structure fields to variables in the local workspace
struct2var(stats);
struct2var(stereo_click,'all',1,'_stereo');

%% define some new variables - this should all be moved into get_covariate_stats
stereo_click_modulated = max_deviation_pval_stereo<pval;
prefers_left_by_deviation = strcmp(pref_click_by_deviation,'left');
prefers_left_by_average = strcmp(pref_click_by_average,'left');
%clicks_side_modulated_by_average = average_differences_pval<pval;
%clicks_side_modulated_by_deviation = clicks_max_deviation_differences_pval<pval;
choice_modulated = choice_MI_pval<pval;
reward_modulated = reward_MI_left_pval<pval | reward_MI_right_pval<pval;
reward_left_more_modulated = (abs(reward_MI_left)>abs(reward_MI_right));
reward_MI(reward_left_more_modulated) = reward_MI_left(reward_left_more_modulated);
reward_MI(~reward_left_more_modulated) = reward_MI_right(~reward_left_more_modulated);
reward_MI_pval(reward_left_more_modulated) = reward_MI_left_pval(reward_left_more_modulated);
reward_MI_pval(~reward_left_more_modulated) = reward_MI_right_pval(~reward_left_more_modulated);

%% reward modulation figure
figure;subplot(1,2,1);
myscatter(reward_MI_left,reward_MI_right);xlabel('reward modulation, left choice trials');ylabel('reward modulation, right choice trials');
subplot(1,2,2);
significance_histogram(reward_MI,reward_MI_pval,'nbins',15,'value_label','reward modulation (max across choices)')

figure;
mean_unsigned = bootstrp(1000,@(x,y)mean(abs(x(y<pval))),reward_MI,reward_MI_pval);
mean_frac_sig = bootstrp(1000,@(y)mean(y<pval),reward_MI_pval);
data=[(mean_unsigned) (mean_frac_sig)];
errorbar([1 2],mean(data),std(data));
set(gca,'xtick',[1 2],'xticklabel',{'Mean Signed Reward MI','Fraction Modulated'},'xlim',[0.5 2.5]);

%% plot reward kernels
%get kernels
count_pos=0;
count_neg=0;
for i=1:length(stats_cat)
    if reward_MI_pval(i)<pval
        if reward_left_more_modulated(i)
            if reward_MI(i)>0
                count_pos=count_pos+1;                
                hit_kernel_pos(count_pos,:) = stats_cat{i}.ws.spoke_left_hit.data;
                miss_kernel_pos(count_pos,:) = stats_cat{i}.ws.spoke_left_miss.data;  
            else
                count_neg = count_neg+1;
                hit_kernel_neg(count_neg,:) = stats_cat{i}.ws.spoke_left_hit.data;
                miss_kernel_neg(count_neg,:) = stats_cat{i}.ws.spoke_left_miss.data;                  
            end
        else
            if reward_MI(i)>0
                count_pos=count_pos+1;                
                hit_kernel_pos(count_pos,:) = stats_cat{i}.ws.spoke_right_hit.data;
                miss_kernel_pos(count_pos,:) = stats_cat{i}.ws.spoke_right_miss.data;  
            else
                count_neg = count_neg+1;
                hit_kernel_neg(count_neg,:) = stats_cat{i}.ws.spoke_right_hit.data;
                miss_kernel_neg(count_neg,:) = stats_cat{i}.ws.spoke_right_miss.data;                  
            end          
        end
    end
end
% do plot
figure;
subplot(2,1,1);
h(1)=shadedErrorBar(tr/1000,exp(hit_kernel_pos),{@mean,@SE},'b');
hold on;h(2)=shadedErrorBar(tr/1000,exp(miss_kernel_pos),{@mean,@SE},'r');
legend([h.mainLine],{'Correct Trial','Error Trial'});
xlabel('Time (s) after reward delivery or omission');
title('Prefer Reward Delivery, n=31');
ylabel('Gain');
subplot(2,1,2);
h(1)=shadedErrorBar(tr/1000,exp(hit_kernel_neg),{@mean,@SE},'b');
hold on;h(2)=shadedErrorBar(tr/1000,exp(miss_kernel_neg),{@mean,@SE},'r');
legend([h.mainLine],{'Correct Trial','Error Trial'});
xlabel('Time (s) after reward delivery or omission');
title('Prefer Reward Omission, n=14');
ylabel('Gain');

%% choice modulation
figure;subplot(1,3,1);
choice_MI_cpoke_out_prepost_contra_ipsi = choice_MI_cpoke_out_prepost;
choice_MI_cpoke_out_prepost_contra_ipsi(strcmp(hemisphere,'left'),:) = -choice_MI_cpoke_out_prepost(strcmp(hemisphere,'left'),:);
myscatter(choice_MI_cpoke_out_prepost_contra_ipsi);xlabel('choice modulation, premovement');ylabel('choice modulation, after movement onset');
subplot(1,3,2);
significance_histogram(choice_MI_cpoke_out_prepost_contra_ipsi(:,1),choice_MI_pval(:,1),'nbins',15,'value_label','choice modulation, premovement')
subplot(1,3,3);
significance_histogram(choice_MI_cpoke_out_prepost_contra_ipsi(:,2),choice_MI_pval(:,2),'nbins',15,'value_label','choice modulation, after movement onset')

figure;
for i=1:2
    mean_unsigned = bootstrp(1000,@(x,y)mean(abs(x(y<pval))),choice_MI_cpoke_out_prepost_contra_ipsi(:,i),choice_MI_pval(:,i));
    mean_frac_sig = bootstrp(1000,@(y)mean(y<pval),choice_MI_pval(:,i));
    data=[(mean_unsigned) (mean_frac_sig)];
    subplot(1,2,1);
    errorbar(i,mean(data(:,1)),std(data(:,i)));hold on;
    xlabel('Mean Signed Choice MI');
    set(gca,'xtick',[1 2],'xticklabel',{'pre-move' 'during move'},'xlim',[0.5 2.5]);        
    subplot(1,2,2);
    errorbar(i,mean(data(:,2)),std(data(:,i)));hold on
    xlabel('Fraction Modulated');
   set(gca,'xtick',[1 2],'xticklabel',{'pre-move' 'during move'},'xlim',[0.5 2.5]);    
end

%% plot choice modulation kernels
figure;
for k=1:2
clear pref_choice_kernel nonpref_choice_kernel
count=0;
for i=1:length(stats_cat)
    if choice_MI_pval(i,k)<5
        count=count+1;                        
        if choice_MI_cpoke_out_prepost(i,k)<0
            pref_choice_kernel(count,:) = stats_cat{i}.ws.cpoke_out_right.data;
            nonpref_choice_kernel(count,:) = stats_cat{i}.ws.cpoke_out_left.data;
        else
            pref_choice_kernel(count,:) = stats_cat{i}.ws.cpoke_out_left.data;
            nonpref_choice_kernel(count,:) = stats_cat{i}.ws.cpoke_out_right.data;        
        end
    end
end
for g=1:5
    outliers{g}= abs(zscore(mean(pref_choice_kernel,2)))>6 | abs(zscore(mean(nonpref_choice_kernel,2)))>6;
    pref_choice_kernel = pref_choice_kernel(~outliers{g},:);
    nonpref_choice_kernel = nonpref_choice_kernel(~outliers{g},:);
end
% %exclude{1} = 23; %ADS exclude 23
% %exclude{2} = 24; 
%  exclude{1} = [4 8]; % TS
%  exclude{2} = exclude{1};
tr=stats_cat{1}.ws.cpoke_out_left.tr;
if k==2
    idx = tr>0;
else
    idx = tr<0;    
end
    subplot(1,2,k);
    boots = bootstrp(1000,@nanmean,exp(pref_choice_kernel(:,idx)));
    h(1)=shadedErrorBar(tr(idx)/1000,(boots),{@mean,@std},'b');
    boots = bootstrp(1000,@nanmean,exp(nonpref_choice_kernel(:,idx)));    
    hold on;h(2)=shadedErrorBar(tr(idx)/1000,(boots),{@mean,@std},'r');    
legend([h.mainLine],{'Preferred Choice Trials','Antipreferred Choice Trials'});
xlabel('Time (s), relative to movement onset');
ylabel('Gain');
end

%% time to peak of stereo and preferred click (for only significant peaks)
a=figure;
histfun = @(data,pvals,varargin)histogram(data(pvals<pval),'NumBins',40,'BinLimits',[0 300],'DisplayStyle','stairs','Normalization','cdf',varargin{:});
for i=1:3
    figure(a);
    if i==1
        h(i) = histfun(max_deviation_time_stereo,max_deviation_pval_stereo,'edgecolor',[0 0 0]);hold on;
        med=median(max_deviation_time_stereo(max_deviation_pval_stereo<pval));
        line([1 1].*med,[1 1.1]);
        text(med,1,sprintf('%g',med));   
        boots = bootstrp(1000,@nanmedian,max_deviation_time_stereo(max_deviation_pval_stereo<pval));
        
    elseif i==2 % these two could be substituted with the allclick version if necessary 
        h(i) = histfun(pref_click_max_deviation_time,pref_click_max_deviation_pval,'edgecolor',[0.4 0.4 0.4]);
        med=median(pref_click_max_deviation_time(pref_click_max_deviation_pval<pval));
        line([1 1].*med,[1 1.1]);
        text(med,1,sprintf('%g',med)); 
        boots = bootstrp(1000,@nanmean,pref_click_max_deviation_time(pref_click_max_deviation_pval<pval));
        
    else % these two could be substituted with the allclick version if necessary
        h(i) = histfun(nonpref_click_max_deviation_time,nonpref_click_max_deviation_pval,'edgecolor',[0.8 0.8 0.8]);
        line([1 1]*median(nonpref_click_max_deviation_time(nonpref_click_max_deviation_pval<pval)),[1 1.1],'color',[0.8 0.8 0.8]);        
        boots = bootstrp(1000,@nanmedian,nonpref_click_max_deviation_time(nonpref_click_max_deviation_pval<pval));
        
    end
    figure(123);
    if i==1
        clf
    end
    errorbar(i,mean(boots),std(boots));hold on;title('time to peak amplitude');
end
legend(h,{'time to peak stereo','time to peak pref click least adapted third','time to peak nonpref click, least adapted third'});



%% time to peak of side selectivity in clicks
figure;
histfun = @(data,pvals)histogram(data(pvals<pval),'NumBins',20,'BinLimits',[0 300],'DisplayStyle','stairs','Normalization','cdf');
for i=1:3
    if i==1
        h(i) = histfun(max_deviation_time_LR_difference,max_deviation_LR_difference_max_MI_pval);hold on;
        med=median(max_deviation_time_LR_difference(max_deviation_LR_difference_max_MI_pval<pval));
        line([1 1].*med,[1 1.1]);
        text(med,1,sprintf('%g',med));
        boots = bootstrp(1000,@nanmean,max_deviation_time_LR_difference(max_deviation_LR_difference_max_MI_pval<pval));
    elseif i==2 %  allclick version 
        h(i) = histfun(max_deviation_time_LR_difference_allclicks,max_deviation_LR_difference_max_MI_allclicks_pval);hold on;
    end
    legend(h,{'Time to peak side selective, least adapted third of clicks','Time to peak side selectivity, average across all clicks'});
end
figure;
    if i==1
        clf
    end
    errorbar(i,mean(boots),std(boots));hold on;title('median of time to peak side selectivity, least adapted third');
    
figure;
significance_histogram(clicks_allclicks_average_LR_MI,clicks_allclicks_average_LR_MI_pval,'nbins',150);



% plots
% adaptation




%% get side click deviation stats
clear max_deviation_to_plot
for i=1:3
    for c=1:length(pref_click_by_average)
        if strcmp(pref_click_by_average(c),'left')
            eval(['curr_clicks = left_clicks',num2str(i)]);
        else
            eval(['curr_clicks = right_clicks',num2str(i)]);
        end
        max_deviation_to_plot(c,i) = curr_clicks(c).max_deviation;        
    end
end



%% plot adaptation

max_deviation_stereo_to_plot = max_deviation_stereo;


max_deviation_is_negative = max_deviation_to_plot(:,3)<1;
max_deviation_to_plot(max_deviation_is_negative,:) = 1./max_deviation_to_plot(max_deviation_is_negative,:);
max_deviation_stereo_to_plot(max_deviation_stereo_to_plot<1) = 1./max_deviation_stereo_to_plot(max_deviation_stereo_to_plot<1);
max_deviation_is_outlier = abs(zscore(max_deviation_to_plot))>5 | max_deviation_to_plot>100;
max_deviation_to_plot(max_deviation_is_outlier)=NaN;
max_deviation_stereo_is_outlier = max_deviation_stereo_to_plot>100 | abs(zscore(max_deviation_stereo_to_plot(:)))>5;
max_deviation_stereo_to_plot(max_deviation_stereo_is_outlier)=NaN;

stereo_modulated = max_deviation_pval_stereo<pval;
max_deviation_to_plot(~stereo_modulated,:)=NaN;
max_deviation_stereo_to_plot(~stereo_modulated)=NaN;

figure;errorbar([1 2 3 4],mean([max_deviation_stereo_to_plot(:) fliplr(max_deviation_to_plot)]),std([max_deviation_stereo_to_plot(:) fliplr(max_deviation_to_plot)])./sqrt(size(max_deviation_to_plot,1)));
%figure;plot([1 2 3 4],([max_deviation_stereo_to_plot(:) max_deviation_to_plot])');
set(gca,'xtick',1:4,'xticklabel',{'Stim Start (stereo click)','most adapted clicks (pref side)','','least adapted clicks (pref side)'});
ylabel('Spike Rate Fold Change at Peak Response');
xlabel('Click Covariate');

%% plot click kernels
count=0;
clear pref_click_kernel nonpref_click_kernel
for i=1:length(stats_cat)
    count=count+1;                            
    if strcmp(pref_click_by_deviation{i},'left') %choice_MI_cpoke_out_prepost(i,1)>0 % clicks_average_LR_MI>0 
        for k=1:3
            pref_click_kernel(count,k,:) = exp(stats_cat{i}.ws.(['left_clicks',num2str(k)]).data);
            nonpref_click_kernel(count,k,:) = exp(stats_cat{i}.ws.(['right_clicks',num2str(k)]).data);            
        end
    else
        for k=1:3
            pref_click_kernel(count,k,:) = exp(stats_cat{i}.ws.(['right_clicks',num2str(k)]).data);
            nonpref_click_kernel(count,k,:) = exp(stats_cat{i}.ws.(['left_clicks',num2str(k)]).data);      
        end        
    end
end
%pref_click_kernel(max_deviation_is_negative,:,:) = 1./pref_click_kernel(max_deviation_is_negative,:,:);
%nonpref_click_kernel(max_deviation_is_negative,:,:) = 1./nonpref_click_kernel(max_deviation_is_negative,:,:);
for i=1:2

bad_idx = abs(zscore(log(mean(mean(pref_click_kernel,3),2))))>1;
pref_click_kernel(bad_idx,:,:)=NaN;
bad_idx = abs(zscore(log(mean(mean(nonpref_click_kernel,3),2))))>1;
nonpref_click_kernel(bad_idx,:,:)=NaN;
end
tr=stats_cat{1}.ws.left_clicks1.tr;
for i=1:1
figure;
if i==2
    idx = max_deviation_is_negative & clicks_average_LR_MI_pval<pval;% pref_click_max_deviation_pval<pval
else
    idx = ~max_deviation_is_negative ;%& clicks_average_LR_MI_pval<pval; %& pref_click_max_deviation_pval<pval
end

for g=1:3
subplot(1,3,g);
%diffboots = bootstrp(1000,@nanmean,squeeze(pref_click_kernel(idx,g,:)-nonpref_click_kernel(idx,g,:)));
%h(1)=shadedErrorBar(tr/1000,diffboots,{@mean,@std},'b');
h(1) = shadedErrorBar(tr/1000,(squeeze(pref_click_kernel(idx,g,:))),{@mean,@SE},'b');
hold on;h(2)=shadedErrorBar(tr/1000,(squeeze(nonpref_click_kernel(idx,g,:))),{@mean,@SE},'r');
legend([h.mainLine],{'Pref Click','Nonpref Click'});
xlabel('Time (s) after click');
ylabel('Gain');
matchylim(gcf);
end
end

%% plot click side modulation (need to fix to hold them to the same tr) did that, but need a better measure of significance. really need that fucking covariance across covariates.
left_right_max_deviation_to_plot = left_right_max_deviation( stereo_click_modulated(:) ,3);
left_right_max_deviation_to_plot(max_deviation_is_negative) = -left_right_max_deviation_to_plot(max_deviation_is_negative);
%left_right_max_deviation_to_plot = left_right_max_deviation_to_plot(max_deviation_pval(: ,3)<1);
pvals = left_right_max_deviation_pval( stereo_click_modulated(:)); 
%figure;h=histogram(left_right_max_deviation_to_plot);h.NumBins=20;
figure;CPhistogram(left_right_max_deviation_to_plot,pvals)


%% plot stereo click
clear stereo_click_kernel
for i=1:length(stats_cat)
    stereo_click_kernel(i,:) = exp(stats_cat{i}.ws.stereo_click.data);
end
%pref_click_kernel(max_deviation_is_negative,:,:) = 1./pref_click_kernel(max_deviation_is_negative,:,:);
%nonpref_click_kernel(max_deviation_is_negative,:,:) = 1./nonpref_click_kernel(max_deviation_is_negative,:,:);
for i=1:1

    
    bad_idx=any(abs(zscore(stereo_click_kernel))>5,2);
    
%bad_idx = abs(zscore((m(stereo_click_kernel,2))))>4;
stereo_click_kernel(bad_idx,:)=NaN;
end
tr=stats_cat{1}.ws.stereo_click.tr;
for i=1:2
figure;
if i==2
    idx1 = max_deviation_stereo<1 & stereo_modulated;% & clicks_average_LR_MI_pval<pval;% pref_click_max_deviation_pval<pval
else
    idx2 =  (max_deviation_stereo>1 ) & stereo_modulated; %& clicks_average_LR_MI_pval<pval; %& pref_click_max_deviation_pval<pval
end



end
h = shadedErrorBar(tr/1000,cat(1,(squeeze(stereo_click_kernel(idx2,:))),(squeeze(stereo_click_kernel(idx2,:)))),{@mean,@SE},'b');
legend(h.mainLine,'Stereo Click');
xlabel('Time (s) after click');
ylabel('Gain');
%diffboots = bootstrp(1000,@nanmean,squeeze(pref_click_kernel(idx,g,:)-nonpref_click_kernel(idx,g,:)));
%h(1)=shadedErrorBar(tr/1000,diffboots,{@mean,@std},'b');


%% clicks modulation and selectivity statistics
figure;
boots = bootstrp(1000,@nanmean,stereo_click_modulated<pval);
boots2 = bootstrp(1000,@nanmean,clicks_max_deviation_LR_MI_pval<pval);
figure;errorbar([1 2],mean([boots boots2]),std([boots boots2]));
set(gca,'xtick',[1 2],'xticklabel',{'fraction click modulated','fraction side selective'});

