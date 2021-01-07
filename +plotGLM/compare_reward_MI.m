function compare_reward_MI(stats1,stats2,params)
    %% parse and validate inputs
    arguments
        stats1 struct
        stats2 struct
        params.color1 (1,3) {mustBeNumeric,mustBeLessThanOrEqual(params.color1,1),...
            mustBeGreaterThanOrEqual(params.color1,0)} = [1 1 1]/2
        params.color2 (1,3) {mustBeNumeric,mustBeLessThanOrEqual(params.color2,1),...
            mustBeGreaterThanOrEqual(params.color2,0)} = spectrumRGB(473)
        params.nbins (1,1) {mustBeNumeric,mustBePositive} = 30      
        params.log_scale (1,1) logical = false
    end
    P = get_parameters;
    inputnames = {strrep(inputname(1),'_','\_') strrep(inputname(2),'_','\_')};
    
    %% compare reward selectivity, pre move
    reward_MI_left{1} = cat(1,stats1.reward_MI_left);
    reward_MI_right{1} = cat(1,stats1.reward_MI_right);
    reward_MI_left{2} = cat(1,stats2.reward_MI_left);
    reward_MI_right{2} = cat(1,stats2.reward_MI_right);
    reward_MI_left_pval{1} = cat(1,stats1.reward_MI_left_pval);
    reward_MI_right_pval{1} = cat(1,stats1.reward_MI_right_pval);
    reward_MI_left_pval{2} = cat(1,stats2.reward_MI_left_pval);
    reward_MI_right_pval{2} = cat(1,stats2.reward_MI_right_pval);    
    for i=1:2
        reward_MI_avg{i} = (reward_MI_left{i} + reward_MI_right{i})/2;
        left_bigger = abs(reward_MI_left{i}) > abs(reward_MI_right{i});
        reward_MI_max{i}(left_bigger) = reward_MI_left{i}(left_bigger);
        reward_MI_max{i}(~left_bigger) = reward_MI_right{i}(~left_bigger);        
    end
    
    
    for i=1:2
        left_sig{i} = reward_MI_left_pval{i}<0.05;
        right_sig{i} = reward_MI_right_pval{i}<0.05;        
        reward_MI_max{i}(left_sig{i} & ~right_sig{i}) = reward_MI_left{i}(left_sig{i} & ~right_sig{i});
        reward_MI_max{i}(~left_sig{i} & right_sig{i}) = reward_MI_left{i}(~left_sig{i} & right_sig{i});    
        reward_MI_max_pval{i}(~left_sig{i} & ~right_sig{i}) = 1;
        reward_MI_max_pval{i}(left_sig{i} | right_sig{i}) = 0; 
        subplot(1,3,i);
        significance_histogram(reward_MI_avg{i}(:),reward_MI_max_pval{i}(:));        
        set(gca,P.axes_properties{:},'xlim',[-1 1]);
        xlabel('Reward MI');
        ylabel('No. Cells');
        title('');
        if i==2
            title('Tagged (Putative D2)');
        else
            title('Untagged');
        end
    end
        

    %% compare reward selectivity for both sides
%     for i=1:2
%         subplot(1,5,i+2);
%         myscatter(reward_MI_left{i},reward_MI_right{i});
%         set(gca,P.axes_properties{:});
%         xlabel('Reward MI, Left Choice Trials');
%         ylabel('Reward MI, Right Choice Trials');
%     end
    
    subplot(1,3,3);
    boots1 = bootstrp(5000,@mean,reward_MI_avg{1});
    boots2 = bootstrp(5000,@mean,reward_MI_avg{2});    
    err=[abs(prctile(boots1,[2.5;97.5])-mean(reward_MI_avg{1})),abs(prctile(boots2,[2.5;97.5])-mean(reward_MI_avg{2}))];
    errorbar([1 2],[mean(reward_MI_avg{1}) mean(reward_MI_avg{2})],err(1,:),err(2,:),'LineStyle','none','CapSize',0.00001,'Marker','o');
    set(gca,P.axes_properties{:},'xtick',[ 1 2],'xlim',[0.5 2.5],'xticklabel',{'untagged' 'tagged'},'ylim',[-0.1 0.1],'ytick',[-0.1:0.05:0.1]);
    ylabel('Reward MI');
    line(xlim,[0 0],'color','k','LineStyle',':');
    
    
    
    
end