function compare_choice_MI(stats1,stats2,params)
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
    inputnames = {inputname(1) inputname(2)};
    
    %% compare choice selectivity, pre move
    choice_MI{1} = cat(1,stats1.choice_MI_cpoke_out_prepost);
    choice_MI{2} = cat(1,stats2.choice_MI_cpoke_out_prepost);    
    choice_MI_pval{1} = cat(1,stats1.choice_MI_pval);
    choice_MI_pval{2} = cat(1,stats2.choice_MI_pval);       
    subplot(3,2,1);    
    h(1)=significance_histogram(choice_MI{1}(:,1),choice_MI_pval{1}(:,1),'bin_limits',[-1 1],'nbins',params.nbins,'FaceColor',params.color1,'Normalization','probability');
    subplot(3,2,2);    
    h(2)=significance_histogram(choice_MI{2}(:,1),choice_MI_pval{2}(:,1),'bin_limits',[-1 1],'nbins',params.nbins,'FaceColor',params.color2,'Normalization','probability');    
    xlabel('Choice MI, pre-MOVE');
    ylabel('Frac. Cells');

    %% compare choice selectivity, post move 
    subplot(3,2,3);    
    h(1)=significance_histogram(choice_MI{1}(:,2),choice_MI_pval{1}(:,2),'bin_limits',[-1 1],'nbins',params.nbins,'FaceColor',params.color1,'Normalization','probability');
    subplot(3,2,4);    
    h(2)=significance_histogram(choice_MI{2}(:,2),choice_MI_pval{2}(:,2),'bin_limits',[-1 1],'nbins',params.nbins,'FaceColor',params.color2,'Normalization','probability');    
    xlabel('Choice MI, post-MOVE');
    ylabel('Frac. Cells');
    
    for i=1:2
        subplot(3,2,4+i);
        myscatter(choice_MI{i}(:,1),choice_MI{i}(:,2));hold on;
        scatter(mean(choice_MI{i}(:,1)),mean(choice_MI{i}(:,2)),'or','filled','SizeData',50);
        set(gca,P.axes_properties{:},'xlim',[-1 1],'ylim',[-1 1]);
        xlabel('Choice MI, pre-MOVE');
        ylabel('Choice MI, post-MOVE');     
        grid on;
        line([ 0 0],ylim,'Color','k');
        line(xlim,[ 0 0],'Color','k'); 
        title('');    
        text(0.6,0.8,{'CONTRA PRE-MOVE','CONTRA POST-MOVE',sprintf('%g%%',round(100*mean(choice_MI{i}(:,1)>0 & choice_MI{i}(:,2)>0)))},'FontSize',8,'HorizontalAlignment','center','Color','r','FontWeight','bold')
        text(0.6,-0.8,{'CONTRA PRE-MOVE','IPSI POST-MOVE',sprintf('%g%%',round(100*mean(choice_MI{i}(:,1)>0 & choice_MI{i}(:,2)<0)))},'FontSize',8,'HorizontalAlignment','center','Color','r','FontWeight','bold')
        text(-0.6,0.8,{'IPSI PRE-MOVE','CONTRA POST-MOVE',sprintf('%g%%',round(100*mean(choice_MI{i}(:,1)<0 & choice_MI{i}(:,2)>0)))},'FontSize',8,'HorizontalAlignment','center','Color','r','FontWeight','bold')
        text(-0.6,-0.8,{'IPSI PRE-MOVE','IPSI POST-MOVE',sprintf('%g%%',round(100*mean(choice_MI{i}(:,1)<0 & choice_MI{i}(:,2)<0)))},'FontSize',8,'HorizontalAlignment','center','Color','r','FontWeight','bold')
    end
    
end