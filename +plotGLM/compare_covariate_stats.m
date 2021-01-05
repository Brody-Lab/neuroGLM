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
    inputnames = {inputname(1) inputname(2)};
    
    %% compare reward selectivity, pre move
    reward_MI_left{1} = cat(1,stats1.reward_MI_left);
    reward_MI_right{1} = cat(1,stats1.reward_MI_right);
    reward_MI_left{2} = cat(1,stats2.reward_MI_left);
    reward_MI_right{2} = cat(1,stats2.reward_MI_right);
    reward_MI_avg{1} = (reward_MI_left{1} + reward_MI_right{1})/2;
    reward_MI_avg{2} = (reward_MI_left{2} + reward_MI_right{2})/2;    
   

    %% compare reward selectivity for both sides
    for i=1:2
        subplot(1,3,i);
        myscatter(reward_MI_left{1},reward_MI_right{2});
        xlabel('Reward MI, Left Choice Trials');
        ylabel('Reward MI, Right Choice Trials');
    end
    
    
    
    
    
    
end