function dspec = build_dspec_for_pbups(dspec,covariates,cellno)
    % Build 'designSpec' which specifies how to generate the design matrix
    % Each covariate to include in the model and analysis is specified.
    
    %% make condition lambda functions
    left_cond = @(trial) (~trial.pokedR);
    right_cond = @(trial) (trial.pokedR);
    hit_cond = @(trial) (trial.ishit);
    error_cond = @(trial) (~trial.ishit);
    left_hit_cond = @(trial) (trial.ishit & ~trial.pokedR);
    right_hit_cond = @(trial) (trial.ishit & trial.pokedR);
    left_error_cond = @(trial) (~trial.ishit & ~trial.pokedR);
    right_error_cond = @(trial) (~trial.ishit & trial.pokedR);
    
    events = fieldnames(dspec.expt.type);
    for i=1:length(events)
        if ~strcmp(dspec.expt.type.(events{i}),'timing')
            is_timing(i)=false;
        else
            is_timing(i)=true;
        end
    end
    timings = events(is_timing);
    if ischar(covariates)
        covariates={covariates};
    end
    covariates=unique(regexprep(covariates,'(.*[a-z])[0-9].*','$1'));     % remove trailing numbers because i search through all matching timing events anyway
    %% define click basis
    % good compromise that can explain posterior (fast) and anterior striatum (slow) click responses. also now with the offset it starts at zero.
    click_basis=basisFactory.makeNonlinearRaisedCos(6,1,[20 399],10);   
    
    %% loop over covariates
    for i=1:length(covariates)
        
        switch covariates{i} 
            
            case 'spike_history' 
                dspec = buildGLM.addCovariateSpiketrain(dspec, 'spike_history', ['sptrain',num2str(cellno)], 'History filter');
                
            case 'cpoke_in'
                bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 2500, 8, dspec.expt.binfun);
                %bs=basisFactory.makeNonlinearRaisedCos(7,1,[-2000 1800],500);    
                dspec = buildGLM.addCovariateTiming(dspec, 'cpoke_in','cpoke_in', 'Acausal filter aligned to cpoke_in', bs,-1250);  
                
            case 'cpoke_in_left'
                if ismember('cpoke_in',covariates)
                    warning('Cpoke_in and cpoke_in_left are redundant covariates.');
                end
                bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 2500, 8, dspec.expt.binfun);    
                %bs=basisFactory.makeNonlinearRaisedCos(7,1,[-2000 1800],350);    
                dspec = buildGLM.addCovariateTiming(dspec, 'cpoke_in_left','cpoke_in', 'Acausal filter aligned to cpoke_in on left choice trials', bs,-1000,left_cond);    

            case 'cpoke_in_right'
                if ismember('cpoke_in',covariates)
                    warning('Cpoke_in and cpoke_in_right are redundant covariates.');
                end                
                bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 2500, 8, dspec.expt.binfun);    
                %bs=basisFactory.makeNonlinearRaisedCos(7,1,[-2000 1800],350);    
                dspec = buildGLM.addCovariateTiming(dspec, 'cpoke_in_right','cpoke_in', 'Acausal filter aligned to cpoke_in on right choice trials', bs,-1000,right_cond);    

            case 'cpoke_out' 
                bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1500, 5, dspec.expt.binfun);        
                %bs=basisFactory.makeNonlinearRaisedCos(5,1,[-1000 1000],200);    
                dspec = buildGLM.addCovariateTiming(dspec, 'cpoke_out','cpoke_out', 'Acausal filter aligned to cpoke_out', bs,-750);

            case 'cpoke_out_left'
                if ismember('cpoke_out',covariates)
                    warning('Cpoke_out and cpoke_out_left are redundant covariates.');
                end                
                bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1500, 5, dspec.expt.binfun);        
                %bs=basisFactory.makeNonlinearRaisedCos(5,1,[-1000 1000],200);    
                dspec = buildGLM.addCovariateTiming(dspec, 'cpoke_out_left','cpoke_out', 'Acausal filter aligned to cpoke_out on left choice trials', bs,-750,left_cond);

            case 'cpoke_out_right'  
                if ismember('cpoke_out',covariates)
                    warning('Cpoke_out and cpoke_out_right are redundant covariates.');
                end                   
                bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1500, 5, dspec.expt.binfun);            
                %bs=basisFactory.makeNonlinearRaisedCos(5,1,[-1500 1500],200);    
                dspec = buildGLM.addCovariateTiming(dspec, 'cpoke_out_right','cpoke_out', 'Acausal filter aligned to cpoke_out on right choice trials', bs,-750,right_cond);    

            case 'stereo_click'        
                dspec = buildGLM.addCovariateTiming(dspec, 'stereo_click','stereo_click', 'Causal filter aligned to clicks_on on left choice trials', click_basis,0);

            case 'stereo_click_left'   
                if ismember('stereo_click',covariates)
                    warning('stereo_click and stereo_click_left are redundant covariates.');
                end                   
                dspec = buildGLM.addCovariateTiming(dspec, 'stereo_click_left','stereo_click', 'Causal filter aligned to clicks_on on left choice trials', click_basis,0,left_cond);

            case 'stereo_click_right'      
                if ismember('stereo_click',covariates)
                    warning('stereo_click and stereo_click_right are redundant covariates.');
                end                    
                dspec = buildGLM.addCovariateTiming(dspec, 'stereo_click_right','stereo_click', 'Causal filter aligned to clicks_on on right choice trials', click_basis,right_cond);

            case 'left_clicks'
                left_click_timings = timings(contains(timings,'left_clicks'));
                for k=1:length(left_click_timings)
                    dspec = buildGLM.addCovariateTiming(dspec, left_click_timings{k},[], left_click_timings{k}, click_basis,0);
                end

            case 'right_clicks'
                right_click_timings = timings(contains(timings,'right_clicks'));
                for k=1:length(right_click_timings)
                    dspec = buildGLM.addCovariateTiming(dspec, right_click_timings{k},[], right_click_timings{k}, click_basis,0);
                end

            case 'all_clicks'
                all_click_timings = timings(contains(timings,'all_clicks'));                
                for k=1:length(all_click_timings)
                    dspec = buildGLM.addCovariateTiming(dspec, all_click_timings{k},[], all_click_timings{k}, click_basis,0);
                end
                
            case 'spoke_left'
                %bs=basisFactory.makeNonlinearRaisedCos(5,1,[-1000 1000],200); 
                bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1000, 5, dspec.expt.binfun);                
                dspec = buildGLM.addCovariateTiming(dspec, 'spoke_left','spoke', 'Acausal filter aligned to spoke on left choice hit trials', bs,-500,left_cond);
    
            case 'spoke_right'
                %bs=basisFactory.makeNonlinearRaisedCos(5,1,[-1000 1000],200); 
                bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1000, 5, dspec.expt.binfun);               
                dspec = buildGLM.addCovariateTiming(dspec, 'spoke_right','spoke', 'Acausal filter aligned to spoke on right choice hit trials', bs,-500,right_cond);                 
                
            case 'spoke_left_hit'
                if ismember('spoke_left',covariates)
                    warning('spoke_left and spoke_left_hit are redundant covariates.');
                end                       
                %bs=basisFactory.makeNonlinearRaisedCos(5,1,[-1000 1000],200); 
                bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1000, 5, dspec.expt.binfun);                
                dspec = buildGLM.addCovariateTiming(dspec, 'spoke_left_hit','spoke', 'Acausal filter aligned to spoke on left choice hit trials', bs,-500,left_hit_cond);
    
            case 'spoke_right_hit'
                if ismember('spoke_right',covariates)
                    warning('spoke_right and spoke_right_hit are redundant covariates.');
                end                    
                %bs=basisFactory.makeNonlinearRaisedCos(5,1,[-1000 1000],200); 
                bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1000, 5, dspec.expt.binfun);               
                dspec = buildGLM.addCovariateTiming(dspec, 'spoke_right_hit','spoke', 'Acausal filter aligned to spoke on right choice hit trials', bs,-500,right_hit_cond); 
    
           case 'spoke_left_miss'
                if ismember('spoke_left',covariates)
                    warning('spoke_left and spoke_left_miss are redundant covariates.');
                end                         
                %bs=basisFactory.makeNonlinearRaisedCos(5,1,[-1000 1000],200);  
                bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1000, 5, dspec.expt.binfun);                
                dspec = buildGLM.addCovariateTiming(dspec, 'spoke_left_miss','spoke', 'Acausal filter aligned to spoke on left choice miss trials', bs,-500,left_error_cond);
    
            case 'spoke_right_miss'
                if ismember('spoke_right',covariates)
                    warning('spoke_right and spoke_right_miss are redundant covariates.');
                end                       
                %bs=basisFactory.makeNonlinearRaisedCos(5,1,[-1000 1000],200);   
                bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 1000, 5, dspec.expt.binfun);                
                dspec = buildGLM.addCovariateTiming(dspec, 'spoke_right_miss','spoke', 'Acausal filter aligned to spoke on right choice miss trials', bs,-500,right_error_cond);   
                
            otherwise
                
                error('%s not a recognized covariate.',covariates{i});
        end
    end
end