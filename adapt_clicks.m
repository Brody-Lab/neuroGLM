function  [left_adapted,right_adapted] = adapt_clicks(left_clicks,right_clicks,stereo_click,phi,tau_phi,within_stream)
    % stereo clicks don't need to be adapted because they are fit with their own
    % kernel anyway
    % assumes left_clicks and right_clicks are sorted
    % cross-stream adaptation by default
    % 1st three inputs must be row vectors
    
    left_adapted=left_clicks;
    right_adapted=right_clicks;
    if nargin<4
        return
    elseif nargin<5
        error('You provided phi but not tau_phi.');
    elseif nargin<6
        within_stream=false;
    end
    if abs(phi-1)<=eps || abs(tau_phi)<=eps
        return
    end
    if within_stream
        left_adapted = adapt_clicks_internal(left_clicks,phi,tau_phi);     
        right_adapted = adapt_clicks_internal(right_clicks,phi,tau_phi);       
    else
        stims=[left_clicks right_clicks stereo_click];
        idx=[ones(1,numel(left_clicks)) -ones(1,numel(right_clicks)) 0]; 
        [stims,sort_idx]=sort(stims);
        idx=idx(sort_idx);        
        adapted = adapt_clicks_internal(stims,phi,tau_phi);     
        left_adapted=adapted(idx==1);
        right_adapted=adapted(idx==-1);        
    end
end

function adapted = adapt_clicks_internal(stims,phi,tau_phi)
    % assumes stims is sorted
    adapted=ones(size(stims));     
    ici = diff(stims);
    expterm = exp(-ici/tau_phi);
    for i=2:length(stims)
        adapted(i) = 1 + expterm(i-1).*(adapted(i-1)*phi-1);     
        if ici(i-1)==0
            adapted(i-1) = 1 + expterm(i-1).*(adapted(i-1)*phi-1);     % if the ici is exactly zero (which happens rarely) don't let the arbitrarily first one win. adapt them both to each other.
        end
    end
end