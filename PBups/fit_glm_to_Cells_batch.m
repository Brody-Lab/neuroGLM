function res = fit_glm_to_Cells_batch(rat,varargin)
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('remake',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('sort_fldr','X:\RATTER\PhysData\Mountainsort\Adrian',@(x)isdir(x));
    p.parse(varargin{:});
    params=p.Results;
    if ~isdir([params.sort_fldr,filesep,rat])
        error('\nNo data in %s for rat %s.',params.sort_fldr,rat);
    end
    possible_sessions = dir([params.sort_fldr,filesep,rat]);
    possible_sessions = {possible_sessions.name};
    possible_sessions = possible_sessions(strncmp(possible_sessions,'20',2));
    for i=9:9
        fprintf('\nWorking on rat %s, folder %s.\n',rat,possible_sessions{i});
        msorted_filename = fixfilesep([params.sort_fldr,filesep,rat,filesep,possible_sessions{i},filesep,'Msorted_',rat,'_',possible_sessions{i},'.mat']); 
        if ~exist(msorted_filename,'file')
            res{i} = 'no msorted';
            fprintf([res{i},'\n']);
            continue
        end
        fprintf('Loading %s.\n',msorted_filename);
        file_contents = load(msorted_filename);
        if ~any(file_contents.Msorted.n_cluster) 
            res{i} = 'no cells';
            fprintf([res{i},'\n']);
            continue
        end
        try
            try
                file_contents.Msorted.Trials = add_shuffled_click_times(file_contents.Msorted.Trials);
                file_contents.Cells = PB_make_PETH(file_contents.Msorted);
            catch
                file_contents.Cells.Trials = add_shuffled_click_times(file_contents.Cells.Trials);
                file_contents.Cells = PB_make_PETH(file_contents.Cells);
            end
            [a,b,c] = fileparts(msorted_filename);                             
            fit_glm_to_Cells(file_contents.Cells,'use_parallel',false,'save',true,varargin{:},'save_path',[a,filesep,b,'_glmfits_3clicks.mat']);            
            res{i} = 'success';                
        catch ME
            res{i} = ME;
            fprintf('Encountered error in glm fitting on rat %s, folder %s.\n',rat,possible_sessions{i});
        end 
    end
end