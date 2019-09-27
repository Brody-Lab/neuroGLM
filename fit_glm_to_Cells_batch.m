function res = fit_glm_to_Msorted_batch(rat,varargin)
    p=inputParser;
    p.KeepUnmatched=true;
    p.addParameter('remake',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    p.addParameter('mountainsort_fldr','X:\RATTER\PhysData\Mountainsort\Adrian',@(x)isdir(x));
    p.parse(varargin{:});
    params=p.Results;
    if ~isdir([params.mountainsort_fldr,filesep,rat])
        error('\nNo data in %s for rat %s.',params.mountainsort_fldr,rat);
    end
    possible_sessions = dir([params.mountainsort_fldr,filesep,rat]);
    possible_sessions = {possible_sessions.name};
    possible_sessions = possible_sessions(strncmp(possible_sessions,'20',2));
    for i=1:length(possible_sessions)
        fprintf('\nWorking on rat %s, folder %s.\n',rat,possible_sessions{i});
        msorted_filename = fixfilesep([params.mountainsort_fldr,filesep,rat,filesep,possible_sessions{i},filesep,'Msorted_',rat,'_',possible_sessions{i},'.mat']); 
        if ~exist(msorted_filename,'file')
            res{i} = 'no msorted';
            fprintf([res{i},'\n']);
            continue
        end
        fprintf('Loading %s.\n',msorted_filename);
        file_contents = load(msorted_filename);
        try
            fit_glm_to_Msorted(file_contents.Msorted,'use_parallel',true,'save',true,varargin{:});
            res{i} = 'success';
        catch ME
            res{i} = ME;
            fprintf('Encountered error in glm fitting on rat %s, folder %s.\n',rat,possible_sessions{i});
        end 
    end
end