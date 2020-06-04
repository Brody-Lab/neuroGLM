function [branch,commit] = return_git_status(repo_path)
    wd=cd;
    cd(repo_path);
    [~,branch] = system('git rev-parse --abbrev-ref HEAD');
    branch = deblank(branch);
    [~,commit] = system('git rev-parse HEAD');
    commit=deblank(commit);
    cd(wd); 
end