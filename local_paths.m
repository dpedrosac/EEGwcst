function out = local_paths
% sets and returns local paths, also for fieldtrip

if strcmp(getenv('username'), 'dpedrosa')
    out.wdir     = 'D:\skripte\lambda\';   
    out.data_dir = 'D:\skripte\lambda\data';   
    out.ftdir    = 'D:\skripte\fieldtrip';    
end
if strcmp(getenv('USER'), 'urs')
    rootdir = '/home/urs/sync/projects/wcst_eeg'
    out.ftdir     = '/opt/fieldtrip/fieldtrip-20210507';
    out.wdir      = strcat(rootdir,'/analysis');    
    out.data_dir  = strcat(rootdir,'/data');    
    addpath(strcat(rootdir,'/analysis/localfunctions'));
end


% add project paths
addpath(genpath(out.wdir));

% set fieldtrip paths using the recommended way
addpath(out.ftdir);
ft_defaults;
