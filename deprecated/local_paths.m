function out_wdir = local_paths
% sets and returns local paths, also for fieldtrip

if strcmp(getenv('username'), 'dpedr')
    out_wdir     = 'D:\skripte\lambda\';   
    out_data_dir = 'D:\skripte\lambda\data';   
    out_ftdir    = 'D:\skripte\fieldtrip';    
end
if strcmp(getenv('USER'), 'urs')
    rootdir = '/home/urs/sync/projects/wcst_eeg';
    out_ftdir     = '/opt/fieldtrip/fieldtrip-20210507';
    out_wdir      = strcat(rootdir,'/analysis');    
    out_data_dir  = strcat(rootdir,'/data');    
    addpath(strcat(rootdir,'/analysis/localfunctions'));
end

% add project paths
addpath(genpath(out_wdir));

% set fieldtrip paths using the recommended way
addpath(out_ftdir);
ft_defaults;
