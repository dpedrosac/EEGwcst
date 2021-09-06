function estimate_erp_averages_generic(all_subj, ROOTDIR, wdir, toi, ending)

%   This function loads all data and estimates the metrics needed, saving
%   them as structures (avg*filename*). As a new feature, there is now a
%   conversion with ft_struct2single(x)

%   Copyright (C) July 2021
%   D. Pedrosa, University Hospital of Gießen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

% toi = {             % Different trial number available; for further details cf: ./paradigm/WCST-Tremorparadigma.sce
%     [21,22,23], ... % Early trials (wo/ alcohol), formerly [2,3]
%     [24,25], ...    % Late trials (wo/ alcohol), formerly [6,7]
%     [102:103], ...  % Early trials (w/ alcohol)
%     [106:107], ...  % Late trials (w/ alcohol)
%     [21:25], ...    % Right trials (wo/ alcohol)
%     [121:125], ...  % Right trials (w/ alcohol)
%     [10, 20], ...   % Shift trials(wo/ alcohol)
%     [110, 120]};    % Shift trials(w/ alcohol)                  %#ok<*NBRAK>

%% Start isolating trials of interest (toi) for all subjects

if nargin > 5
    fprinft('\neither trials of interest (toi) or ending or both missing, please double-check!');
    return
end

inputdir = fullfile(wdir, 'data_final');
outdir = fullfile(ROOTDIR, 'data');

tic;
avg = cell(1,numel(all_subj));
parfor proc = 1:numel(all_subj)
    fprintf('\nthe subject being processed is: %s \n', ...
        strcat('S', num2str(all_subj(proc))));
    filename_data = fullfile(inputdir, ...
        ['data_final_erp_S', num2str(all_subj(proc)),'.mat']);
    partemp = load(filename_data);
    
    % Extract trials of interest (toi)
    cfg = [];
    cfg.trials = find(ismember(partemp.data_final.trialinfo, toi));
    data_temp = ft_selectdata(cfg, partemp.data_final);
    if isempty(data_temp.trialinfo), continue; end
    
    cfg = [];
    cfg.bpfiler = 'yes';
    cfg.bpfreq = [.1 30];
    data_temp = ft_preprocessing(cfg, data_temp);
    
    cfg = [];
    cfg.keeptrials = 'yes';
    data_temp = ft_timelockanalysis(cfg, data_temp);
    data_temp.avg = nanmean(data_temp.trial, 1);
    data_temp = ft_struct2single(data_temp);                                % to save space, data is converted to "single"
    avg{proc} = data_temp;
end
toc

%% Save averages for all trials of interest (toi) defined above
save(fullfile(outdir, sprintf('avg_%s_erp.mat', ending)), 'avg', '-v7.3')

end