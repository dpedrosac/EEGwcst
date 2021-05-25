function results_estimation(wdir, subj, metric)

%   This function estimates metrics for every subject, that may be defined
%   before

%   Copyright (C) March 2018
%   D. Pedrosa, University Hospital of Gieﬂen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

% general settings

% trls_oi = {2:3, 6:7, [10,20], 21:25, 102:103, 106:107, [110,120], 121:125}; % trials of interest
% names = {'earlyWO', 'lateWO', 'wrongWO', 'rightWO', ...
%     'earlyALC', 'lateALC', 'wrongALC', 'rightALC'};

trls_oi = {[10,20], 20, 21:25, [10,20:25], ...
    [110,120], 120, 121:125, [110, 120:125]}; % trials of interest
names = {'wrongWO1', 'wrongWO2', 'rightWO', 'allWO', ...
    'wrongALC1', 'wrongALC2', 'rightALC', 'allALC'};

%% Start esimation of /metric)
for k = 1:numel(names) % loop through all different events of interest
    for dId = 1:numel(subj) % loop through subjects available
        fprintf('\nthe subject being processed is: %s \n', strcat('S', num2str(subj(dId))));
        
        switch metric % switch which makes it possible to estrimate the ERP or the TFR
            case 'ERP'
                filename_save = fullfile(wdir, 'results', ...
                    strcat(names(k), '_ERP.mat'));
                %                 if exist(filename_save{1}, 'file'); continue; end
                clear data_final
                filename = strcat(wdir, 'S', num2str(subj(dId)), ...                % filename which is loaded in the neyt line
                    '\data_final_S', num2str(subj(dId)), '.mat'); load(filename);
                ERP{dId} = ftERP_estimate(data_final, trls_oi{k}, wdir, subj(dId), names{k});     % estimate the ERP from the data using only the trials defined in (trls_oi)
            case 'TFR'
                filename_save = fullfile(wdir, 'results', ...
                    strcat(names(k), '_TFR.mat'));
                %                 if exist(filename_save{1}, 'file'); continue; end
                clear data_final
                filename = strcat(wdir, 'S', num2str(subj(dId)), ...                % filename which is loaded in the neyt line
                    '\data_final_TFR_S', num2str(subj(dId)), '.mat'); load(filename);
                TFR{dId} = ftTFR_estimate(data_final, trls_oi{k}, wdir, subj(dId), names{k});     % estimate the TFR from the data using only the trials defined in (trls_oi)
        end
        
    end
end

end

function ERP_avg = ftERP_estimate(data_final, trlsnum, wdir, subj, name)
%   This function estimates ERP using the fieldtrip toolbox with the
%   settings defined below and saves the data in a file

%   Copyright (C) March 2018
%   D. Pedrosa, University Hospital of Gieﬂen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.


cfg = [];
cfg.trials      = find(ismember(data_final.trialinfo, trlsnum));            % use only the
cfg.feedback    = 'no';
data_select     = ft_selectdata(cfg,data_final);

cfg = [];
cfg.keeptrials          = 'yes';
cfg.covariance          = 'yes';
cfg.covariancewindow    = [-inf -.05];
cfg.feedback            = 'no';
ERP_avg                 = ft_timelockanalysis(cfg, data_select);

filename_ERP = fullfile(wdir, strcat(upper('s'), num2str(subj)), ...
    strcat('ERP_', name , strcat('_', upper('s'), num2str(subj)), '.mat'));
save(filename_ERP, 'ERP_avg', '-v7.3');
end

function TFR_avg = ftTFR_estimate(data_final, trlsnum, wdir, subj, name)
%   This function estimates TFR using the fieldtrip toolbox with the
%   settings defined below and saves the data in a file

%   Copyright (C) March 2018
%   D. Pedrosa, University Hospital of Gieﬂen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

neighbours = define_neighbours(data_final, wdir, 35, 0);                 % runs the function in which definition of neighbours is
ch_2plot = 'Pz';
idx_ch = neighbours(find(strcmp({neighbours(:).label}, ch_2plot))).neighblabel; %#ok<FNDSB>

cfg = [];
cfg.trials = find(ismember(data_final.trialinfo, trlsnum));
cfg.channel = idx_ch;
data_select = ft_selectdata(cfg, data_final);

cfg = [];
cfg.pad         = 'nextpow2';
cfg.output      = 'pow';
cfg.channel     = 'all';
cfg.keeptrials  = 'yes';
cfg.method      = 'mtmconvol';
cfg.foi         = 6:1:90;
cfg.t_ftimwin   = 7./cfg.foi;
cfg.tapsmofrq   = .3*cfg.foi;
cfg.toi         = 'all';%-.5:.05:1.5;
TFR_avg         = ft_freqanalysis(cfg, data_select);

filename_TFR = fullfile(wdir, strcat(upper('s'), num2str(subj)), ...
    strcat('TFR', name , strcat('_', upper('s'), num2str(subj)), '.mat'));
save(filename_TFR, 'TFR_avg', '-v7.3');
end