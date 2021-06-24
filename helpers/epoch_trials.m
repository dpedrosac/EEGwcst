<<<<<<< HEAD
function epoch_trials(filename_clean, processed_subj, datadir, idx_file)
=======
function epoch_trials(filename_preproc, filename_save, processed_subj, datadir, idx_file)
>>>>>>> beab179d0fc15269759e57e1d637962ca00d14c4

%   This function epoches data and marks the "incomplete runs" that is
%   where errors ocurred or anticipations were traceable

%   Copyright (C) June 2021
%   D. Pedrosa, University Hospital of Gießen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

conds = {'WO', 'ALC'};
data_dir = fullfile(datadir, 'data_clean');
<<<<<<< HEAD
data_epoched = cell(1,2);
for c = 1:2 % loop through both conditions, WO and ALC
    load(fullfile(datadir, 'header_and_events', ...
        sprintf('events_%s_%s.mat', conds{c}, processed_subj)));            % loads events with marked "wrong" runs
    load(fullfile(data_dir, filename_clean{c}));                            % this line loads the cleaned data into workspace (for details see clean_data.m)
=======
data_epoched = cell(1,numel(conds));
frsp = 5000/200;

for c = 1:2 % loop through both conditions, WO and ALC
    load(fullfile(datadir, 'header_and_events', ...
        sprintf('events_%s_%s.mat', conds{c}, processed_subj)));            % loads events with marked "wrong" runs
    load(fullfile(data_dir, filename_preproc{c}));                            % this line loads the cleaned data into workspace (for details see clean_data.m)
>>>>>>> beab179d0fc15269759e57e1d637962ca00d14c4
    load(fullfile(data_dir, sprintf('trialdef_%s_%s.mat', ...
        processed_subj, conds{c})));                                        % this line loads the trial definitions so that data may be cut into chunks to be processed
    data_clean = data_preproc{idx_file};                        %#ok<USENS> % selects whether erp (1) or tfr (2) data is used
    %% ==> here a mechanism selecting bad runs and adding specific codes is required (subfunction!!)

    %% Start selecting trials and merge data
    cfg = [];
    cfg.trl         = trialdef.trl;
<<<<<<< HEAD
    cfg.trl(:,1:2)  = round(cfg.trl(:,1:2)./frsp);      % the next two lines intend to "resample" the trialdefinition in a
    cfg.trl(:,3)    = round(cfg.trl(:,3)./frsp);        % admittedly not very elegant way, but it works!
=======
    cfg.trl(:,1:2)  = round(cfg.trl(:,1:2)./frsp);                          % the next two lines intend to "resample" the trialdefinition in a
    cfg.trl(:,3)    = round(cfg.trl(:,3)./frsp);                            % admittedly not very elegant way, but it works!
>>>>>>> beab179d0fc15269759e57e1d637962ca00d14c4
    idx_low         = find(cfg.trl(:,1) - ...
        cfg.trl(:,2) == -499);
    idx_high         = find(cfg.trl(:,1) - ...
        cfg.trl(:,2) == -501);
    if ~isempty(idx_low)
        cfg.trl(idx_low,2) = cfg.trl(idx_low,2)+1;
    end
    
    if ~isempty(idx_high)
        cfg.trl(idx_high,2) = cfg.trl(idx_high,2)+1;
    end
    
<<<<<<< HEAD
    cfg.minlength   = 'maxperlen';                          % this ensures all resulting trials are equal length
=======
    cfg.minlength   = 'maxperlen';                                          % this ensures all resulting trials are equal length
>>>>>>> beab179d0fc15269759e57e1d637962ca00d14c4
    data_epoched{c} = ft_redefinetrial(cfg, data_clean);
    
    if c == 2 % in the second condition, data from the ALC condition is changed to make it identifiable and data is merged to one file
        data_epoched{c}.trialinfo = data_epoched{c}.trialinfo+100;          % adds arbitrary value to trialinfo data (ALC), to make it unambiguos
        
        cfg = []; % merges both conditions to a single file
        data_merged = ft_appenddata(cfg, data_epoched{1}, ...
            data_epoched{2}); clear data_epoched
<<<<<<< HEAD
        save(strcat(outdir, filename_epoched), 'data_merged', '-v7.3');     % saves merged EEG and acc data -to one file
=======
        save(strcat(outdir, filename_save), 'data_merged', '-v7.3');        % saves merged EEG and acc data -to one file
>>>>>>> beab179d0fc15269759e57e1d637962ca00d14c4
    end
end