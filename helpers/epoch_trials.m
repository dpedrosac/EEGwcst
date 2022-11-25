function epoch_trials(filename_preproc, filename_save, processed_subj, ...
    wdir, ROOTDIR, metadata, type_calc)

%   This function epoches data and marks the "incomplete runs" that is
%   where errors ocurred or anticipations were traceable

%   Copyright (C) June 2021, modified July 2021 and October 2022
%   D. Pedrosa, University Hospital of Gießen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

conds           = {'WO', 'ALC'};
data_epoched    = cell(1,numel(conds));                                     % pre-allocate space
frsp            = 5000/200;                                                 % sampling rate (cf. ./read_data.m)
idx_type        = 1;                                                        % default = 1 is 'erp'
if strcmp(type_calc, 'tfr'); idx_type=2; end

for c = 1:2 % loop through both conditions, WO and ALC
    ev_all = load(fullfile(ROOTDIR, 'data', 'header_and_events', ...
        sprintf('events_%s_%s.mat', conds{c}, processed_subj)));%#ok<*LOAD> % loads events with marked "wrong" runs
    ev_all = ev_all.events; %clear events;
    load(fullfile(wdir, 'data_preprocessed', filename_preproc{c}));         % this line loads the cleaned data into workspace (for details see clean_data.m)
    load(fullfile(wdir, 'trialdef', sprintf('trialdef_%s_%s.mat', ...
        processed_subj, conds{c})));                                        % this line loads the trial definitions so that data may be cut into chunks to be processed
    clear data_clean
    data_clean = data_preproc{idx_type};                        %#ok<USENS> % selects whether erp (1) or tfr (2) data is used

    %if ~isfield(ev_all, 'incomplete')
        [tmp, ~] = select_runs(processed_subj, ROOTDIR);                   % writes the erroneous runs into the events-file
    %end
    
    %% Start selecting trials and merge data
    cfg = [];
    cfg.trl         = trialdef.trl;
    cfg.trl(:,1:2)  = round(cfg.trl(:,1:2)./frsp);                          % the next two lines intend to "resample" the trialdefinition in a
    cfg.trl(:,3)    = round(cfg.trl(:,3)./frsp);                            % admittedly not very elegant way, but it works!
    
    idx_low         = find(cfg.trl(:,1) - cfg.trl(:,2) == -499);
    idx_high        = find(cfg.trl(:,1) - cfg.trl(:,2) == -501);
    if ~isempty(idx_low); cfg.trl(idx_low,2) = cfg.trl(idx_low,2)+1; end    % the next lines prevent
    if ~isempty(idx_high); cfg.trl(idx_high,2) = cfg.trl(idx_high,2)+1; end % problems at the "edges"
    
    cfg.minlength   = 'maxperlen';                                          % this ensures all resulting trials are equal length
    data_epoched{c} = ft_redefinetrial(cfg, data_clean);
    
    % Change trial number if run was marked as 'wrong'/'incomplete'
    incomplete_events = [ev_all.incomplete].';
    incomplete_events = incomplete_events(3:end);
    if length(incomplete_events) ~= length(data_epoched{c}.trialinfo)
        incomplete_events = incomplete_events(1:end-1);
    end
    
    idx_incomplete = find(incomplete_events~=0);
    idx_incomplete = idx_incomplete(idx_incomplete<=...                     % this avoids that trials at the end of the run crash the script
        length(data_epoched{c}.trialinfo));

    % data_epoched{c}.trialinfo(idx_incomplete) = ...
    %    data_epoched{c}.trialinfo(idx_incomplete) + 1000;
    
    % Here, *badtrials* saved in the metadata file are removed""
    
    cfg = [];
    cfg.trials = setdiff(1:length(data_epoched{c}.trialinfo), ...
        metadata.badtrials{c});
    data_epoched{c} = ft_selectdata(cfg, data_epoched{c});
    
    if c == 2 % in the second condition, data from the ALC condition is changed to make it identifiable and data is merged to one file
        data_epoched{c}.trialinfo = data_epoched{c}.trialinfo+100;          % adds arbitrary value to trialinfo data (ALC), to make it unambiguos
        
        cfg = []; % merges both conditions to a single file
        cfg.keepsampleinfo  = 'no';
        data_merged         = ft_appenddata(cfg, data_epoched{1}, ...
            data_epoched{2}); clear data_epoched
        %data_merged = ft_struct2single(data_merged);
        save(filename_save, 'data_merged', '-v7.3');        % saves merged EEG and acc data -to one file
    end
end