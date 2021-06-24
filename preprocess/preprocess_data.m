function preprocess_data(subj, ROOTDIR, type)

%   This function does all necessary preprocessing steps in order to get
%   all further analyses ready

%   Copyright (C) December 2017
%   D. Pedrosa, University Hospital of Gieﬂen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

%% General settings
cd(ROOTDIR);
loaddir     = [ROOTDIR '\data\'];
load([loaddir '\patdat.mat']);                                   %#ok<LOAD> % this file loads the meta data
if strcmp(type, 'p'); tolom = subj{2}; else tolom = subj{1}; end            % selects whether (p) pateints or controls (c) are analysed (see (type))

steps2apply = 3;                                                            % three steps available: (1):
hpf         = [.1 3];                                                       % high-pass filter frequency
lpf         = 30;
frsp        = 5000/200;                                                     % factor at which data was resampled
type_calc   = 'erp';
outdir = fullfile(loaddir, 'data_merged');                                  % directory at which data will be saved
if ~exist(outdir, 'dir'); mkdir(outdir); end

for np = tolom % loop through all subjects of one group
    if strcmp(type, 'p')
        disp(['the patient actually being computed is ' patient(np).name]); % this line provides the patient currently being analyzed in the command window
        code_participant = upper(patient(np).code);                         % relevant information for later in the next few lines
        bad_trials = patient(np).badtrials;
    else
        disp(['the subject actually being computed is ' control(np).name]); % this line provides the patient currently being analyzed in the command window
        code_participant = upper(control(np).code);                         % relevant information for later in the next few lines
        bad_trials = control(np).badtrials;
    end
    
    %% Starts with the different steps available (see general settings)
    for steps = steps2apply % loop through all steps to to (see general settings)
        switch (steps)
            case 1 % start filtering data
                filename_clean = {strcat('datclean_', code_participant, '_WO.mat'), ...
                    strcat('datclean_', code_participant, '_ALC.mat')};
                filename_preproc = {strcat('datpreproc_', code_participant, '_WO.mat'), ...
                    strcat('datpreproc_', code_participant, '_ALC.mat')};
                
                for c = 1:2 % loop through both conditions
                    if ~exist(strcat(outdir, filename_preproc{c}), 'file')     % the next few lines check if data is already present and skips further processing if so to avoid redundancy
                        fprintf('\npre-preprocessing for %s already done, continuing with next step ...\n', code_participant )
                        continue
                    else
                        load(strcat(outdir, filename_clean{c})); %#ok<LOAD> % this line loads the cleaned data into workspace (for details see clean_data.m)
                        data_clean.elec = ...                               % this block assigns the standard electrode positions to the data
                            ft_read_sens([ROOTDIR 'wcst\brainamp.sfp']);
                        
                        [data_clean.label,I] = sort(data_clean.label);
                        data_clean.trial{1,1} = data_clean.trial{1,1}(I,:);
                        
                        % Do the pre-processing in two steps: (1)
                        % preprocess for the ERP estimation, (2) preprocess
                        % for the TRF computation
                        
                        cfg = [];
                        cfg.reref       = 'yes';                            % the next few lines are intended to create average-referenced data
                        cfg.refchannel  = 'EEG';                            % to average reference data, the refchannel is set to 'all'
                        cfg.detrend     = 'yes';
                        cfg.bpfilter    = 'yes';                            % band-pass filter between .1 and 30Hz
                        cfg.bpfilttype  = 'firws';                          % Filter-type: FIR (window-scaled), for visualisation set cfg.plotfiltresp to 'yes'
                        cfg.bpfreq      = [hpf(1) lpf];                     % filter frequencies
                        cfg.plotfiltresp= 'no';
                        data_preproc{1} = ft_preprocessing(cfg, data_clean);
                        
                        cfg = [];
                        cfg.reref       = 'yes';                            % the next few lines are intended to create average-referenced data
                        cfg.refchannel  = 'EEG';                            % to average reference data, the refchannel is set to 'all'
                        cfg.detrend     = 'yes';
                        cfg.hpfilter    = 'yes';                            % in a second step, data is filtered in order to make TFR estimations possible, therefore, only hpf is used (lpf may be applied later)
                        cfg.hpfilttype  = 'firws';                          % Filter-type: FIR (window-scaled), for visualisation set cfg.plotfiltresp to 'yes'
                        cfg.hpfreq      = hpf(2);                           % cutoff frequency
                        cfg.plotfiltresp= 'no';
                        data_preproc{2} = ft_preprocessing(cfg, data_clean);
                        
                        save(strcat(outdir, filename_preproc{c}), ...
                            'data_preproc', '-v7.3');                           % saves the preprocessed EEG to a file
                    end
                end
                
            case 2 % this step rereferences the data to an average electrode
                % Defines the filename from where data is loaded (filename_clean) and
                % where/what filename to save
                filename_clean = {strcat('datpreproc_', code_participant, '_WO.mat'), ...
                    strcat('datpreproc_', code_participant, '_ALC.mat')};
                filename_trialdef = {strcat('trialdef_', code_participant, '_WO.mat'), ...
                    strcat('trialdef_', code_participant, '_ALC.mat')};
                
                switch type_calc
                    case 'erp'
                        filename_epoched = strcat('data_merged_', code_participant, '.mat');
                        idx_file = 1;
                    case 'tfr'
                        filename_epoched = strcat('data_merged_TFR_', code_participant, '.mat');
                        idx_file = 2;
                end
                
                for c = 1:2 % loop through both conditions
                    if ~exist(strcat(outdir, filename_epoched), 'file')     % the next few lines check if data is already present and skips further processing if so to avoid redundancy
                        fprintf('\nepoching for %s already done, continuing with next step ...\n', code_participant )
                        continue
                    else
                        load(strcat(outdir, filename_clean{c}));            % this line loads the cleaned data into workspace (for details see clean_data.m)
                        load(strcat(outdir, filename_trialdef{c}));         % this line loads the trial definitions so that data may be cut into chunks to be processed
                        
                        data_clean = data_preproc{idx_file};
                        
                        cfg             = [];
                        cfg.trl         = trialdef.trl;
                        cfg.trl(:,1:2)  = round(cfg.trl(:,1:2)./frsp);      % the next two lines intend to "resample" the trialdefinition in a
                        cfg.trl(:,3)    = round(cfg.trl(:,3)./frsp);        % admittedly not very elegant way, but it works!
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
                        
                        cfg.minlength   = 'maxperlen';                          % this ensures all resulting trials are equal length
                        data_epoched{c} = ft_redefinetrial(cfg, data_clean);
                        
                        if c == 2 % in the second condition, data from the ALC condition is changed to make it identifiable and data is merged to one file
                            data_epoched{c}.trialinfo = ...
                                data_epoched{c}.trialinfo+100;              % adds a n arbitrary value to the trialinfo data, in order to make it unambiguos later
                            
                            cfg = [];                                       % merges both conditions to a single file
                            data_merged = ...
                                ft_appenddata(cfg, data_epoched{1}, ...
                                data_epoched{2}); clear data_epoched
                            save(strcat(outdir, filename_epoched), ...
                                'data_merged', '-v7.3');                    % saves the merged EEG and acc data -to one file
                        end
                    end
                end
                
            case 3 % this step extracts the badtrials and removes them from data; in case not yet identified, data is plotted
                
                % Defines the filename from where data is loaded (filename_epoched) and
                % where/what filename to save

                switch type_calc
                    case 'erp'
                        filename_final = strcat('data_final_', code_participant, '.mat');
                        filename_epoched = strcat('data_merged_', code_participant, '.mat');
                        idx_file = 1;
                    case 'tfr'
                        filename_final = strcat('data_final_TFR_', code_participant, '.mat');
                        filename_epoched = strcat('data_merged_TFR_', code_participant, '.mat');
                        idx_file = 2;
                end

                if exist(strcat(outdir, filename_final), 'file')     % the next few lines check if data is already present and skips further processing if so to avoid redundancy
                    fprintf('\nremoving badtrials for %s already done, continuing with next step ...\n', code_participant )
                    continue
                else
                    load(fullfile(outdir, filename_epoched));            % this line loads the cleaned data into workspace (for details see clean_data.m)
                    
                    % The next few lines are intended to remove
                    % anticipations or "bad blocks"
                    idx_fst = [1, 101];                                     % trialinfo that marks the beginning of a new block
                    idx_fstrl = []; clear tmp_trials                        % initialise variable to fill with content during loop
                    for c = 1:numel(idx_fst) % loop through both "conditions"
                        idx_fstrl = [idx_fstrl; ...                         % creates a vector of "first trials"
                            find(data_merged.trialinfo == idx_fst(c))]; %#ok<AGROW>
                    end
                    idx_bad = remove_badblocks(data_merged, idx_fstrl, ...  % this sctipt reads out the anticipation errors and discards the entire block of these
                        code_participant);                                  % furthermore trial errors are read out and indexed (in total a variable (idx_bad) is created
                    
                    idx9131 = [];                                           % every experiment started with the combination of code 91 and 31, so this is identified and later added to the idx_bad
                    idx9131 = [idx9131; strfind(data_merged.trialinfo.', [91 31]); ...
                        [strfind(data_merged.trialinfo.', [91 31])+1]]; %#ok<AGROW,NBRAK>
                    idx9131 = [idx9131; strfind(data_merged.trialinfo.', [191 131]); ...
                        [strfind(data_merged.trialinfo.', [191 131])+1]]; %#ok<AGROW,NBRAK>
                    
                    % Remove blocks with anticipation or errors from merged
                    % file
                    cfg = [];
                    cfg.trials = setdiff(1:length(data_merged.trialinfo), ...
                        [idx9131; idx_bad]);
                    data_merged = ft_selectdata(cfg, data_merged);
                    
                    if isempty(bad_trials) % if the bad trials are not yet identified, this part of the code plots all trials in order to search for bad ones
                        
                        idx = [{1:7}, {101:107}];                               % code at which a new trial was initiaded (from which ERPs will be estimated later)
                        idx_trl = []; clear tmp_trials                          % initialise variable to fill with content during loop
                        for c = 1:numel(idx)
                            tmp_trials = ismember(data_merged.trialinfo, idx{c});
                            idx_trl = [idx_trl; find(tmp_trials)];
                            clear tmp_trials
                        end
                        
                        %% Start plotting
                        cfg = [];
                        cfg.trials = idx_trl;                               % select only the trials that are useful for detecting artifacts
                        data_plot = ft_selectdata(cfg,data_merged);
                        
                        cfg = [];                           % and after interpolation of artefact channels
                        cfg.channel                 = 'EEG';
                        cfg.viewmode                = 'vertical';
                        cfg.preproc.bpfilter        = 'yes';                        % filtering and pre-processing of the data
                        cfg.preproc.bpfreq          = [4 70];
                        cfg.preproc.bsfilter        = 'yes';                        % filtering and pre-processing of the data
                        cfg.bsfreq                  = [48 52];
                        cfg.preproc.demean          = 'yes';
                        cfg.preproc.detrend         = 'yes';
                        cfg.layout                  = 'eeg1005.lay';                % specifies the layout file that should be used for plotting
                        ft_databrowser(cfg, data_plot)
                        keyboard
                        % at this point, the bad trials should be entered and
                        % aves to the patdat.m file
                        
                        clear data_plot
                    end
                    
                    trial_count = 1:3:numel(data_merged.trialinfo);         % because the badtrials are checked on the reduced list, this vector is necessary to remove the "right" trials in a later step
                    
                    % Remove blocks with anticipation or errors from merged file
                    cfg = [];
                    cfg.trials = setdiff(1:length(data_merged.trialinfo), ...
                        trial_count(bad_trials));
                    data_final = ft_selectdata(cfg, data_merged); %#ok<NASGU>
                    
                    save(strcat(outdir, filename_final), ...
                        'data_final', '-v7.3');                    % saves the merged EEG and acc data -to one file
                    
                end
        end
    end
end