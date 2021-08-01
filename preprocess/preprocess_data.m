function preprocess_data(subj, ROOTDIR, wdir, type)

%   This function does all necessary preprocessing steps in order to get
%   all further analyses ready, that is filtering, epoching and removing
%   bad trials

%   Copyright (C) December 2017, modified July 2021
%   D. Pedrosa, University Hospital of Gieﬂen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

%% General settings
cd(ROOTDIR);
loaddir     = [ROOTDIR '\data\'];
load([loaddir '\patdat.mat']);                                  %#ok<LOAD>  % this file loads the meta data
if strcmp(type, 'p'); tolom = subj{2}; else; tolom = subj{1}; end           % selects whether (p) pateints or controls (c) are analysed (see (type))
steps2apply = 3;                                                            % three steps available: (1): filter data (2): epoching and merging data (3): removing bad trials and 'runs'
type_calc   = 'erp';

for np = tolom % loop through all subjects of one group
    temp = control; seq = 'subject';
    if strcmp(type, 'p'); temp = patient; seq = 'patient'; end
    fprintf('the %s actually being computed is %s\n', seq, temp(np).name)
    code_parti = upper(temp(np).code);                                      % relevant information for later in the next few lines
    bad_trials = temp(np).badtrials;
    
    %% Starts with different steps available (see General settings above)
    for steps = steps2apply
        switch (steps)
            case 1 % Filter data, output is struct with filtered data for {1} erp- and {2} tfr-analyses
                hpf     = [.1 3];                                           % high-pass filter frequency
                lpf     = 30;
                inputdir= fullfile(wdir, 'cleaned');
                outdir = fullfile(wdir, 'data_cleaned');                    % directory at which data will be saved
                if ~exist(outdir, 'dir'); mkdir(outdir); end
                
                filename_clean = ...
                    {strcat('datclean_', code_parti, '_WO.mat'), ...
                    strcat('datclean_', code_parti, '_ALC.mat')};
                filename_preproc = ...
                    {strcat('datpreproc_', code_parti, '_WO.mat'), ...
                    strcat('datpreproc_', code_parti, '_ALC.mat')};
                
                for c = 1:2 % loop through both conditions
                    if exist(strcat(outdir, filename_preproc{c}), 'file')   % the next few lines check if data is already present and skips further processing if so to avoid redundancy
                        fprintf('\npre-preprocessing for %s already done, continuing with next step ...\n', code_parti )
                        continue
                    else
                        filter_rawdata(fullfile(inputdir, filename_clean{c}), ...
                            filename_preproc{c}, hpf, lpf, ROOTDIR, outdir)
                    end
                end
                
            case 2 % cut data to relevant chunks and merge conditions                
                outdir = fullfile(wdir, 'data_merged');                     % directory at which data will be saved
                if ~exist(outdir, 'dir'); mkdir(outdir); end
                
                filename_clean = ...
                    {strcat('datpreproc_', code_parti,'_WO.mat'), ...
                    strcat('datpreproc_', code_parti, '_ALC.mat')};
                filename_save = fullfile(outdir, ...
                    sprintf('data_merged_%s_%s.mat', type_calc, code_parti));
                
                if exist(filename_save, 'file')  % the next few lines check if data is already present and skips further processing if so to avoid redundancy
                    fprintf('\nepoching for %s already done, continuing with next step ...\n', code_parti )
                    continue
                else
                    epoch_trials(filename_clean, filename_save, ...
                        code_parti, wdir, ROOTDIR, type_calc)
                end
                
            case 3 % remove badtrials (after identifying them if necessary)
                inputdir = fullfile(wdir, 'data_merged');                   % directory data will be loaded from
                outdir = fullfile(wdir, 'data_final');                      % directory at which data will be saved
                filename_final = sprintf('data_final_%s_%s.mat', ...
                    type_calc, code_parti);
                filename_epoched = sprintf('data_merged_%s_%s.mat', ...
                    type_calc, code_parti);
                
                if exist(fullfile(outdir, filename_final), 'file')     % the next few lines check if data is already present and skips further processing if so to avoid redundancy
                    fprintf('\nremoving badtrials for %s already done, continuing with next step ...\n', code_parti )
                    continue
                else
                    load(fullfile(inputdir, filename_epoched)); %#ok<LOAD>  % this line loads the cleaned data into workspace (for details see clean_data.m)
                    
                    idx_fst = [1, 101];                                     % trialinfo that marks the beginning of a new block
                    idx_fstrl = []; clear tmp_trials                        % initialise variable to fill with content during loop
                    for c = 1:numel(idx_fst) % loop through both "conditions"
                        idx_fstrl = [idx_fstrl; ...                         % creates a vector of "first trials"
                            find(data_merged.trialinfo == idx_fst(c))];
                    end
                    
                    % Remove bad trials and 91/31 trials
                    idx_bad = [];
                    start_signal = [91, 191];                               % per default, every recording starts with a combination of '91/31 trials'
                    for b = 1:numel(start_signal)
                        idx_bad = [idx_bad; ...
                            find(data_merged.trialinfo==start_signal(b))];
                        idx_bad = [idx_bad; ...
                            find(data_merged.trialinfo==start_signal(b))+1];% necessary as code: '31' is ambiguous (start code *and* code of certain trials in MCST)
                    end
                    idx_bad = [idx_bad; find(data_merged.trialinfo>1000)];  %#ok<*AGROW>
                    
                    cfg = [];
                    cfg.trials = setdiff(1:length(data_merged.trialinfo), ...
                        idx_bad);
                    data_merged = ft_selectdata(cfg, data_merged);
                    
                    % Remove bad trials after eyeballing (if necessary)
                    if isempty(bad_trials) % if bad trials not yet identified, this plots all trials in order to check 'visually'
                        [bc, bad_trials] = plot_relevant_trials(data_merged);
                        if strcmp(type, 'p')
                            patient(np).badtrials = bad_trials;
                            patient(np).badchannels{1} = ...
                            unique([patient(np).badchannels{1}, bc{1:end-1}]);
                        else
                            control(np).badtrials = bad_trials;
                            control(np).badchannels{1} = ...
                            unique([control(np).badchannels{1}, bc{1:end-1}]);
                        end
                        save(fullfile(ROOTDIR, 'data', 'patdat.mat'),...    % save new patdat.mat file with metadata
                            'control', 'patient', '-v7.3');
                    end
                    
                    cfg = [];
                    cfg.trials = setdiff(1:length(data_merged.trialinfo), ...
                        bad_trials);
                    data_final = ft_selectdata(cfg, data_merged);
                    data_final = ft_struct2single(data_final);
                    save(fullfile(outdir, filename_final), ...
                        'data_final', '-v7.3');                             % saves the cleaned EEG data to one file
                end
        end
   end
end

end

function [bc, bt] = plot_relevant_trials(data_merged)
%% This function is just a helper to cut data and plot everything according
% to the trials saved in the trialdef file
frsp = 5000/data_merged.fsample;

idx = [{21:27}, {121:127}];                           % code at which a new trial was initiated (from which ERPs will be estimated later)
idx_trl = []; clear tmp_trials                      % initialise variable to fill with content during loop
for c = 1:numel(idx)
    tmp_trials = ismember(data_merged.trialinfo, idx{c});
    idx_trl = [idx_trl; find(tmp_trials)]; clear tmp_trials
end

cfg = [];
cfg.trials = idx_trl;                                                       % select only the trials that are useful for detecting artifacts
data_plot = ft_selectdata(cfg,data_merged);
data_plot.sampleinfo(:,1) = [1:length(data_plot.sampleinfo)].'.*1000;       % as two different conditions were "merged", a new definition of
data_plot.sampleinfo(:,2) = data_plot.sampleinfo(:,1)+500;                  % sampleinfo is included here so that ft_databrowser doesn't crash

data_merged.sampleinfo(:,1) = [1:length(data_merged.sampleinfo)].'.*1000;       % as two different conditions were "merged", a new definition of
data_merged.sampleinfo(:,2) = data_merged.sampleinfo(:,1)+500;                  % sampleinfo is included here so that ft_databrowser doesn't crash


% Start plotting data stacked 'vertically'
cfg = [];
cfg.channel         = 1:42;%'EEG';                                          % only EEG channels are plotted, to obtain approximately equal number of channels, the first 42 are selected
cfg.viewmode        = 'vertical';                                           % display data vertically
cfg.yaxis           = [-1 1].*40;                                           % scale of y-axis (arbitrary)
cfg.preproc.bpfilter= 'yes';                                                % band-pass filter settings in the next two lines
cfg.preproc.bpfreq  = [4 80];
cfg.preproc.bsfilter= 'yes';                                                % band-stop filter settings in the next two lines (notch filter for 50Hz noise)
cfg.preproc.bsfreq  = [48 52];
cfg.preproc.demean  = 'yes';                                                % de-mean data
cfg.preproc.detrend = 'yes';                                                % detrend data
%cfg.blocksize       = 25;                                                   % no. of seconds to display
cfg.layout          = 'eeg1005.lay';                                        % specifies the layout file that should be used for plotting

cfg.datafile = data_plot;
ft_databrowser(cfg, data_plot);

cfg.datafile = data_merged;
cfg.blocksize = 2.5;
ft_databrowser(cfg, data_merged);

% Displayed/select channels with artefacts according to some metric
cfg = [];
cfg.method      = 'summary';
cfg.metric      = 'zvalue';
cfg.channel     = 'EEG';
cfg.keepchannel = 'nan';                                                    % replacing "bad channels" with nan makes it easier to idetify them later
cfg.keeptrial   = 'nan';                                                    % replacing "bad channels" with nan makes it easier to idetify them later
dummy           = ft_rejectvisual(cfg, data_merged);

% Select bad trials according to 'ft_rejectvisual routine and save them
bt = find(cell2mat(arrayfun(@(q) all(isnan(dummy.trial{q}(:))), ...
    1:numel(dummy.trial), 'Un', 0)));
bc_select = ...
    find(any(isnan(cat(2,dummy.trial{setdiff(1:numel(dummy.trial), bt)})),2));
bc = {dummy.label{bc_select}}; %#ok<FNDSB>

prompt  = sprintf('Please take a minute to verify the results;\nbad trials:\t%s,\nbad channels:\t%s', ...
    regexprep(num2str(bt),'\s+',','), strjoin(bc,', '));
waitfor(warndlg(prompt, 'Warning'));
keyboard
close all
end
