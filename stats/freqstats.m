function freqstats(subj, ROOTDIR, wdir, type)

%   This function estuimates frequency stats for the TFR representation
%   of the final data and estimates the latter in cases where inexistant

%   ## Version 1.0 # first version

%   Copyright (C) December 2021
%   D. Pedrosa
%   University Hospital of Gießen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

%% General settings
[wdir, ROOTDIR] = EEGwcst_defaults(0);
addpath(fullfile(ROOTDIR, 'othercolor'))

%% Load necessary data into workspace
cd(wdir)
load(fullfile(wdir, 'patdat.mat'));                             %#ok<LOAD>  % this file loads the meta data
if strcmp(type, 'p'); tolom = subj{2}; else; tolom = subj{1}; end           % selects whether (p) pateints or controls (c) are analysed (see (type))
debug           = 0;                                            %#ok<NASGU>

% Define and create necessary input/output directories
outdir_tfr          = fullfile(wdir, 'tfr_data');                           % directory in which data will be saved after processing 'raw MRI'
if ~exist(outdir_tfr, 'dir'), mkdir(outdir_tfr); end                        % create directory, if not present already

% Start preparations for TFR statistics
for np = tolom
    temp = control; seq = 'subject';                                        % get metadata correctly along with next line
    if strcmp(type, 'p'); temp = patient; seq = 'patient'; end
    
    fprintf('the %s actually being computed is %s\n', seq, temp(np).name)
    code_parti = upper(temp(np).code);                                      % relevant information for later in the next few lines
    
    %% Create TFR representation if not existent
    filename_tfr_pre = fullfile(outdir_tfr, sprintf('tfr_%s', ...           % first part of filename which is completed later
        upper(code_parti)));
    if ~exist(sprintf('%s_memoryWO.mat', filename_tfr_pre), 'file')         % checks for existence of one of the files
        estimate_TFR(filename_tfr_pre, wdir, code_parti)
    end
    
end


end

function estimate_TFR(filename_tfr_pre, wdir, code_parti)

%   This function estimoates TF representations for preprocessed data using
%   multitaper methods

%   ## Version 1.0 # first version

%   Copyright (C) December 2021
%   D. Pedrosa
%   University Hospital of Gießen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

dir_preprocessed    = fullfile(wdir, 'data_final');                         % input directory
data_final = fullfile(dir_preprocessed, ...
    sprintf('data_final_tfr_%s.mat', code_parti));                          % loads preprocessed TFR files
load(data_final);                                                           % loads the final, preprocessed data into workspace

% compute time-freq representation TFR for every condition of interest
cfg = [];
cfg.output      = 'pow';
cfg.channel     = 'EEG';
cfg.method      = 'mtmconvol';
cfg.foi         = logspace(log10(5),log10(90),35); %3:1:90;
cfg.tapsmofrq   = 0.4 *cfg.foi;%(cfg.foi*2^(3/4/2) - cfg.foi*2^(-3/4/2))/1.5;
cfg.t_ftimwin   = 5./cfg.foi;%

cfg.toi         = -0.5:1/data_final.fsample:2;
cfg.taper       = 'dpss';
cfg.pad         = 'maxperlen'; %"nextpow2";
cfg.keeptrials = 'no';

conditions = {'memoryWO', 'efficientErrorWO', ...
    'memoryALC', 'efficientErrorALC'};

trls = {[21:25], [10,20], [110, 120], [121:125]};
for k = 1:numel(conditions)
    cfg.trials = find(ismember(data_final.trialinfo, trls{k}));
    tfreq      = ft_freqanalysis(cfg, data_final);
    
    save(sprintf('%s_%s.mat', filename_tfr_pre, conditions{k}), ...
        'tfreq', '-v7.3');
end

end