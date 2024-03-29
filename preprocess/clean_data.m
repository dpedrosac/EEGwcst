function clean_data(subj, ROOTDIR, wdir, type)

%   This function loads the data and performs all 'cleaning' steps with
%   the fieldtrip toolbox, i.e. resampling, correction of bad channels,
%   etc.

%   Copyright (C) January 2018, modified August 2021
%   D. Pedrosa, University Hospital of Gie�en and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.


%% General settings
load(fullfile(wdir, 'patdat.mat'));                                         % this file loads the meta data
if strcmp(type, 'p'); tolom = subj{2}; else; tolom = subj{1}; end           % selects whether (p) pateints or controls (c) are analysed (see (type))

steps2apply     = 1:3;                                                        % four steps available: (1): load data (2): downsample (3): read HDR and events (4): artifact removal
flag_check      = 0;                                                        % this option is intended for debgugging purposes, that is to visualise steps
cond            = {'WO', 'ALC'};                                            % different conditions, later important for saving
fx_transpose    = @(x) x.';

for np = tolom                                                             % provides a list of subjects to be analysed
    temp = control; seq = 'subject';
    if strcmp(type, 'p'); temp = patient; seq = 'patient'; end
    fprintf("\n========================")
    fprintf('the %s actually being computed is %s\n', seq, temp(np).name)
    code_parti = upper(temp(np).code);                                      % relevant information for later in the next few lines
    bad_channels = temp(np).badchannels;
    bad_trials = temp(np).badtrials;
    
    %% Starts with the different steps available (see general settings)
    for steps = steps2apply                                                 % for better control, this loop enables to simply run some of the analyses
        switch (steps)
            case 1 %% removes components corresponding to blink artifacts identified from ICA
                
                inputdir= fullfile(wdir, 'rawEEG');
                
                % Define and create necessary output directories (outdir)
                outdir = fullfile(wdir, 'data_noica');                      % directory at which data will be saved
                if ~exist(outdir, 'dir'), mkdir(outdir); end                % create directory, if not present already
                
                filename_rspl = ...
                    {strcat('datrsp_', code_parti, '_WO.mat'), ...
                    strcat('datrsp_', code_parti, '_ALC.mat')};
                
                for c = 1:numel(cond) % loop through both conditions (WO) and (ALC)
                    filename_noica = ...
                        {strcat('datnoica_', code_parti, '_WO.mat'), ...
                        strcat('datnoica_', code_parti, '_ALC.mat')};
                    
                    if exist(fullfile(outdir, filename_noica{c}), 'file')  % the next few lines check if data is already present and skips further processing if so to avoid redundancy
                        fprintf('\nartifact interpolation for %s already done, continuing with next step ... \n', code_parti)
                        continue
                    else
                        load(fullfile(inputdir, filename_rspl{c}));         % this line loads the resampled data into workspace
                        data_rsp = sorted_data(data_rsp, 1);                % to make FT recognize 'Iz', it is necessary to rename it; besides, elec labels are sorted alphabetically for consistency
                        
                        % Filter and detrend data to remove drifts
                        cfg = [];                                           % this steps intends to remve drifts from the continous data
                        cfg.channel = 'EEG';                                % for ICA being more precise
                        cfg.hpfilter= 'yes';
                        cfg.detrend = 'yes';
                        cfg.hpfreq  = 0.5;
                        data_rsp    = ft_preprocessing(cfg, data_rsp);
                        
                        % Estimates the components available in the data
                        cfg = [];                                           % this steps runs an independent component analysis on the continous data
                        cfg.method  = 'runica';                             % to determine the blink artifacts and later reject them
                        cfg.channel = 'EEG';                                % select only EEG channels
                        data_comp  = ft_componentanalysis(cfg, data_rsp);
                        
                        % The next few line plot the Fpx channels in order
                        % to compare with the ICA later (if needed)
                        cfg = [];
                        cfg.channel         = {'Fp1', 'Fp2', 'Fpz'};        % channels to be plotted
                        cfg.viewmode        = 'butterfly';                  % data is plotted either 'vertical' or as a 'butterfly'
                        cfg.preproc.bpfilter= 'yes';                        % filtering and pre-processing of the data
                        cfg.preproc.bpfreq  = [1 40];
                        cfg.preproc.demean  = 'yes';
                        cfg.preproc.detrend = 'yes';
                        cfg.blocksize       = 10;                           % no. of seconds to display
                        cfg.layout          = 'EEG1005.lay';                % specifies the layout file that should be used for plotting
                        ft_databrowser(cfg, data_rsp)
                        
                        % The next lines plot the components available in order to
                        % selects the ones which correpsond to blink artefacts
                        cfg = [];
                        cfg.viewmode = 'component';
                        cfg.blocksize= 10;                                  % no. of seconds to display
                        cfg.layout   = 'EEG1005.lay';                       % specifies the layout file that should be used for plotting
                        try ft_databrowser(cfg, data_comp); catch; end
                        
                        done = 0;
                        while done == 0
                            prompt  = sprintf('Enter ICA-component to remove. In case of several components,use semicoli "";"" in between');
                            dlgtitle= 'ICA-Components to remove';                                   % GUI title
                            dims    = [1, 60];                                                      % GUI size
                            definput= "1;2";                                                        % GUI default input
                            keyboard;
                            
                            answer = inputdlg(prompt,dlgtitle, dims, definput);
                            x = fx_transpose(str2double(split(answer{1}, ';')));% extract the components to remove
                            done = 1; %#ok<NASGU>
                            if length(x) > 4
                                answer = questdlg('>4 components selected. Continue?', ...
                                    'Warning: Many components to be removed', 'Yes','No','No');
                                if strcmp(answer, 'Yes'); done = 1; else; done = 0; end
                            else
                                done = 1;
                            end
                        end
                        
                        % Data now to be reconstructed, excluding selected components
                        cfg = [];
                        cfg.component   = x;
                        data_noica      = ...                                   % the next line removes the selected components
                            ft_rejectcomponent(cfg, data_comp, data_rsp);
                        
                        %  The next few lines plot the difference between
                        %  before the ICA removal and after
                        if flag_check == 1
                            cfg = [];
                            cfg.datafile        = data_rsp;                 % first file to be plotted
                            cfg.channel         = 'EEG';                    % channels to plot
                            cfg.viewmode        = 'vertical';
                            cfg.preproc.bpfilter= 'yes';                    % defines the preprocessing thta happens with the data, that is
                            cfg.preproc.bpfreq  = [.1 30];                   % band-pass filtering
                            cfg.preproc.demean  = 'yes';                    % de-meaning
                            cfg.preproc.detrend = 'yes';                    % de-trending
                            cfg.blocksize       = 10;                       % no. of seconds to display
                            cfg.ylim           = [-1 1].*8;                                           % scale of y-axis (arbitrary)
                            cfg.layout          = 'EEG1005.lay';            % specifies the layout file that should be used for plotting
                            ft_databrowser(cfg, data_rsp)                   % before bad channel interpolation
                            
                            cfg.datafile = data_noica;                      % second file to be plotted
                            ft_databrowser(cfg, data_noica);                % after all steps of interpolation/ICA and rejection
                            keyboard
                        end
                    end
                    close all
                    save(fullfile(outdir, filename_noica{c}), ...
                        'data_noica', '-v7.3');
                end
                
            case (2) %%   This step determines 'bad_channels' and saves ..
                %   them in the metadat-file, spline interpolates those
                %   data in second step and, finally, does preprocessing
                
                inputdir = fullfile(wdir, 'data_noica');
                
                % Define and create necessary output directories (outdir)
                outdir = fullfile(wdir, 'data_clean');                      % directory at which data will be saved
                if ~exist(outdir, 'dir'), mkdir(outdir); end                % create directory, if not present already
                
                filename_noica = ...
                    {strcat('datnoica_', code_parti, '_WO.mat'), ...
                    strcat('datnoica_', code_parti, '_ALC.mat')};
                filename_clean = ...
                    {strcat('datclean_', code_parti, '_WO.mat'), ...
                    strcat('datclean_', code_parti, '_ALC.mat')};
                
                for c = 1:numel(cond) % loop through both conditions (WO) and (ALC)
                    if exist(fullfile(outdir, filename_clean{c}), 'file')   % the next few lines check if data is already present and skips further processing if so to avoid redundancy
                        fprintf('\ncleaning data for %s already done, continuing with next step ...', code_parti)
                        continue
                    else
                        load(fullfile(inputdir, filename_noica{c}));        %#ok<*LOAD> % this line loads the resampled data into workspace
                        
                        try load(fullfile(ROOTDIR, 'neighbours.mat'));      % loads the neighbours definition, and creates such a definition in case there is none
                        catch
                            dist = 40;                                      % distance neighbouring information is taken from (see define_neighbours)
                            neighbours = define_neighbours(data_noica, ...
                                wdir, dist, flag_check);                    % runs the function in which definition of neighbours is
                        end
                        
                        data_noica.elec = ...                               % this block assigns the electrode positions according to BRAINVISION to the data
                            ft_read_sens(fullfile(ROOTDIR, 'brainamp2.sfp'));
                        data_noica.elec = ...
                            ft_convert_units(data_noica.elec, 'mm');
                        data_noica = sorted_data(data_noica, 0);            % to make FT recognize 'Iz', it is necessary to rename it; besides, elec labels are sorted alphabetically for consistency
                        
                        % Remove bad trials after eyeballing (if necessary)
                        if isempty(bad_channels) || isempty(bad_trials{c})  % if bad trials not yet identified, this plots all trials in order to check 'visually'
                            [bc, bt] = ...
                                plot_relevant_trials(data_noica, ...
                                cond{c}, code_parti, wdir, bad_channels{1});
                            if strcmp(type, 'p')
                                patient(np).badtrials{c} = bt;
                                patient(np).badchannels{1} = ...
                                    unique([patient(np).badchannels{1}, ...
                                    bc{1:end-1}]);
                            else
                                control(np).badtrials{c} = bt;
                                control(np).badchannels{1} = ...
                                    unique([control(np).badchannels{1}, ...
                                    bc{1:end-1}]);
                            end
                            save(fullfile(ROOTDIR, 'data', 'patdat.mat'),...% save new patdat.mat file with metadata
                                'control', 'patient', '-v7.3');
                        end
                    end
                    
                    bc_all = unique(union(temp(np).badchannels{:}));
                    
                    %% In the next part 'bad channels are interpolated using
                    % spherical spline method (Perrin et al. 1989)
                    cfg = [];                                               % Spherical spline interpolation according to
                    cfg.method      = 'spline';                             % Perrin et al. 1989 is used as method
                    cfg.neighbours  = neighbours;
                    cfg.trials      = 'all';
                    %                     cfg.order           = 3;
                    %                     cfg.lambda          = 1e-5;                             % regularisation parameter
                    data_noarti     = data_noica;
                    cfg.badchannel  = reshape(bc_all, numel(bc_all),1);
                    data_noarti     = ft_channelrepair(cfg, data_noarti);
                    data_noarti     = sorted_data(data_noarti,0);
                    
                    if flag_check == 1 % plots the difference of interpolated channels
                        cfg = [];
                        cfg.channel         = bc_all;                       % only EEG channels are plotted, to obtain approximately equal number of channels, the first 42 are selected
                        cfg.ylim            = [-1 1].*12;
                        cfg.viewmode        = 'vertical';                   % display data vertically
                        cfg.preproc.bpfilter= 'yes';                        % band-pass filter settings in the next two lines
                        cfg.preproc.bpfreq  = [4 80];
                        cfg.preproc.bsfilter= 'yes';                        % band-stop filter settings in the next two lines (notch filter for 50Hz noise)
                        cfg.preproc.bsfreq  = [48 52];
                        cfg.preproc.demean  = 'yes';                        % de-mean data
                        cfg.preproc.detrend = 'yes';                        % detrend data
                        cfg.blocksize       = 5;                            % no. of seconds to display
                        cfg.layout          = 'EEG1005.lay';                % specifies the layout file that should be used for plotting
                        cfg.datafile = data_noica;
                        ft_databrowser(cfg, data_noica);                    % after all steps of interpolation/ICA and rejection
                        
                        cfg.datafile = data_noarti;
                        ft_databrowser(cfg, data_noarti);                   % after all steps of interpolation/ICA and rejection
                        keyboard
                        close all
                    end
                    clear data_noica
                    
%                     %% Rereference and preprocess data
%                     cfg = [];
%                     cfg.reref           = 'yes';                            % the next few lines are intended to create average-referenced data
%                     cfg.refchannel      = 'all';                            % to average reference data, the refchannel is set to 'all'
%                     %cfg.implicitref     = 'FCz';
%                     data_clean = ft_preprocessing(cfg, data_noarti);
                    data_clean = data_noarti;

                    save(fullfile(outdir, filename_clean{c}), ...
                        'data_clean', '-v7.3');
                end
        end
    end
end
end

function [bc, bt] = plot_relevant_trials(data_rsp, cond, subj, wdir, bc)
%% This function is just a helper to cut data and plot everything according
% to the trials saved in the trialdef file

frsp = 5000/200;
load(fullfile(wdir, 'trialdef', ...
    strcat(sprintf('trialdef_%s_%s.mat', subj, cond))));                    % this line loads the trial definitions so that data may be cut into chunks to be processed

%% Start selecting trials and merge data
cfg = [];
cfg.trl         = trialdef.trl;
cfg.trl(:,1:2)  = round(cfg.trl(:,1:2)./frsp);                              % the next two lines intend to "resample" the trialdefinition in a
cfg.trl(:,3)    = round(cfg.trl(:,3)./frsp);                                % admittedly not very elegant way, but it works!

idx_low         = find(cfg.trl(:,1) - cfg.trl(:,2) == -499);
idx_high        = find(cfg.trl(:,1) - cfg.trl(:,2) == -501);
if ~isempty(idx_low); cfg.trl(idx_low,2) = cfg.trl(idx_low,2)+1; end        % the next lines prevent
if ~isempty(idx_high); cfg.trl(idx_high,2) = cfg.trl(idx_high,2)+1; end     % problems at the "edges"

cfg.minlength   = 'maxperlen';                                              % this ensures all resulting trials are equal length
data_plot = ft_redefinetrial(cfg, data_rsp);

%%
cfg = [];
cfg.channel         = 1:42;                                                 % only EEG channels are plotted, to obtain approximately equal number of channels, the first 42 are selected
cfg.viewmode        = 'vertical';                                           % display data vertically
cfg.yaxis           = [-1 1].*40;                                           % scale of y-axis (arbitrary)
cfg.preproc.bpfilter= 'yes';                                                % band-pass filter settings in the next two lines
cfg.preproc.bpfreq  = [4 80];
cfg.preproc.bsfilter= 'yes';                                                % band-stop filter settings in the next two lines (notch filter for 50Hz noise)
cfg.preproc.bsfreq  = [48 52];
cfg.preproc.demean  = 'yes';                                                % de-mean data
cfg.preproc.detrend = 'yes';                                                % detrend data
cfg.layout          = 'EEG1005.lay';                                        % specifies the layout file that should be used for plotting

cfg.datafile = data_plot;
ft_databrowser(cfg, data_plot);

%% Displayed/select channels with artefacts according to some metric

cfg = [];
cfg.method      = 'summary';
cfg.metric      = 'zvalue';
cfg.channel     = data_plot.label(~ismember(data_plot.label, bc));
cfg.keepchannel = 'nan';                                                    % replacing "bad channels" with nan makes it easier to idetify them later
cfg.keeptrial   = 'nan';                                                    % replacing "bad channels" with nan makes it easier to idetify them later
% cfg.trials = setdiff(1:numel(data_plot.trial), bad_trials{c})             % possible option that could be implemented 
dummy           = ft_rejectvisual(cfg, data_plot);

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
