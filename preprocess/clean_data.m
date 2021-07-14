function clean_data(subj, wdir, type)

%   This function loads the data and performs all 'cleaning' steps with
%   the fieldtrip toolbox, i.e. resampling, correction of bad channels,
%   etc.

%   Copyright (C) January 2018, modified from a script of a different project
%   D. Pedrosa, University Hospital of Gieﬂen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

cd(wdir);
load([wdir 'wcst\patdat.mat']);                                             % this file loads the meta data

if strcmp(type, 'p'); tolom = subj{2}; else; tolom = subj{1}; end           % selects whether (p) pateints or controls (c) are analysed (see (type))

% General settings
steps2apply = 2;                                                            % four steps available: (1): load data (2): downsample (3): read HDR and events (4): artifact removal
flag_check  = 1;                                                            % this option is intended for debgugging purposes, that is to visualise steps
loaddir = [wdir 'wcst\'];
% loaddir = 'Z:\wcst\';
cond = {'WO', 'ALC'};

for np = tolom                                                             % provides a list of subjects to be analysed
    if strcmp(type, 'p')
        disp(['the patient actually being computed is ' patient(np).name]); % this line provides the patient currently being analyzed in the command window
        code_participant    = upper(patient(np).code);
    else
        disp(['the subject actually being computed is ' control(np).name]); % this line provides the patient currently being analyzed in the command window
        code_participant= upper(control(np).code);
    end
    outdir = strcat(loaddir, code_participant, '\');                        % directory at which data will be saved
    
    %% Starts with the different steps available (see general settings)
    for steps = steps2apply                                                 % for better control, this loop enables to simply run some of the analyses
        switch (steps)
            case 1 %% this step removes the components corresponding to blink artifacts identified from ICA
                
                % Defines rge filename from where data is loaded and
                % where/what filename to save
                filename_rspl = {strcat('datrsp_', code_participant, '_WO.mat'), ...
                    strcat('datrsp_', code_participant, '_ALC.mat')};
                
                filename_noica = {strcat('datnoica_', code_participant, '_WO.mat'), ...
                    strcat('datnoica_', code_participant, '_ALC.mat')};
                
                for c = 1:2 % loop through both conditions (WO) and (ALC)
                    if ~exist(strcat(outdir, filename_noica{c}), 'file')     % the next few lines check if data is already present and skips further processing if so to avoid redundancy
                        disp(['artifact interpolation for ' code_participant ' already ' ...
                            'done, continuing with next step ...'])
                        continue
                    else
                        load(strcat(outdir, filename_rspl{c}));             % this line loads the resampled data into workspace
                        
                        try                                                 % in order to make FT recognize 'Iz', it is necessary to rename it; besides, elec labels are sorted alphabetically for consistency
                            data_rsp.label{strcmp(data_rsp.label, 'IZ')} = 'Iz';
                            [data_rsp.label,I] = sort(data_rsp.label);
                            data_rsp.trial{1,1} = data_rsp.trial{1,1}(I,:);
                        catch
                        end
                        
                        % Filter and detrend data to remove drifts
                        cfg = [];                                           % this steps intends to remve drifts from the continous data
                        cfg.channel = 'EEG';                                % for ICA being more precise
                        cfg.hpfilter= 'yes';
                        cfg.detrend = 'yes';
                        cfg.hpfreq  = 0.5;
                        data_rsp    = ft_preprocessing(cfg, data_rsp);
                        
                        % Estimates the different components available in the data
                        cfg = [];                                           % this steps runs an independent component analysis on the continous data
                        cfg.method  = 'runica';                             % to determine the blink artifacts and later reject them
                        cfg.channel = 'EEG';                                % select only EEG channels
                        data_comp  = ft_componentanalysis(cfg, data_rsp);
                        
                        % The next few line plot the Fpx channels in order to
                        % compare with the ICA later (if needed)
                        cfg = [];
                        cfg.channel         = {'Fp1', 'Fp2', 'Fpz'};        % channels to be plotted
                        cfg.viewmode        = 'butterfly';                  % data is plotted either 'vertical' or as a 'butterfly'
                        cfg.preproc.bpfilter= 'yes';                        % filtering and pre-processing of the data
                        cfg.preproc.bpfreq  = [1 40];
                        cfg.preproc.demean  = 'yes';
                        cfg.preproc.detrend = 'yes';
                        cfg.blocksize       = 10;                           % no. of seconds to display
                        cfg.layout          = 'eeg1005.lay';                % specifies the layout file that should be used for plotting
                        ft_databrowser(cfg, data_rsp)
                        
                        % The next lines plot the components available in order to
                        % selects the ones which correpsond to blink artefacts
                        cfg = [];
                        cfg.viewmode = 'component';
                        cfg.blocksize= 10;                                  % no. of seconds to display
                        cfg.layout   = 'eeg1005.lay';                       % specifies the layout file that should be used for plotting
                        try ft_databrowser(cfg, data_comp); catch; end
                        
                        % The original data is now be reconstructed, excluding those components
                        cfg = [];
                        prompt          = 'Please note the components to be removed, e.g. [1,4,65]? ';
                        x               = input(prompt);
                        cfg.component   = x;
                        keyboard
                        data_noica   = ...                                   % the next line removes the selected components
                            ft_rejectcomponent(cfg, data_comp, data_rsp);
                        
                        %  The next few lines plot the difference between
                        %  before the ICA removal and after
                        if flag_check == 1
                            cfg = [];
                            cfg.datafile        = data_rsp;                 % first file to be plotted
                            cfg.channel         = 'EEG';                    % channels to plot
                            cfg.viewmode        = 'vertical';
                            cfg.preproc.bpfilter= 'yes';                    % defines the preprocessing thta happens with the data, that is
                            cfg.preproc.bpfreq  = [1 40];                   % band-pass filtering
                            cfg.preproc.demean  = 'yes';                    % de-meaning
                            cfg.preproc.detrend = 'yes';                    % de-trending
                            cfg.blocksize       = 10;                       % no. of seconds to display
                            cfg.layout          = 'eeg1005.lay';            % specifies the layout file that should be used for plotting
                            ft_databrowser(cfg, data_rsp)                   % before bad channel interpolation
                            
                            cfg.datafile = data_noica;                      % second file to be plotted
                            ft_databrowser(cfg, data_noica);                % after all steps of interpolation/ICA and rejection
                        end
                    end
                    save(strcat(outdir, filename_noica{c}), 'data_noica', '-v7.3');
                end
                
            case (2) %%   This step determines the 'bad_channels' and saves ..
                %   them in the metadat-file, spline interpolates those
                %   data in a second step and, finally, does the
                %   preprocessing
                
                % Defines filename from where data is loaded and
                % where/what filename to save
                filename_noica = {strcat('datnoica_', code_participant, '_WO.mat'), ...
                    strcat('datnoica_', code_participant, '_ALC.mat')};
                
                filename_clean = {strcat('datclean_', code_participant, '_WO.mat'), ...
                    strcat('datclean_', code_participant, '_ALC.mat')};
                
                for c = 1:2 % loop through both conditions (WO) and (ALC)
                    if ~exist(strcat(outdir, filename_clean{c}), 'file')     % the next few lines check if data is already present and skips further processing if so to avoid redundancy
                        disp(['cleaning data for ' code_participant ' already ' ...
                            'done, continuing with next step ...'])
                        continue
                    else
                        load(strcat(outdir, filename_noica{c}));             % this line loads the resampled data into workspace
                        
                        try
                            load(strcat(loaddir, 'neighbours.mat'));        % loads the neighbours definition, and creates such a definition in case there is none
                        catch
                            dist = 40;                                      % distance neighbouring information is taken from (see define_neighbours)
                            neighbours = define_neighbours(data_noica, ...
                                loaddir, dist, flag_check);                 % runs the function in which definition of neighbours is
                        end
                        
                        data_noica.elec = ...                               % this block assigns the electrode positions according to BRAINVISION to the data
                            ft_read_sens([loaddir 'brainamp2.sfp']);
                        data_noica.elec = ...
                            ft_convert_units(data_noica.elec, 'mm');
                        
                        [data_noica.label,I] = sort(data_noica.label);      % sorts the channels for consistency
                        data_noica.trial{1,1} = data_noica.trial{1,1}(I,:);
                        
                        if (strcmp(type, 'p') && isempty(patient(np).badchannels{c})) ...
                                || (strcmp(type, 'c') && isempty(control(np).badchannels{c}))
                            % if the 'badchannels are not available in the
                            % metadata, this part of the script is run
                            
                            cfg = [];
                            cfg.channel         = 1:42;%'EEG';              % only EEG channels are plotted, to obtain approximately equal number of channels, the first 42 are selected
                            cfg.viewmode        = 'vertical';               % display data vertically
                            cfg.preproc.bpfilter= 'yes';                    % band-pass filter settings in the next two lines
                            cfg.preproc.bpfreq  = [1 80];
                            cfg.preproc.bsfilter= 'yes';                    % band-stop filter settings in the next two lines (notch filter for 50Hz noise)
                            cfg.preproc.bsfreq  = [48 52];
                            cfg.preproc.demean  = 'yes';                    % de-mean data
                            cfg.preproc.detrend = 'yes';                    % detrend data
                            cfg.blocksize       = 25;                       % no. of seconds to display
                            cfg.layout          = 'eeg1005.lay';            % specifies the layout file that should be used for plotting
                            cfg.datafile = data_noica;
                            ft_databrowser(cfg, data_noica);                  % after all steps of interpolation/ICA and rejection
                            
                            % In the next few lines, the channels with
                            % artifacts are selected, according to the
                            % z-value of the data
                            cfg = [];
                            cfg.method = 'summary';
                            cfg.metric = 'zvalue';
                            cfg.channel = 'EEG';
                            cfg.keepchannel = 'yes';
                            dummy = ...
                                ft_rejectvisual(cfg, data_noica);
                            
                            % badchannels are saved in the metafile of the
                            % participants
                            prompt  = 'Please note the badchannels after visual inspection, e.g. {''Fp1'', ''T7'', ''T4''} ';
                            if strcmp(type, 'p')
                                patient(np).badchannels{c} = input(prompt);
                            else
                                control(np).badchannels{c} = input(prompt);
                            end
                            clear dummy
                            
                            save(strcat(loaddir, 'patdat.mat'), 'control', 'patient', '-v7.3');
                        end
                    end
                    
                    if strcmp(type, 'p')                                    % this creates a unique list of 'badchannels' per subject, that is
                        %                         bc_all = unique(cat(2,patient(np).badchannels{:}));
                        bc_all = union(patient(np).badchannels{:});         % 'bad channels'are merged for both conditions to aoid bias.
                    else
                        %                         bc_all = unique(cat(2,control(np).badchannels{:}));
                        bc_all = union(control(np).badchannels{:});
                    end
                    
                    % the next few lines create an index of badchannels and
                    % their locations (not needed necessarily)
                    fx_find = @(x) cellfun(@(c) strcmp(c, bc_all{x}), ...
                        data_noica.label, 'Un', 0);
                    idx_all = [];                               %#ok<NASGU>
                    idx_all = arrayfun(@(i) fx_find(i), ...
                        1:numel(bc_all), 'Un', 0); idx_all = ...
                        cell2mat(horzcat(idx_all{:}));
                    [idx,~] = find(idx_all==1); clear idx_all;
                    
                    %% In the next part, the 'bad channels are interpolated
                    % using the sopherical spline method as suggested by
                    % Perrin et al. 1989
                    cfg = [];                                               % Spherical spline interpolation according to
                    cfg.method          = 'spline';                         % Perrin et al. 1989 is used as method
                    cfg.neighbours      = neighbours;
                    cfg.trials          = 'all';
                    %                     cfg.order           = 3;
                    %                     cfg.lambda          = 1e-5;                             % regularisation parameter
                    data_noarti         = data_noica;
                    cfg.badchannel = reshape(bc_all, numel(bc_all),1);
                    data_noarti    = ...
                        ft_channelrepair(cfg, data_noarti);
                    [data_noarti.label,I]  = sort(data_noarti.label);       % sorts elec labels alphabetically for consistency
                    data_noarti.trial{1,1} = data_noarti.trial{1,1}(I,:);
                    
                    if flag_check == 1 % plots the difference of interpolated channels
                        cfg = [];
                        cfg.channel         = idx.';                        % only EEG channels are plotted, to obtain approximately equal number of channels, the first 42 are selected
                        cfg.viewmode        = 'vertical';                   % display data vertically
                        cfg.preproc.bpfilter= 'yes';                        % band-pass filter settings in the next two lines
                        cfg.preproc.bpfreq  = [1 80];
                        cfg.preproc.bsfilter= 'yes';                    % band-stop filter settings in the next two lines (notch filter for 50Hz noise)
                        cfg.preproc.bsfreq  = [48 52];
                        cfg.preproc.demean  = 'yes';                    % de-mean data
                        cfg.preproc.detrend = 'yes';                    % detrend data
                        cfg.blocksize       = 25;                       % no. of seconds to display
                        cfg.layout          = 'eeg1005.lay';            % specifies the layout file that should be used for plotting
                        cfg.datafile = data_noica;
                        ft_databrowser(cfg, data_noica);                  % after all steps of interpolation/ICA and rejection
                        
                        cfg.datafile = data_noarti;
                        ft_databrowser(cfg, data_noarti);                  % after all steps of interpolation/ICA and rejection
                    end
                    clear data_noica
                    
                    %% Rereference and preprocess data
                    cfg = [];
                    cfg.reref           = 'yes';                            % the next few lines are intended to create average-referenced data
                    cfg.refchannel      = 'all';                            % to average reference data, the refchannel is set to 'all'
                    cfg.implicitref     = 'FCz';
                    data_reref = ft_preprocessing(cfg, data_noarti);
                    
                    %% Add the peripheral data
                    load(strcat(loaddir, '\', code_participant, '\', ...
                        'datper_', code_participant, '_', cond{c}, '.mat'));% loads the peripheral (accelerometer) data
                    
                    % Preprocess the peripheral data
                    cfg = [];                                               % the next few lines do the preprocessing of the
                    cfg.channel         = 'all';                            % peripheral data
                    cfg.preproc.detrend = 'yes';
                    cfg.preproc.hpfilter= 'yes';
                    cfg.preproc.hpfreq  = 1;
                    data_peripheral     = ft_preprocessing(cfg, data_per);
                    data_peripheral = rmfield(data_peripheral, {'topo', 'unmixing', 'topolabel'});
                    
                    % Append data together
                    data_clean = ft_appenddata([],data_peripheral, ...
                        data_reref);
                    
                    save(strcat(outdir, filename_clean{c}), ...
                        'data_clean', '-v7.3');
                end
                
            case (3) % In this step, accelerometer data is separated and preprocessed (PCA)
                
                % Defines the filename from where data is loaded and
                % where/what filename to save
                filename_rspl = {strcat('datrsp_', code_participant, ...
                    '_WO.mat'), ...
                    strcat('datrsp_', code_participant, '_ALC.mat')};
                filename_peripheral = {strcat('datper_', code_participant, ...
                    '_WO.mat'), ...
                    strcat('datper_', code_participant, '_ALC.mat')};
                
                for c = 1:2 % loops through both conditions
                    if exist(strcat(outdir, 's', filename_peripheral{c}), ...
                        'file')
                        disp(['preprocessing the peripheral data for ' ...
                            code_participant ' already ' ...
                            'done, continuing with next step'])
                        clear filename*
                        continue
                    else
                        load(strcat(outdir, filename_rspl{c}));
                        
                        % Preprocess the Accelerometer data
                        cfg = [];                                           % this step reduces accelerometer data to principal component
                        cfg.method          = 'pca';                        % determine the principal component of the accelerometer data
                        cfg.channel         = ...
                            {'AccXLinks', 'AccYLinks', 'AccZLinks'};        % selectall three accelerometer channels
                        cfg.numcomponent    = 1;                            % 'boil down' to most prominent component
                        data_per              = ...
                            ft_componentanalysis(cfg, data_rsp);
                        data_per.label{1}     = 'acc';                      % renames accelerometer data to 'acc'
                        save(strcat(outdir, filename_peripheral{c}), ...    % saves data as data_per, which is later merged with the EEG data
                            'data_per', '-v7.3');
                    end
                end           
        end
    end
end
