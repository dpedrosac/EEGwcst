function read_data(subj, wdir, type)

%   This function loads the (raw) data and performs a resampling with the
%   fieldtrip toolbox

%   Copyright (C) December 2017
%   D. Pedrosa, University Hospital of Gieﬂen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.


loadgen = 'Z:\EEG_raw\';                                                    % folder, in which all subjects are saved

%% General settings
load([wdir '\wcst\patdat.mat']);                                            % this file loads the meta data
rspl_freq       = 200;                                                      % frequency at which the data will be resampled
% outdir_branch   = 'Z:\wcst\';                                               % drive with principal directory in which data will be saved
outdir_branch   = 'C:\Users\David\projekte\wcst';                                               % drive with principal directory in which data will be saved
prefix          = 'EEG-Tremorstudie_';                                      % prefix of BRAINVISION data
steps2apply     = 1:2;
if strcmp(type, 'p'); tolom = subj{2}; else; tolom = subj{1}; end           % selects whether (p) pateints or controls (c) are analysed (see (type))

for np = tolom                                                             % provides a list of subjects to be analysed
    if strcmp(type, 'p')
        disp(['the patient actually being computed is ' patient(np).name]); % this line provides the patient currently being analyzed in the command window
        code_participant    = upper(patient(np).code);
        file_suffix = {patient(np).file, patient(np).file2};                % Because of nonuniform nomencalture, the filename needs to be provided at this point for the WO ...
        % and the ALCOHOL condition, respectively (see metadata)
        nam = patient(np).name;
    else
        disp(['the subject actually being computed is ' control(np).name]); % this line provides the patient currently being analyzed in the command window
        code_participant= upper(control(np).code);
        file_suffix     = {control(np).file, control(np).file2};            % Because of nonuniform nomencalture, the filename needs to be provided at this point for the WO ...
        % and the ALCOHOL condition, respectively (see metadata)
        nam = control(np).name;
    end
    load_dir = strcat(loadgen, code_participant);                           % defines and changes the working directory to where data is stored
    savedir = strcat(outdir_branch, '\', code_participant);
    
    %% Starts with the different steps available (see general settings)
    for steps = steps2apply                                                 % for better control, this loop enables to simply run some of the analyses
        switch (steps)
            case 1 %% this step reads data and saves it a format compatible with fieldtrip
                if ~exist(savedir, 'dir'), mkdir(savedir); end                          % create directory, if not present already
                
                try cd(load_dir);
                catch
                    fprintf('inexisting folder for subj:  %s \n', nam);
                    continue
                end
                
                cond = {'WO', 'ALC'};                                                   % different conditions, later important for saving
                for c = 1:2                                                             % loop through differentr conditions, i.e. 1 = (WO), 2 = (ALC) condition
                    filename_final = ...
                        strcat('datrsp_', code_participant, '_', cond{c}, '.mat');   % filename under which data will be saved in the (outdir) directiry
                    if exist(strcat(savedir, '\', filename_final), 'file')                               % if filename has already been processed, subject is skipped to save time
                        if c == 2; fprintf('subj: %s has already been processed, next subj... \n', nam); end
                        continue
                    else
                        clear filename1 data_all data_rsp;                              % clears data that will be used from now on in loops
                        filename1 = strcat(prefix, code_participant, '_', ...           % filename, so that data may be loaded
                            file_suffix{c}, '.eeg');
                        
                        if ~exist(filename1, 'file')
                            fprintf('problem with data from %s, skipping. \n', nam);
                            %                 continue;
                            keyboard;                                 % sloppy fix, as some files have really odd names, this is introduced
                        end                                                             % at this point to manually adapt names
                        
                        % Reads data from original files taking specific channels into account
                        cfg = [];
                        cfg.dataset = filename1;                                        % reads data from filename as defined before
                        cfg.channel = {'all', '-Flex*'};                                % Extract all data except some of the EMG
                        data_all = ft_preprocessing(cfg);
                        
                        % Resamples data to lower frequency
                        cfg = [];
                        cfg.resamplefs  = rspl_freq;                                    % resampling frequency
                        cfg.detrend     = 'no';
                        data_rsp        = ft_resampledata(cfg, data_all); %#ok<NASGU>
                        clear data_all;                                                 % clears data_all to save space on workspace
                        
                        % Saves data to pre-specified folder
                        save(strcat(savedir, '\', filename_final), 'data_rsp', '-v7.3');
                    end
                end
                
            case (2) % the next part extracts the header information and creates both header and event files
                try cd(load_dir);                                           % changes directory to subject-specfic folder
                catch
                    fprintf('inexisting folder for subj:  %s \n', nam);
                    continue
                end
                
                cond = {'WO', 'ALC'};                                       % different conditions, later important for saving
                for c = 1:2   % loop through different conditions, i.e. 1 = (WO), 2 = (ALC) condition
                    filename_final = ...
                        strcat('header_', code_participant, '_', cond{c}, '.mat');   % filename under which data will be saved in the (outdir) directiry
                    if exist(strcat(savedir, '\', filename_final), 'file')                               % if filename has already been processed, subject is skipped to save time
                        if c == 2; fprintf('subj: %s has already been processed, headers already available ... \n', nam); end
                        continue
                    else
                        clear filename_raw;                                 % clears data that will be used from now on in loops
                        filename_raw = strcat(loadgen, code_participant, ...
                            '\', prefix, code_participant, '_', ...         % filename, so that data may be loaded from the original BRAINvision file
                            file_suffix{c}, '.eeg');
                        
                        if ~exist(filename_raw, 'file')                     %this part of the code serves to avoid breaks, in case the file has a different name
                            fprintf('problem with data from %s, skipping. \n', nam);
%                             continue;
                            keyboard;                                       % sloppy fix, as some files have really odd names, this is introduced
                        end                                                 % at this point to manually adapt names
                        
                        % Reads data from original files taking specific channels into account
                        filename_header = ...
                            strcat('header_', code_participant, '_', cond{c}, '.mat');
                        
                        hdr = ft_read_header(filename_raw); %#ok<NASGU>
                        save(strcat(savedir, '\', filename_header), 'hdr', '-v7.3');
                        
                        events = ft_read_event(filename_raw); %#ok<NASGU>
                        save(strcat(savedir, '\', 'events_', cond{c}, '_', code_participant, ...
                            '.mat'), 'events', '-v7.3');
                        
                        % start getting the data to later epoch the file
                        cfg                         = [];
                        cfg.datafile                = filename_raw;
                        cfg.headerfile              = ...
                            strcat(filename_raw(1:end-4), '.vhdr');               % this just defines the orignal header file
                        cfg.trialfun                = 'ft_trialfun_general';    % this is the default
                        cfg.trialdef.eventtype      = 'Stimulus';
                        cfg.trialdef.prestim        = .5;    % in seconds before time = 0;
                        cfg.trialdef.poststim       = 2;  % in seconds after time = 0;
                        trialdef = ft_definetrial(cfg); %#ok<NASGU>
                        save(strcat(savedir, '\', 'trialdef_', ...
                            code_participant,'_', cond{c}, '.mat'), 'trialdef', '-v7.3'); % saves trial data
                        
                        filename_excl = strcat(outdir_branch, '\events_', cond{c}, '.xlsx' );
                        table = output_excel(events, hdr);
                        writetable(table, filename_excl,'Sheet',code_participant)
                    end
                end
        end
    end
end