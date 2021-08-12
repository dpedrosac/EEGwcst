function read_data(subj, wdir, type)

%   This function loads the (raw) data and performs a resampling with the
%   fieldtrip toolbox

%   Copyright (C) December 2017, modified August 2021
%   D. Pedrosa, University Hospital of Gieﬂen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.


raw_folder = '/media/storage/EEGraw';                                       % folder, in which all subjects are saved (temporarily)

%% General settings
load(fullfile(wdir, 'patdat.mat'));                              %#ok<LOAD> % this file loads the meta data
rspl_freq       = 200;                                                      % frequency at which the data will be resampled
prefix          = 'EEG-Tremorstudie_';                                      % prefix of BRAINVISION data
steps2apply     = 1:2;
if strcmp(type, 'p'); tolom = subj{2}; else; tolom = subj{1}; end           % selects whether (p) pateints or controls (c) are analysed (see (type))
cond            = {'WO', 'ALC'};                                            % different conditions, later important for saving

% Start reading data one channel at a time
for np = tolom
    temp = control; seq = 'subject';
    if strcmp(type, 'p'); temp = patient; seq = 'patient'; end
    fprintf("\n========================")
    fprintf('the %s actually being computed is %s\n', seq, temp(np).name)
    code_participant = upper(temp(np).code);
    file_suffix = {temp(np).file, temp(np).file2};                          % Because of nonuniform nomencalture, the filename needs to be provided at this point for the WO ...
                                                                            % and the ALCOHOL condition, respectively (see metadata)
    nam = temp(np).name;
    read_dir = fullfile(raw_folder, code_participant);                      % defines directory to where data is stored
    
    %% Starts with the different steps available (see general settings)
    for steps = steps2apply                                                 % for better control, this loop enables to simply run some of the analyses
        switch (steps)
            case 1 %% this step reads data and saves it a format compatible with fieldtrip
                try cd(read_dir);
                catch
                    fprintf('inexisting folder for subj:  %s \n', nam);
                    continue
                end
                
                % Define and create necessary output directories (outdir)
                outdir = fullfile(wdir, 'rawEEG');                          % directory in which data will be saved after processing 'raw MRI'
                if ~exist(outdir, 'dir'), mkdir(outdir); end                % create directory, if not present already
                
                for c = 1:numel(cond)                                       % loop through differentr conditions, i.e. 1 = (WO), 2 = (ALC) condition
                    filename_final = ...
                        strcat('datrsp_', code_participant, '_', ...
                        cond{c}, '.mat');   % filename under which data will be saved in the (outdir) directiry
                    if exist(fullfile(outdir, filename_final), 'file')      % if filename has already been processed, subject is skipped to save time
                        if c == 2; fprintf('\n\tsubj: %s has already been processed, next subj... \n', nam); end
                        continue
                    else
                        clear filename1 data_all data_rsp;                              % clears data that will be used from now on in loops
                        filename1 = strcat(prefix, code_participant, '_', ...           % filename, so that data may be loaded
                            file_suffix{c}, '.eeg');
                        
                        if ~exist(fullfile(read_dir, filename1), 'file')
                            fprintf('\nproblem with data from subj %s @ %s condition please select file manually \n', nam, cond{c});
                            cd(read_dir)                                 % change directory and select file manually
                            [file,path] = uigetfile('*.eeg');
                            if isequal(file,0)
                                fprintf('No file selected, skipping to next subject ...\n');
                                continue
                            else
                                filename1 = file;
                                fprintf('The selected file for the %s condition is: %s\n', cond{c}, fullfile(path,file));
                            end
                            cd(wdir)
                        end
                        
                        % Reads data from original files taking specific channels into account
                        cfg = [];                                           % cfg is used in the for-loop for reading and downsampling data channelwise
                        cfg.resamplefs = rspl_freq;
                        %cfg.feedback    = 'text';
                        cfg.detrend     = 'no';
                        
                        singlechan = cell(1, 127);                          % pre-allocate space
                        for i = 1:127 % for loop necessary because of > 2GB size
                            cfg_temp = [];
                            cfg_temp.channel= i;                           % you can use a number as well as a string
                            cfg_temp.dataset= ...
                                fullfile(read_dir, filename1);              % reads data from filename as defined before
                            temp            = ft_preprocessing(cfg_temp);
                            fprintf('\nprocessing channel no. %s ', temp.label{1});
                            singlechan{i}   = ft_resampledata(cfg, temp);
                        end % for all channels
                        
                        cfg = []; % append data to one structure
                        data_rsp = ft_appenddata(cfg, singlechan{:});
                        
                        % Saves data to pre-specified folder
                        save(fullfile(outdir, filename_final), ...
                            'data_rsp', '-v7.3');
                    end
                end
                
            case (2) % extract header information & create header/event
                try cd(read_dir);                                           % changes directory to subject-specfic folder
                catch
                    fprintf('inexisting folder for subj:  %s \n', nam);
                    continue
                end
                
                % Define and create necessary output directories (outdir)
                outdir_header = fullfile(wdir, 'header_and_events');        % directory in which data will be saved after processing 'raw MRI'
                if ~exist(outdir_header , 'dir'), mkdir(outdir_header); end % create directory, if not present already
                
                for c = 1:numel(cond) % loop through conditions, i.e. 1 = (WO), 2 = (ALC)
                    filename_final = ...
                        strcat('header_', code_participant, '_', ...
                        cond{c}, '.mat');                                   % filename under which data will be saved in the (outdir) directiry
                    if exist(fullfile(outdir_header, ...
                            filename_final), 'file')                        % if filename has already been processed, subject is skipped to save time
                        if c == 2; fprintf('subj: %s has already been processed, headers already available ... \n', nam); end
                        continue
                    else
                        clear filename1 data_all data_rsp;                  % clears data that will be used from now on in loops
                        filename1 = strcat(prefix, code_participant, ...
                            '_', file_suffix{c}, '.eeg');                   % filename, so that data may be loaded
                        
                        if ~exist(fullfile(read_dir, filename1), 'file')
                            fprintf('\nproblem with data from subj %s @ %s condition please select file manually \n', nam, cond{c});
                            cd(read_dir)                                 % change directory and select file manually
                            [file,path] = uigetfile('*.eeg');
                            if isequal(file,0)
                                fprintf('No file selected, skipping to next subject ...\n');
                                continue
                            else
                                filename1 = file;
                                fprintf('The selected file for the %s condition is: %s\n', cond{c}, fullfile(path,file));
                            end
                        end
                        
                        % Reads header and preprocess data
                        filename_header = ...
                            strcat('header_', code_participant, '_', ...
                            cond{c}, '.mat');
                        
                        hdr = ft_read_header(fullfile(read_dir, filename1));
                        save(fullfile(outdir_header, filename_header), ...
                            'hdr', '-v7.3');
                        
                        events = ...
                            ft_read_event(fullfile(read_dir, filename1));
                        save(fullfile(outdir_header, strcat('events_', ...
                            cond{c}, '_', code_participant, ...
                            '.mat')), 'events', '-v7.3');
                        
                        % start getting the data to later epoch the file
                        cfg                         = [];
                        cfg.datafile                = ...
                            fullfile(read_dir, filename1);
                        cfg.headerfile              = ...
                            strcat(filename1(1:end-4), '.vhdr');            % this just defines the orignal header file
                        cfg.trialfun                = 'ft_trialfun_general';% this is the default
                        cfg.trialdef.eventtype      = 'Stimulus';
                        cfg.trialdef.prestim        = .5; % in seconds before time = 0;
                        cfg.trialdef.poststim       = 2;  % in seconds after time = 0;
                        trialdef = ft_definetrial(cfg);
                        save(fullfile(outdir_header, strcat('trialdef_', ...
                            code_participant,'_', cond{c}, '.mat')), ...
                            'trialdef', '-v7.3');                           % saves trial data
                        
                        filename_excl = fullfile(outdir_header, ...
                            strcat('events_', cond{c}, '.xlsx'));
                        table_output = output_excel(events, hdr);
                        writetable(table_output, filename_excl, ...
                            'Sheet',code_participant)
                    end
                end
        end
    end
end