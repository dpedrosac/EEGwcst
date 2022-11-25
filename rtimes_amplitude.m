function rtimes_amplitude(subj1proc, subj2proc, ROOTDIR, toi_time, bsl_time, ...
    channel, file_appendix, neighbours)

%   This function gathers response times and amplitude of ERP at a specific
%   toi per subject

%   Copyright (C) September 2022
%   D. Pedrosa, University Hospital of Gießen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

if nargin < 7
    fprintf('\nNot enough arguments, please provide at least 7 arguments')
    return
elseif nargin < 8
    neighbours = 'no';                                                      % only when explicitly called, the surrounding channels will be used as well 
end

%% General settings
event_dir       = fullfile(ROOTDIR, 'data', 'header_and_events');           % folder in which the events are saved
erp_dir         = fullfile(ROOTDIR, 'data', 'data_final');                  % folder with the final ERP data
leg             = {'ET-patients', 'CTRL-subjects'};                         % the two groups analyzed
conds           = {'WO', 'ALC'};                                            % two different conditions available
trls2est        = [21:25];                                                  % different trials to estimate; thereby 21:25 is the code for a right answer
lims_outliers   = [.3 8];                                                   % limits at which reaction times are considered wrong/artifact
fx_concat       = @(x) x(:);
fx_bslsub_all   = ...
    @(x,z) bsxfun(@minus, x, nanmean(x(:,z(1):z(2)),2));

%% Extract response times and ERP amplitude according to trial of interest

data_all = nan(2*numel([subj1proc, subj2proc]), 13); iter = 0;              % pre-allocate space
for g = 1:2 % loop through groups (ET-patients, CTRL-subj)
    if g == 1; sjts = subj1proc; else; sjts = subj2proc; end                % gets a list for either CTRL-subj. or ET-patients, depending on the variable (g)
    for c = 1:numel(conds) % loop through conditions
        for s = 1:numel(sjts) % loop through subjects
            iter = iter + 1;
            
            fprintf("Processing subj: %s (%d/%d)", ...
                strcat('S', num2str(sjts(s))), iter, ...
                2*numel([subj1proc, subj2proc]))
            
            filename = fullfile(event_dir, strcat('events_', conds{c}, '_', ... 
                strcat('S', num2str(sjts(s))), '.mat')); load(filename);% the last three lines define a filename and load the event data into workspace
            
            data_all(iter, 1) = sjts(s);                                    % subject number
            data_all(iter, 2) = g;                                          % group
            data_all(iter, 3) = c;                                          % condition
            
            for e = 1:numel(trls2est) % loop over trials of interest
                idx     = [];                                               % sets the indices to an empty value to fill it later
                ev_tmp  = strcat('S', {' '}, num2str(trls2est(e)));         % finds trials coded as respective option (cf. trls2est)
                
                % Start extracting reaction times
                idx     = [idx; find(strcmp({events.value}, ev_tmp)).'];    % add the indices as indicated by the trial number
                if strcmp(ev_tmp{1}, 'S 25'); idx = idx(1:end-1); end       % prevents error because of last trial
                tmp = [([events(idx+2).sample]./5000 - ...                  % time between display of the trial S1, S2, ..., S7 and the repsonse, when considering the previous answer
                    [events(idx+1).sample]./5000).', ...
                    ones(length(idx),1) *sjts(s)];
                
                % Remove outliers according to limits set before
                if max(tmp(:,1)) > lims_outliers(2) || ...                  % removes outliers if necessary
                        min(tmp(:,1)) < lims_outliers(1)
                    idx_outl = find(tmp(:,1) > lims_outliers(2) | ...       % provides the indices where response time is larger or smaller than
                        min(tmp(:,1) < lims_outliers(1)));                  % limits provided before and removes them from list in order to prevent outliers
                    idx_val = setdiff(1:size(tmp,1), idx_outl);
                    tmp = tmp(idx_val,:); clear idx_out1 idx_val;           % changes tmp, to a list wo/ outliers
                end
                
                data_all(iter, e+3) = nanmean(tmp(:,1));                    % mean reaction time for all trials {e} in subj {s} and condition {¢}
                 
                % Start extracting ERP data
                filename = fullfile(erp_dir, strcat('data_final_tfr_', ...
                    strcat('S', num2str(sjts(s))), '.mat')); load(filename);% the last lines define a filename and load the final data into worksapce for ERP ampitude extraction
                
                if iter == 1 % needs definition only once
                    cfg = []; % provides neighour data from which neighbouring channels may be selected
                    cfg.method          = 'distance';
                    cfg.feedback        = 'no';                             % can be set to 'yes' to see results
                    cfg.neighbourdist   = .03;
                    neighbours_data     = ...
                        ft_prepare_neighbours(cfg, data_final);
                end
                
                etemp = trls2est(e);                                        % trial of interest to extract
                if c == 2; etemp = trls2est(e) + 100; end                   % in case of ALC condition {c=2}, 100 mus be added to get correct trials
                idx_trial = data_final.trialinfo==etemp;                    % logical indexing of trials of interest
                bsl = dsearchn(data_final.time{1}', bsl_time');             % indices for baseline period as defined when calling function
                toi = dsearchn(data_final.time{1}', toi_time');             % indices for time of interest as defined when calling function
                
                switch neighbours
                    case ('no')
                        channel_of_interest = channel;
                    case ('yes')
                        chtemp = strcmp({neighbours_data(:).label}, channel);
                        channel_of_interest = ...
                            [channel, neighbours_data(chtemp).neighblabel]; % add the channels of interest to the original channel selected when calling function                        
                end
                
                cfg = []; % select data of interest with channel and trial defined before
                cfg.channel     = channel_of_interest;
                cfg.trials      = idx_trial;
                cfg.feedback    = 'no';
                dat_temp        = ft_selectdata(cfg, data_final);
                
                dat_tot = cat(1, dat_temp.trial{:});                        % concatenate everything to ch(1-x) x time matrix
                try
                    dat_tot = fx_bslsub_all(dat_tot, bsl); clear data_temp % subtract baseline from all trials
                    data_all(iter, e+8) = ...
                        nanmean(fx_concat(dat_tot(:, toi(1):toi(2))));      % estimates mean from time of interest and fills data table                    
                catch
                    data_all(iter, e+8) = NaN;
                end
            end
            clear dat_temp 
        end
    end
end

%% Save data to csv-file
data_table = table(data_all(:,1), data_all(:,2), data_all(:,3), ...
    data_all(:,4), data_all(:,5), data_all(:,6), data_all(:,7), ...
    data_all(:,8), data_all(:,9),data_all(:,10),data_all(:,11), ...
    data_all(:,12), data_all(:,13), ...
    'VariableNames', {'ID', 'group', 'condition', ...
    'rt21', 'rt22', 'rt23', 'rt24', 'rt25', ...
    'ERP21', 'ERP22', 'ERP23', 'ERP24','ERP25'});

writetable(data_table, fullfile(ROOTDIR, 'results', ...
    sprintf('ERPampl_reactionTime_%s_%s.csv', channel, file_appendix)))


