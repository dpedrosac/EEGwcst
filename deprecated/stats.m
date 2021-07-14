function stats(subj, wdir, type)

%   This function loads the events, sorts them and prints a summary of how
%   many different events are available

%   Copyright (C) January 2018
%   D. Pedrosa, University Hospital of Gieﬂen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.


%% General settings
load([wdir '\wcst\patdat.mat']);                                            % this file loads the meta data
outdir_branch   = 'C:\Users\David\projekte\wcst\';                                               % drive with principal directory in which data will be saved
steps2apply     = 2;
if strcmp(type, 'p'); tolom = subj{2}; else; tolom = subj{1}; end           % selects whether (p) pateints or controls (c) are analysed (see (type))
events_all{1} = {};
events_all{2} = {};

for np = tolom                                                             % provides a list of subjects to be analysed
    if strcmp(type, 'p')
        disp(['the patient actually being computed is ' patient(np).name]); % this line provides the patient currently being analyzed in the command window
        code_participant    = upper(patient(np).code);
        % and the ALCOHOL condition, respectively (see metadata)
        nam = patient(np).name;
    else
        disp(['the subject actually being computed is ' control(np).name]); % this line provides the patient currently being analyzed in the command window
        code_participant= upper(control(np).code);
        % and the ALCOHOL condition, respectively (see metadata)
        nam = control(np).name;
    end
    load_dir = strcat(outdir_branch, code_participant);                           % defines and changes the working directory to where data is stored
    savedir = strcat(outdir_branch, '\');
    
    %% Starts with the different steps available (see general settings)
    for steps = steps2apply                                                 % for better control, this loop enables to simply run some of the analyses
        switch (steps)
            case (1)
                
                
            case (2) %% this step reads data and saves it a format compatible with fieldtrip
                
                try cd(load_dir);
                catch
                    fprintf('inexisting folder for subj:  %s \n', nam);
                    continue
                end
                
                cond = {'WO', 'ALC'};                                                   % different conditions, later important for saving
                for c = 1:2                                                             % loop through differentr conditions, i.e. 1 = (WO), 2 = (ALC) condition
                    filename_final = ...
                        strcat('events_', cond{c}, '_', code_participant, '.mat');   % filename under which data will be saved in the (outdir) directiry
                    clear filename1 data_all data_rsp;                              % clears data that will be used from now on in loops
                    
                    if ~exist(filename_final, 'file')
                        fprintf('problem with data from %s, skipping. \n', nam);
                        continue;
                    end                                                             % at this point to manually adapt names
                    
                    load(filename_final)
                    [events(:).subj] = deal(sscanf(code_participant, 'S%d'));
                    events(end).difs = nan;
                    for k = 1:numel(events)-1
                    events(k).difs = [events(k+1).sample - events(k).sample]./5000;
                    end
                    events_all{c}{end+1} = events;
                end
        end
        
    end
end


isString            = cellfun('isclass', {events(:).value}, 'char');
isBUA(isString)     = ~cellfun(@isempty, regexp(reshape({events(isString).value}, numel({events(isString).value}),1),'BUA'));
[~, idx_tot] = find(isString == 1 & isBUA ==0);
ev_unique = unique(cellstr({events(idx_tot).value}));

datall = cat(2,events_all{2}{:});
counts_ev =  arrayfun(@(x) strcmp(ev_unique{x}, {datall.value}), 1:numel(ev_unique), 'Un', 0);
counts_evtot = cat(1,counts_ev{:}).';
numbers_tot = sum(counts_evtot,1)


counts_events = arrayfun(@(x) strcmp(ev_unique{x}, {events.value}), 1:numel(ev_unique), 'Un', 0)
counts_events = cat(1, ans{:});
numbers_tot = sum(counts_events,2);