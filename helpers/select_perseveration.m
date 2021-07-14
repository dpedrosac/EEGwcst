function perseverative_error = select_perseveration(code, ROOTDIR)

%   This function extracts the perseverative errors

%   Copyright (C) July 2021
%   D. Pedrosa, University Hospital of Gießen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

%% General settings
events = cell(2,1);                                                         % pre-allocate space
load_dir = fullfile(ROOTDIR, 'data', 'header_and_events');
conds = {'WO', 'ALC'};
warning('off','MATLAB:strrep:InvalidInputType');                            % gets rid of the warning concerning the strrep function
perseverative_error = nan(1,2);

for c = 1:2 % loops through the two different conditions
    filename_events = fullfile(load_dir, strcat('events_', conds{c}, '_', ...
        upper(code), '.mat'));
    load(filename_events)   %#ok<LOAD> % loads the data
    
    %% Description: =========================================   %%
    % According to Barceló et al. (1999/2003), trials without
    % category switch at the third trial are defined as
    % perseveration;
    
    incorrect = zeros(numel(events),1);                                     % create new value fror structure
    S3 = find(strcmp(cellfun(@(x) strrep(x,' ', ''), ...
        {events(:).value}, 'Uniform', 0), 'S3'));                           % selects/identifies the third trial of every repetition
    iter = 0;
    for m = 1:numel(S3)                                                     % run though one run andlook for errors
        if events(S3(m)).incomplete == 1; continue; end                     % skip runs identified as incoplete (cf. select_runs.m)
        events_temp = {events([S3(m)-4, S3(m)-1, S3(m), S3(m)+2]).value};   % identifies a specific sequencde with two inefficient errors at S1/2 (10, 20) and a perseveration (50) at S3
        
        try
            if strcmp(cell2mat(events_temp), 'S 10S 20S  3S 50')
                iter = iter + 1;
                incorrect(S3(m)) = iter;
            end
        catch
        end
        temp = num2cell(incorrect);
        [events.incomplete] = temp{:};     
    end
    perseverative_error(1,c) = max(incorrect);
end