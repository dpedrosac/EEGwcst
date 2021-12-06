function [wrong_runs, complete_runs, set_loss] = select_runs(code, ROOTDIR)

%   This function selects the complete runs and the ones with some sort of
%   error; data is appended to the events_xx.mat file in the .\data
%   directory

%   Copyright (C) June 2021, modified July 2021
%   D. Pedrosa, Urs Kleinholdermann University Hospital of Gießen and Marburg
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
wrong_runs = nan(1,2);
complete_runs = nan(1,2);
set_loss = nan(1,2);
for c = 1:2 % loops through the two different conditions
    filename_events = fullfile(load_dir, strcat('events_', conds{c}, '_', ...
        upper(code), '.mat'));
    load(filename_events)   %#ok<LOAD> % loads the data
    
    %% Description: =========================================   %%
    %   According to Barceló et al. (1999/2003), only completed
    %   trials should be used for ERP analyses. Therefore it is
    %   necessary to mark trials in which some error occurred;
    
    incorrect = zeros(numel(events),1);                                     % create new value for structure
    S1 = find(strcmp(cellfun(@(x) strrep(x,' ', ''), ...
        {events(:).value}, 'Uniform', 0), 'S1'));                           % gets the first trial of every repetition
    temp_set_loss = nan(numel(S1), numel(S1));
    iter = 0;
    for m = 1:numel(S1)                                                     % run though one run andlook for errors
        if m < numel(S1)
            events_temp = {events(S1(m):S1(m+1)).value};
        else
            events_temp = {events(S1(m):end).value};
        end
        error = find(strcmp(cellfun(@(x) strrep(x,' ', ''), ...
            events_temp, 'Uniform', 0), 'S40') | ...
            strcmp(cellfun(@(x) strrep(x,' ', ''), ...
            events_temp, 'Uniform', 0), 'S50'), 1);
        if ~isempty(error)
            iter = iter + 1;
            if nansum(strcmp(cellfun(@(x) strrep(x,' ', ''), ...
                    events_temp, 'Uniform', 0), 'S40')) > 0
                temp_set_loss(m, iter) = 0;
            elseif nansum(strcmp(cellfun(@(x) strrep(x,' ', ''), ...
                    events_temp, 'Uniform', 0), 'S40')) == 0 && ...
                    all(find(strcmp(cellfun(@(x) strrep(x,' ', ''), ...
                    events_temp, 'Uniform', 0), 'S50')) > 12)
                temp_set_loss(m, iter) = 1;
            end
            if m < numel(S1)
                incorrect(S1(m):S1(m+1)) = iter;
            else
                incorrect(S1(m):end) = iter;
            end
        end
        temp = num2cell(incorrect);
        [events.incomplete] = temp{:};
    end
    save(filename_events, 'events', '-v7.3') % saves updated data
    wrong_runs(1,c) = max(incorrect);
    complete_runs(1,c) = numel(S1);
    set_loss(1,c) = nansum(temp_set_loss(:));
end