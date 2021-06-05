function rt_all = rt_gen(code, steps, paths)

%   This function extracts the data for the estimation of different possible
%   metrics related to the errors. Therefore, only the code of the subject
%   is needed in order to laod the file. Other possible options (steps) are:

%   (1) rt_all - finds a all the events and extracts the duration between
%   the presentation and the button press
%   (2) rt_time - provides a list with all response times for every trial
%   in the dataset
%   (3) error - provides a list with the number of errors and the
%   corresponding trial number to determine the amount of random/efficient
%   and total errors

%   Copyright (C) January 2018, modified June 2021
%   D. Pedrosa, Urs Kleinholdermann University Hospital of Gieﬂen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

events = cell(2,1);                                                         % pre-allocate space
load_dir = fullfile(paths.data_dir, 'header_and_events');
conds = {'WO', 'ALC'};
fx_tp = @(x) x.';                                                           % formula to transpose a vector
warning('off','MATLAB:strrep:InvalidInputType');                            % gets rid of the warning concerning the strrep function

for c = 1:2 % loops through the two different conditions
    load(fullfile(load_dir, strcat('events_', conds{c}, '_', ...
        upper(code), '.mat')))   % loads the data
    
    switch steps
        case 'rt_all' % all events
            idx_all = find(strcmp({events(:).value}, 'S  1') | ...          % finds the indices of the trial number; the event after
                strcmp({events(:).value}, 'S  2') | ...                     % the trial # is thereby the "error code", so that the difference
                strcmp({events(:).value}, 'S  3') | ...                     % between both constitutes the total response time
                strcmp({events(:).value}, 'S  4') | ...
                strcmp({events(:).value}, 'S  5') | ...
                strcmp({events(:).value}, 'S  6') | ...
                strcmp({events(:).value}, 'S  7'));             %#ok<NODEF>
            idx_all = idx_all(1:end-1);
            rt_all{c} = nanmean(([events(idx_all+1).sample] - ...
                [events(idx_all).sample])./5000);                           % estimates the difference between the presentation and the feedback, that is the response time
            
        case 'rt_time' % estimates the total time for the repetition trials
            idx_longest =  find(strcmp(cellfun(@(x) strrep(x,' ', ''), ...
                {events(:).value}, 'Uniform', 0), strcat('S25')));
            rt_time = nan(length(idx_longest),5);
            for k = 1:5 % loop through different responses S21-S25
                num = strcat('S2', num2str(k));
                idx =  find(strcmp(cellfun(@(x) strrep(x,' ', ''), ...      % finds the indices for the different responses
                    {events(:).value}, 'Uniform', 0), num));
                try
                    dattemp = fx_tp(([events(idx-1).sample] - ...           % estimates the time between the button press (S31-S34) and the start of the new trial (S1-S7)
                        [events(idx-2).sample])./5000);
                    rt_time(1:length(dattemp),k) = dattemp;
                catch
                    idx = idx(1:end-1);
                    dattemp = fx_tp(([events(idx-1).sample] - ...
                        [events(idx-2).sample])./5000);
                    rt_time(1:length(dattemp),k) = dattemp;
                end
            end
            rt_all{c} = rt_time;
            
        case 'error' % estimates the different errors available
            %% Description: =========================================   %%
            % In the paradigm described by Barcelo et al. (1999), different
            % errors are possible:
            %   -> Effective errors, when the rule changes and the
            %   participant tries to change the rule according to the
            %   feedback (coded as S10 and S20)
            %   -> Memory errors, when subjects "forget the rule", that is
            %   when the set is changed but someone mistakenly presses the
            %   wrong button but returns to the right rule afterwards
            %   (coded as S50 iff. it is not trial no. 3)
            %   -> perseverative errors, that is the rule is not updated
            %   and subjects maintain the wrong rule or randomly click in
            %   trial no. 3 (also S50).
            %   In order to discriminate all these, not only the error is
            %   saved but also the event before, that is the trial number
            
            indices = {'S10', 'S20', 'S40', 'S50'};                         % all potential errors (cf. mcst.tem and WCST-Tremorparadigma.sce in ~/paradigm/)
            idx = []; clear idx_temp;
            for k = 1:numel(indices) % loops through the different errors
                idx_temp =  strcmp(cellfun(@(x) strrep(x,' ', ''), ...
                    {events(:).value}, 'Uniform', 0), indices{k});
                idx = [idx; find(fx_tp(idx_temp))];                   %#ok<AGROW> % add the index list to the available indices
            end
            % Remove anticipations, that is when someone foresaw
            % (accidentally) the rule change and remove the entire set
            S40 = find(strcmp(cellfun(@(x) strrep(x,' ', ''), ...
                {events(:).value}, 'Uniform', 0), 'S40'));
            
            rem = 1; idx_bad = [];
            if rem == 1 && ~isempty(S40) % removes the entire block
                for m = 1:numel(S40) % loop through the events
                    try
                        lims1 = find(strcmp(cellfun(@(x) strrep(x,' ', ''), ...
                            {events(S40(m):S40(m)+30).value}, 'Uniform', 0), 'S1'));
                        lims1 = lims1-2;
                    catch
                        lims1 = find(strcmp(cellfun(@(x) strrep(x,' ', ''), ...
                            {events(S40(m):numel(events)).value}, 'Uniform', 0), 'S1'));
                        lims1 = lims1-2;
                        if isempty(lims1)
                            lims1 = numel(events)-S40(m);
                        end
                    end
                    lim = [S40(m)-2, S40(m)+lims1];
                    
                    try idx_bad = [idx_bad, lim(1):lim(2)]; catch; keyboard; end
                end
                idx_good = ~ismember(idx, idx_bad);
                idx = idx(idx_good);
            end
    end
    %% Incorporate information about surrounding events & remove parts
    errors = events(idx);
    rt_all{c} = rmfield(errors, {'type', 'duration', 'offset'});%#ok<AGROW> % remove this part as not needed and append some additional information
    for l = 1:length(idx)
        rt_all{c}(l).prev = fx_tp([events(idx(l)-2).value]);
        rt_all{c}(l).code = code;
        rt_all{c}(l).idx = idx(l);
    end
end