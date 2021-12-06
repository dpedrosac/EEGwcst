function rt_subj = rtimes_trials(subj1, subj2, ROOTDIR)

%   This function estimates the response times after a repeat trial or a
%   "shift trial" plots the results and provides a table to run ANOVAE if
%   wanted

%   Copyright (C) March 2018, modified June 2021 and July 2021
%   D. Pedrosa, University Hospital of Gießen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

%% General settings
event_dir = fullfile(ROOTDIR, 'data', 'header_and_events');
conds = {'WO', 'ALC'};                                                      % two different conditions available
trls2est = {[10, 20], 21:25};                                               % different trials to estimate; thereby 10/20 is the code for the tone indicating a shifting error and 21:25 is the code for a right answer
lims_outliers = [.3 12.5];                                                    % limits at which data is considered wrong/artifact
flag_plot = 0;
%% Start estimating response times according to trial of interest
for g = 1:2 % loop through both groups (ET-patients, CTRL-subj)
    if g == 1; sjts = subj1; else; sjts = subj2; end                        % gets a list for either CTRL-subj. or ET-patients, depending on the variable (g)
    for k = 1:numel(trls2est); rt_tmp{1,k} = {[], []}; end                  % pre-allocate space, so that it may be filled with data later
    
    for s = 1:numel(sjts) % loop through subjects
        for c = 1:numel(conds) % loop through available conditions
            filename = fullfile(event_dir, strcat('events_', conds{c}, '_', ...
                strcat('S', num2str(sjts(s))), '.mat')); load(filename);    % the last three lines define a filename and load the event data into workspace
            
            for e = 1:numel(trls2est) % loop through different errors
                idx = [];                                                   % sets the indices to an empty value to fill it later
                for i = 1:numel(trls2est{e}) % loop through the distinct marked events to index for error (e)
                    ev_tmp = strcat('S', {' '}, num2str(trls2est{e}(i)));
                    idx = [idx; find(strcmp({events.value}, ev_tmp)).'];    % add the indices as indicated by the trial number
                end
                
                idx = idx(1:end-1); clear tmp;                              % the last index MUST be removed because only to response time after an answer is of interest
                tmp = [([events(idx+2).sample]./5000 - ...                  % time between display of the trial S1, S2, ..., S7 and the repsonse, when considering the previous answer
                    [events(idx+1).sample]./5000).', ones(length(idx),1) *sjts(s)];
                
                if max(tmp(:,1)) > lims_outliers(2) || ...                  % removes outliers if necessary
                        min(tmp(:,1)) < lims_outliers(1)
                    idx_outl = find(tmp(:,1) > lims_outliers(2) | ...       % provides the indices where response time is larger or smaller than
                        min(tmp(:,1) < lims_outliers(1)));                  % limits provided before and removes them from list in order to prevent outliers
                    idx_val = setdiff(1:size(tmp,1), idx_outl);
                    tmp = tmp(idx_val,:); clear idx_out1 idx_val;           % changes tmp, to a list wo/ outliers
                end
                
                rt_tmp{1,e}{c} = cat(1, [rt_tmp{1,e}{c}], tmp);             % concatenate data to a list with all response times
            end
        end
    end
    if g == 1; rt_et = rt_tmp; else; rt_ctrl = rt_tmp; end                  % saves data to either CTRL or ET structure.
end

%% Plot results
switch flag_plot
    case (1)
        figure(101); hold on;
        p = figure_params_gen;                                                      % load general parameters for plots
        fx_plots = {@(x) nanmean(x), ...                                            % plots of mean and SEM
            @(x) (nanmean(x) - 1.96*nanstd(x)./sqrt(size(x,1))), ...
            @(x) (nanmean(x) + 1.96*nanstd(x)./sqrt(size(x,1)))};
        fx_outliers = @(x) x(x > lims_outliers(1) & x < lims_outliers(2));
        
        % general settings for the plot
        offset  = .07;                                                               % length of the whiskers
        sb      = [1,1.5,3,3.5];                                                     % plot order for the data structure
        clrs    = {p.colors{4}, p.colors{1}};
        leg     = {'ET-patients', 'CTRL-subjects'};
        
        % Start plotting
        for c = 1:numel(conds) % loop through conditions and create subplots for both
            subplot(1,2,c); hold on; iter = 0;
            for k = 1:numel(trls2est) % loop through different trials
                for g = 1:2 % loop through both groups
                    iter = iter +1;
                    if g == 1; rt_all = rt_et; else; rt_all = rt_ctrl; end
                    dat{g} = fx_outliers(rt_all{1,k}{c}(:,1)); %#ok<*AGROW>
                    
                    sjts_idx = unique(rt_all{1,k}{c}(:,2));                         % number of subjects, needed to aggregate means of data
                    dat_bar = nan(numel(sjts_idx),1);                               % pre-allocate space in order to fill it in the next few lines
                    for s = 1:numel(sjts_idx) % loop through all subjects in this group
                        idx_s = find(rt_all{1,k}{c}(:,2) == sjts_idx(s));           % finds the indices for specific subj(s)
                        dat_bar(s,1) = nanmean(dat{g}(idx_s));                      % fills (dat_bar) with data
                    end
                    dat_ttest{iter} = dat_bar; %#ok<*NASGU>
                    b(iter) = bar(sb(iter), nanmean(dat_bar), ...                   % plots bars
                        'FaceColor', clrs{g}, 'FaceAlpha', .7);
                    wh_low(iter) = plot([sb(iter)-offset sb(iter)+offset], ...      % plots upper and lower whiskers
                        fx_plots{2}(dat_bar).*[1 1], 'Color', p.greys{2}, ...
                        'LineWidth', p.lnsize(3));
                    wh_high(iter) = plot([sb(iter)-offset sb(iter)+offset], ...
                        fx_plots{3}(dat_bar).*[1 1], 'Color', p.greys{2}, ...
                        'LineWidth', p.lnsize(3));
                    ln(iter) = plot(sb(iter).*[1 1], ...                            % plots the line in the middle between the whiskers
                        [fx_plots{2}(dat_bar) fx_plots{3}(dat_bar)], ...
                        'Color', p.greys{2}, 'LineWidth', p.lnsize(3));
                end
                
                set(gca, 'FontName', p.ftname, 'FontSize', p.ftsize(2), ...
                    'XTick', [mean(sb(1:2)), mean(sb(3:4))], 'XTickLabel', ...
                    {'set-shift trials', 'memory trials'});
                ylim([0.4 5]); xlim([.25 4.25]);
                
                if c == 1; ylabel(sprintf('Reaction time [in s.]')); end
            end
        end
        legend(b([1,2]), leg)                                                       % adds a legend to the plot
end
%% Put data together as a table
dat = []; iter = 0;
for g = 1:2 % loop through groups
    if g == 1; rtall = rt_et; else; rtall = rt_ctrl; end
    sjts = unique(rtall{1}{1,c}(:,2));
    for s = 1:numel(sjts)
        iter = iter +1;
        for c = 1:2 % loop through conditions
            for k = 1:2 % loop through errors
                idx = find(rtall{k}{1,c}(:,2) == sjts(s));
                tmp = [sjts(s), iter, g, nanmean(rtall{1,k}{1,c}(idx,1)), c, k];
                dat = [dat; tmp];
            end
        end
    end
end

tbl_anova = table(dat(:,1), dat(:,2),dat(:,3),dat(:,4),dat(:,5),dat(:,6));
tbl_anova.Properties.VariableNames = ...
    {'ID', 'subj', 'group', 'rt', 'cond', 'error'};
writetable(tbl_anova, fullfile(ROOTDIR, 'data', ...
    'anova_rtimes_error.txt'), 'Delimiter', '\t');

%% Create data for different groups
dattemp = {};
fx_outliers = @(x) x(x>lims_outliers(1) & x < lims_outliers(2));

for g = 1:2 % loop through both groups
    iter = 0;
    if g == 1; rt_all = rt_et; else; rt_all = rt_ctrl; end
    for k = 1:numel(trls2est) % loop through different trials
        for c = 1:2 % loop through conditions and create matrix in cell
            iter = iter +1;
            dattemp{g} = fx_outliers(rt_all{1,k}{c}(:,1));
            
            sjts_idx = unique(rt_all{1,k}{c}(:,2));                         % number of subjects, needed to aggregate means of data
            dat_bar = nan(numel(sjts_idx),1);                               % pre-allocate space in order to fill it in the next few lines
            for s = 1:numel(sjts_idx) % loop through all subjects in this group
                idx_s = find(rt_all{1,k}{c}(:,2) == sjts_idx(s));           % finds the indices for specific subj(s)
                dat_bar(s,1) = nanmean(dattemp{g}(idx_s));                  % fills (dat_bar) with data
            end
            rt_subj{g}(:,iter) = dat_bar;
        end
    end
end
