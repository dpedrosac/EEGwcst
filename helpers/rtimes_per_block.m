function rtimes_per_block(subj1, subj2, ROOTDIR)

%   This function estimates the response times blockwise, that is per trial

%   Copyright (C) October 2021
%   D. Pedrosa, University Hospital of Gießen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

%% General settings
event_dir       = fullfile(ROOTDIR, 'data', 'header_and_events');
leg             = {'ET-patients', 'CTRL-subjects'};
conds           = {'WO', 'ALC'};                                            % two different conditions available
trls2est        = [10, 20:25];                                              % different trials to estimate; thereby 10/20 is the code for the tone indicating a shifting error and 21:25 is the code for a right answer
lims_outliers   = [.3 8];                                                   % limits at which data is considered wrong/artifact

%% Start extracting response times according to trial of interest
for g = 1:2 % loop through groups (ET-patients, CTRL-subj)
    if g == 1; sjts = subj1; else; sjts = subj2; end                        % gets a list for either CTRL-subj. or ET-patients, depending on the variable (g)
    rt_tmp{1,g} = cell(1,2); %#ok<AGROW>                                    % pre-allocate space
    for c = 1:numel(conds) % loop through conditions
        rt_tmp{1,g}{c} = nan(numel(sjts), numel(trls2est));
        
        for s = 1:numel(sjts) % loop through subjects
            filename = fullfile(event_dir, strcat('events_', conds{c}, '_', ...
                strcat('S', num2str(sjts(s))), '.mat')); load(filename);    % the last three lines define a filename and load the event data into workspace
            
            for e = 1:numel(trls2est)
                idx = [];                                                   % sets the indices to an empty value to fill it later
                ev_tmp = strcat('S', {' '}, num2str(trls2est(e)));          % finds trials coded as respective option (cf. trls2est)
                idx = [idx; find(strcmp({events.value}, ev_tmp)).'];        % add the indices as indicated by the trial number
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
                
                rt_tmp{1,g}{1,c}(s,e) = nanmean(tmp(:,1));                  % save data
            end
        end
    end
end

%% Plot results of grouped data stlye according to:
%   https://github.com/umeshmohan/TufteStyle

% grouped_data = []; % not used anymore!
% for k = 1:numel(trls2est)
%     for g = 1:2
%         grouped_data = [ grouped_data,  rt_tmp{1,g}{1,1}(:,k)];
%     end
% end

figure(103); hold on;
p = figure_params_gen;                                                      % load general parameters for plots

% General settings for the plot
fx_plots = {@(x) nanmean(x), ...                                            % plots of mean and SEM
    @(x) (nanmean(x) - 1.96*nanstd(x)./sqrt(size(x,1))), ...
    @(x) (nanmean(x) + 1.96*nanstd(x)./sqrt(size(x,1)))};
offset      = .03;                                                          % length of whiskers
sb          = {[1,1.25], [2,2.25], [3,3.25], [4,4.25], ...
    [5,5.25], [6,6.25], [6.5, 6.75], 8,8.75};                               % values for ticks of x-axis                                                     % plot order for the data structure
x_ticks     = {'Shift2D', 'Stay1', 'Stay2', 'Stay3', 'Last'};
marker      = {{'^', ones(3,1)}, {'o', zeros(1,3)}};                        % two markers one for each group
idx_xaxis   = [2:5, 7];
idx_columns = {[1,2], [3,4], [5,6], [7,8], [9,10], [11,12], [13,14]};
ttest_values= nan(1,length(idx_xaxis));
h           = ttest_values ;

iter = 0; c = 2;
for k = idx_xaxis
    for g = 1:2
        iter = iter+1;
        dat_bar = rt_tmp{1,g}{1,c}(:,k);
        wh_low(iter) = plot([sb{k}(g)-offset sb{k}(g)+offset], ...          % plots upper and lower whiskers
            fx_plots{2}(dat_bar).*[1 1], 'Color', p.greys{2}, ...
            'LineWidth', p.lnsize(3));
        wh_high(iter) = plot([sb{k}(g)-offset sb{k}(g)+offset], ...
            fx_plots{3}(dat_bar).*[1 1], 'Color', p.greys{2}, ...
            'LineWidth', p.lnsize(3));
        ln(iter) = plot(sb{k}(g).*[1 1], ...                                % plots the line in the middle between the whiskers
            [fx_plots{2}(dat_bar) fx_plots{3}(dat_bar)], ...
            'Color', p.greys{2}, 'LineWidth', p.lnsize(3));
        b(iter) = scatter(sb{k}(g), nanmean(dat_bar), marker{g}{1}, ...
            'MarkerEdgeColor', zeros(1,3), ...
            'MarkerFaceColor', marker{g}{2}, ...
            'LineWidth', p.lnsize(3));
        b(iter).SizeData = 1000;
        
        if g == 2 % run statistics iff g == 2;
            [ttest_values(k), h(k)] = ranksum(rt_tmp{1,g-1}{1,c}(:,k), ...
                rt_tmp{1,g}{1,c}(:,k), 'method','exact');
            if ttest_values(k) < .05
                text(mean(sb{k}), 3.5, '*', 'FontSize', p.ftsize(3), ...
                    'HorizontalAlignment', 'center')
                plot([sb{k}], [1 1].*3.6, '-k')
            end
        end
    end
end

set(gca, 'FontName', p.ftname, 'FontSize', p.ftsize(3), ...
    'XTick', [2.125,3.125,4.125,5.125,6.675], ...
    'XTickLabel', x_ticks, ...
    'YTick', [8:8:40]./10);
ylim([0.8 4]); xlim([2 8]);
% right_boundary = get(gca, 'XLim'); left_boundary = get(gca, 'YLim');
annotation('textbox', [0.6, 0.7, 0.1, 0.1], 'String', ...
    {'Response times within block', 'between groups [in sec.]'}, ...
    'FontSize', p.ftsize(3)*2)

% text(right_boundary(2), left_boundary(2), sprintf('Response times within one block\n per group [in sec.]'), ...
%     'FontSize', p.ftsize(3), 'HorizontalAlignment', 'center')

[~, objh] =legend([b(1), b(length(idx_xaxis)+1)], leg, ...
    'FontSize', p.ftsize(3), 'Location', 'SouthEast', 'box', 'off');
objhl = findobj(objh, 'type', 'patch'); %// objects of legend of type line
set(objhl, 'Markersize', 25, 'LineWidth', p.lnsize(3)); %// set marker size as desired

% Tufte style according to (cf reference above):
axes_handle_original = gca;
figure_handle = get(axes_handle_original, 'Parent');
original_axes_position = get(axes_handle_original,'Position');
x_axis_position = original_axes_position - [0 0.01 0 0];
x_axis_position(4) = eps;
x_axis_handle = axes('Parent', figure_handle,...
    'Position', x_axis_position,...
    'YTick', [], 'YTickLabel', [],...
    'box', 'off', 'tickdir', 'out',...
    'XLabel', get(axes_handle_original, 'XLabel'),...
    'XTickMode', get(axes_handle_original, 'XTickMode'),...
    'XTick', get(axes_handle_original, 'XTick'),...
    'XTickLabelMode', get(axes_handle_original, 'XTickLabelMode'),...
    'XTickLabel', get(axes_handle_original, 'XTickLabel'),...
    'XTickLabelRotation', get(axes_handle_original, 'XTickLabelRotation'), ...
    'FontSize', get(axes_handle_original, 'FontSize'));
linkaxes([axes_handle_original x_axis_handle],'x');

y_axis_position = original_axes_position - [0.01 -0.05 0 0];
y_axis_position(3) = eps;
y_axis_handle = axes('Parent', figure_handle,...
    'Position', y_axis_position,...
    'XTick', [], 'XTickLabel', [],...
    'box', 'off', 'tickdir', 'in',...
    'YLabel', get(axes_handle_original, 'YLabel'),...
    'YTickMode', get(axes_handle_original, 'YTickMode'),...
    'YTick', get(axes_handle_original, 'YTick'),...
    'YTickLabelMode', get(axes_handle_original, 'YTickLabelMode'),...
    'YTickLabel', get(axes_handle_original, 'YTickLabel'),...
    'YTickLabelRotation', get(axes_handle_original, 'YTickLabelRotation'), ...
    'FontSize', get(axes_handle_original, 'FontSize'));
linkaxes([axes_handle_original y_axis_handle],'y');
set(axes_handle_original,'Visible','off');

set(y_axis_handle, 'Ylim', [0.8 4])
set(x_axis_handle, 'Xlim', [1.75 7])

%%
diff = 0;
switch diff
    case (1)
        figure(104); hold on;
        p = figure_params_gen;                                                      % load general parameters for plots
        
        % General settings for the plot
        fx_plots = {@(x) nanmean(x), ...                                            % plots of mean and SEM
            @(x) (nanmean(x) - 1.96*nanstd(x)./sqrt(size(x,1))), ...
            @(x) (nanmean(x) + 1.96*nanstd(x)./sqrt(size(x,1)))};
        offset      = .03;                                                          % length of whiskers
        sb          = {[1,1.25], [2,2.25], [3,3.25], [4,4.25], ...
            [5,5.25], [6,6.25], [6.5, 6.75], 8,8.75};                               % values for ticks of x-axis                                                     % plot order for the data structure
        x_ticks     = {'Shift2D', 'Stay1', 'Stay2', 'Stay3', 'Last'};
        marker      = {{'^', ones(3,1)}, {'o', zeros(1,3)}};                        % two markers one for each group
        idx_xaxis   = [2:5, 7];
        idx_columns = {[1,2], [3,4], [5,6], [7,8], [9,10], [11,12], [13,14]};
        ttest_values= nan(1,length(idx_xaxis));
        h           = ttest_values ;
        
        iter = 0; c = 2;
        for k = idx_xaxis
            for g = 1:2
                iter = iter+1;
                dat_bar = rt_tmp{1,g}{1,1}(:,k) - rt_tmp{1,g}{1,2}(:,k);
                wh_low(iter) = plot([sb{k}(g)-offset sb{k}(g)+offset], ...          % plots upper and lower whiskers
                    fx_plots{2}(dat_bar).*[1 1], 'Color', p.greys{2}, ...
                    'LineWidth', p.lnsize(3));
                wh_high(iter) = plot([sb{k}(g)-offset sb{k}(g)+offset], ...
                    fx_plots{3}(dat_bar).*[1 1], 'Color', p.greys{2}, ...
                    'LineWidth', p.lnsize(3));
                ln(iter) = plot(sb{k}(g).*[1 1], ...                                % plots the line in the middle between the whiskers
                    [fx_plots{2}(dat_bar) fx_plots{3}(dat_bar)], ...
                    'Color', p.greys{2}, 'LineWidth', p.lnsize(3));
                b(iter) = scatter(sb{k}(g), nanmean(dat_bar), marker{g}{1}, ...
                    'MarkerEdgeColor', zeros(1,3), ...
                    'MarkerFaceColor', marker{g}{2}, ...
                    'LineWidth', p.lnsize(3));
                b(iter).SizeData = 1000;
                
                if g == 2 % run statistics iff g == 2;
                    [ttest_values(k), h(k)] = ...
                        ranksum(rt_tmp{1,g-1}{1,1}(:,k) - rt_tmp{1,g-1}{1,2}(:,k), ...
                        rt_tmp{1,g}{1,1}(:,k) - rt_tmp{1,g}{1,2}(:,k));
                    if ttest_values(k) < .05
                        text(mean(sb{k}), 4.3, '*', 'FontSize', p.ftsize(3), ...
                            'HorizontalAlignment', 'center')
                        plot([sb{k}], [1 1].*4.2, '-k')
                    end
                end
            end
        end
        
        set(gca, 'FontName', p.ftname, 'FontSize', p.ftsize(3), ...
            'XTick', [2.125,3.125,4.125,5.125,6.675], ...
            'XTickLabel', x_ticks, ...
            'YTick', [5:5:45]./10);
        
        ylim([0.0 1.2]); xlim([2 8]);
        [~, objh] =legend([b(1), b(length(idx_xaxis)+1)], leg, ...
            'FontSize', p.ftsize(3), 'Location', 'SouthEast', 'box', 'off');
        objhl = findobj(objh, 'type', 'patch'); %// objects of legend of type line
        set(objhl, 'Markersize', 25, 'LineWidth', p.lnsize(3)); %// set marker size as desired
        
        % Tufte style according to (cf reference above):
        axes_handle_original = gca;
        figure_handle = get(axes_handle_original, 'Parent');
        original_axes_position = get(axes_handle_original,'Position');
        x_axis_position = original_axes_position - [0 0.01 0 0];
        x_axis_position(4) = eps;
        x_axis_handle = axes('Parent', figure_handle,...
            'Position', x_axis_position,...
            'YTick', [], 'YTickLabel', [],...
            'box', 'off', 'tickdir', 'out',...
            'XLabel', get(axes_handle_original, 'XLabel'),...
            'XTickMode', get(axes_handle_original, 'XTickMode'),...
            'XTick', get(axes_handle_original, 'XTick'),...
            'XTickLabelMode', get(axes_handle_original, 'XTickLabelMode'),...
            'XTickLabel', get(axes_handle_original, 'XTickLabel'),...
            'XTickLabelRotation', get(axes_handle_original, 'XTickLabelRotation'), ...
            'FontSize', get(axes_handle_original, 'FontSize'));
        linkaxes([axes_handle_original x_axis_handle],'x');
        
        y_axis_position = original_axes_position - [0.01 -0.05 0 0];
        y_axis_position(3) = eps;
        y_axis_handle = axes('Parent', figure_handle,...
            'Position', y_axis_position,...
            'XTick', [], 'XTickLabel', [],...
            'box', 'off', 'tickdir', 'in',...
            'YLabel', get(axes_handle_original, 'YLabel'),...
            'YTickMode', get(axes_handle_original, 'YTickMode'),...
            'YTick', get(axes_handle_original, 'YTick'),...
            'YTickLabelMode', get(axes_handle_original, 'YTickLabelMode'),...
            'YTickLabel', get(axes_handle_original, 'YTickLabel'),...
            'YTickLabelRotation', get(axes_handle_original, 'YTickLabelRotation'), ...
            'FontSize', get(axes_handle_original, 'FontSize'));
        linkaxes([axes_handle_original y_axis_handle],'y');
        set(axes_handle_original,'Visible','off');
        
        %set(y_axis_handle, 'Ylim', [0.5 4.5])
        set(x_axis_handle, 'Xlim', [1.75 7])
end
