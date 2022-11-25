function results_2(control, patient, subj1, subj2, ROOTDIR, wdir)

%   Generate ERP between both groups with several available options

%   Copyright (C) July 2021
%   D. Pedrosa University Hospital of Gießen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

cd(ROOTDIR)
print_legend

%% Define groups
all_subj = [subj1, subj2];
excl_subj = [ 16];                                                     % excluded subjects due to technical reasons
subj1proc = subj1(~ismember(subj1, excl_subj));
subj2proc = subj2(~ismember(subj2, excl_subj));

idx_group       = {find(ismember(all_subj, subj1proc)), ...                     % get 'group indices'
    find(ismember(all_subj, subj2proc))};

%% Description of the analyses
% In what follows, different analyses aim at displaying the results of the
% study; first, the set-shifting differences between groups are displayed
% (figure 14) before and after intale of alcohol. Differences in the
% late p300. There is a cross-interaction in the effect of alcohol which is
% additionally plotted using R routines. Secondly, the "memory" condition
% are shown. For that a separate extraction is performed to get reaction
% times and ERP amplitudes at both conditions for both groups.


%% All shift trials at /wo condition (10,20) - cue-locked
if ~exist(fullfile(ROOTDIR, 'data', 'avg_shift_erp.mat'), 'file')
    estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [10, 20], 'shift')
end
load(fullfile(ROOTDIR, 'data', 'avg_shift_erp.mat'))

%% ET-patients vs control subjects shift trials ([10, 20]) - cue-locked
avg_shift1 = cell(1,2);
for k   = 1:2; avg_shift1{k} = avg(idx_group{k}); end                       % separate average in groups
ch      = {'Fz', 'Cz', 'Pz'};                                               % channels to plot
toi     = [-.4 1];                                                          % time of interest (x-axis)
bsl     = [-.15 0];                                                         % baseline to subtract

lgnd = {'ET-patients', 'CTRL-subjects'};
fignum = 14;
nCol = 3;
mcp = 'cluster_ft';
tit = 'ET-patients vs. CTRL-subj on shift trials';
plot_ERPcomparisions(avg_shift1, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)

%% All shift trials at /ALC condition (110,120) - cue locked
if ~exist(fullfile(ROOTDIR, 'data', 'avg_shiftALC_erp.mat'), 'file')
    estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [110, 120], 'shiftALC')
end
load(fullfile(ROOTDIR, 'data', 'avg_shiftALC_erp.mat'))

avg_shift2 = cell(1,2);
for k   = 1:2; avg_shift2{k} = avg(idx_group{k}); end                       % separate average in groups
clear avg
% Plotting not necessary, as this is performed in the R scripts (see below)

%% Estimate differences at specific time points and between conditions,
%   that is repeated measures ANOVA

%% TODO: avg_shift2  must be turned to one cell (as avg before!!)
plot_erp_amplitudes_boxplot(cat(2, avg_shift1{:}), cat(2, avg_shift2{:}), all_subj, ...
    {[1:19], [20:39]}, 'Pz', [.6 .8], 'p600-800_shift')

%% TODO: avg_shift2  must be turned to one cell (as avg before!!)
plot_erp_amplitudes_boxplot(cat(2, avg_shift1{:}), cat(2, avg_shift2{:}), all_subj, ...
    {[1:19], [20:39]}, 'Fz', [.25 .4], 'p300_shift')


%% Insert results from ANOVA (cross-over interaction) in the right upper corner
% of figure 14; for that, R code es available (./Rscripts/ERPcomparisons_Pz.r)


%% All memory trials at wo/ alcohol condition (21:25) - cue-locked
if ~exist(fullfile(ROOTDIR, 'data', 'avg_memory_erp.mat'), 'file')
    estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, ...
        [21:25], 'memory')
end
load(fullfile(ROOTDIR, 'data', 'avg_memory_erp.mat'))

%% ET-patients vs control subjects (all) memory trials
avg_memory=cell(1,2);
for k   = 1:2; avg_memory{k} = avg(idx_group{k}); end                   % separate average in groups
ch      = {'Fz', 'Cz', 'Pz'};                                               % channels to plot
toi     = [-.15 1];                                                         % time of interest (x-axis)
bsl     = [-.15 0];                                                         % baseline to subtract

lgnd = {'ET-patients', 'CTRL-subjects'};
fignum = 34;
nCol = 3;
mcp = 'cluster_ft';
tit = 'ET-patients vs. CTRL-subj on memory trials';
plot_ERPcomparisions(avg_memory, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)

%% Plot reaction time vs. amplitude in frontal P3b
%rtimes_amplitude(subj1proc, subj2proc, ROOTDIR, [.55 .75], [-.15 0], 'Pz', 'p3B', 'yes')

%% Insert results from ANOVA (cross-over interaction) in the right upper corner
% of figure 14; for that, R code es available (./Rscripts/ERPcomparisons_Pz.r)


%%%%%%%%%%%%%%%%%%%


%% All memory trials at /wo condition (2:7), cue-locked
if ~exist(fullfile(ROOTDIR, 'data', 'avg_mem-cuelocked_erp.mat'), 'file')
    estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [2:7], ...
        'mem-cuelocked', [10,20,110,120])
end
load(fullfile(ROOTDIR, 'data', 'avg_mem-cuelocked_erp.mat'))

%% ET-patients vs control subjects memory trials cue locked
avg_memoryCue=cell(1,2);
for k   = 1:2; avg_memoryCue{k} = avg(idx_group{k}); end                    % separate average in groups
ch      = {'Fz', 'Cz', 'Pz'};                                               % channels to plot
toi     = [-.15 1];                                                         % time of interest (x-axis)
bsl     = [-.15 0];                                                         % baseline to subtract

lgnd = {'ET-patients', 'CTRL-subjects'};
fignum = 54;
nCol = 3;
mcp = 'cluster_ft';
tit = 'ET-patients vs. CTRL-subj on all memory trials (cue locked)';
plot_ERPcomparisions(avg_memoryCue, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)


%% All memory trials at ALC condition (102:107)
if ~exist(fullfile(ROOTDIR, 'data', 'avg_mem-cuelockedALC_erp.mat'), 'file')
    estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [102:107], ...
        'mem-cuelockedALC', [10,20,110,120])
end
load(fullfile(ROOTDIR, 'data', 'avg_mem-cuelockedALC_erp.mat'))

%% ET-patients vs control subjects memory trials cue locked early ALC
avg_memoryCueALC=cell(1,2);
for k   = 1:2; avg_memoryCueALC{k} = avg(idx_group{k}); end              % separate average in groups
ch      = {'Fz', 'Cz', 'Pz'};                                               % channels to plot
toi     = [-.15 1];                                                         % time of interest (x-axis)
bsl     = [-.15 0];                                                         % baseline to subtract

lgnd = {'ET-patients', 'CTRL-subjects'};
fignum = 55;
nCol = 3;
mcp = 'cluster_ft';
tit = 'ET-patients vs. CTRL-subj on memory trials (cue locked) with ALCOHOL';
plot_ERPcomparisions(avg_memoryCueALC, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)


%% TODO: avg_memoryCue must be turned to one cell (as avg before!!)
plot_erp_amplitudes_boxplot(cat(2, avg_memoryCue{:}), cat(2, avg_memoryCueALC{:}), all_subj, ...
    {[1:19], [20:39]}, 'Fz', [.6 .8], 'p600_memory')


%% Early memory trials at /wo condition (3:5), cue-locked
if ~exist(fullfile(ROOTDIR, 'data', 'avg_Earlymemory-cuelocked_erp.mat'), 'file')
    estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [3:5], ...
        'avg_Earlymemory-cuelocked', [10,20,110,120])
end
load(fullfile(ROOTDIR, 'data', 'avg_Earlymemory-cuelocked_erp.mat'))

%% ET-patients vs control subjects memory trials cue locked
avg_EarlyMemoryCue = avg; avg_shift=cell(1,2);
for k   = 1:2; avg_shift{k} = avg_EarlyMemoryCue(idx_group{k}); end         % separate average in groups
ch      = {'Fz', 'Cz', 'Pz'};                                               % channels to plot
toi     = [-.15 1];                                                         % time of interest (x-axis)
bsl     = [-.15 0];                                                         % baseline to subtract

lgnd = {'ET-patients', 'CTRL-subjects'};
fignum = 74;
nCol = 3;
mcp = 'cluster_ft';
tit = 'ET-patients vs. CTRL-subj on early memory trials (cue locked)';
plot_ERPcomparisions(avg_shift, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)


%% Early memory trials at /wo condition (3:5), cue-locked
if ~exist(fullfile(ROOTDIR, 'data', 'avg_Earlymemory-cuelocked_erp.mat'), 'file')
    estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [103:105], ...
        'Earlymemory-cuelocked', [10,20,110,120])
end
load(fullfile(ROOTDIR, 'data', 'avg_Earlymemory-cuelocked_erp.mat'))

%% ET-patients vs control subjects memory trials cue locked
avg_EarlyMemoryCueALC = avg; avg_earlymemoryALC=cell(1,2);
for k   = 1:2; avg_earlymemory{k} = avg_EarlyMemoryCueALC(idx_group{k}); end   % separate average in groups
ch      = {'Fz', 'Cz', 'Pz'};                                               % channels to plot
toi     = [-.15 1];                                                         % time of interest (x-axis)
bsl     = [-.15 0];                                                         % baseline to subtract

lgnd = {'ET-patients', 'CTRL-subjects'};
fignum = 75;
nCol = 3;
mcp = 'cluster_ft';
tit = 'ET-patients vs. CTRL-subj on early memory trials (cue locked)';
plot_ERPcomparisions(avg_earlymemoryALC, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)


%% All shift trials at /wo condition (10,20)
if ~exist(fullfile(ROOTDIR, 'data', 'avg_MemoryCueLateALC_erp.mat'), 'file')
    estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [106:107], 'MemoryCueLateALC', [10,20,110,120])
end
load(fullfile(ROOTDIR, 'data', 'avg_MemoryCueLateALC_erp.mat'))

%% ET-patients vs control subjects memory trials cue locked
avg_LateMemoryCue  = avg; avg_shift=cell(1,2);
for k   = 1:2; avg_shift{k} = avg_LateMemoryCue (idx_group{k}); end         % separate average in groups
ch      = {'Fz', 'Cz', 'Pz'};                                               % channels to plot
toi     = [-.15 1];                                                         % time of interest (x-axis)
bsl     = [-.15 0];                                                         % baseline to subtract

lgnd = {'ET-patients', 'CTRL-subjects'};
fignum = 156;
nCol = 3;
mcp = 'cluster_ft';
tit = 'ET-patients vs. CTRL-subj on late memory trials (cue locked)';
plot_ERPcomparisions(avg_shift, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)



rest_of_script = 0; % From here on, many different other tests and comparisons were tried, that are not relevant to this work
switch rest_of_script
    case(1)
        %% All shift trials at alcohol condition (110,120)
        if ~exist(fullfile(ROOTDIR, 'data', 'avg_shiftALC_erp.mat'), 'file')
            estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [110, 120], 'shiftALC')
        end
        load(fullfile(ROOTDIR, 'data', 'avg_shiftALC_erp.mat'))
        
        %% ET-patients vs control subjects shift trials
        ch      = {'Fz', 'Cz', 'Pz'};                                               % channels to plot
        toi     = [-.15 1];                                                         % time of interest (x-axis)
        bsl     = [-.15 0];                                                         % baseline to subtract
        
        lgnd    = {'ET-patients', 'CTRL-subjects'};
        fignum  = 114;
        nCol    = 3;
        mcp     = 'cluster_ft';
        tit     = 'ET-patients vs. CTRL-subj on shift trials (ALC)';
        plot_ERPcomparisions(avg_shift2, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)
        
        
        %% Topolots at 600-800msec.
        avg_topoCTRL = cell(1,2);
        avg_topoCTRL{1} = avg_shift{2}; avg_topoCTRL{2} = avg_shiftALC{2};
        lgnd = {'wo/ alcohol', 'alcohol', 'Group difference'};
        fignum = 214;
        tit = 'Topographical distribution of 600-800ms potential on shift-trials (CTRL)';
        plot_topoplots(avg_topoCTRL, fignum, lgnd, tit, [.6 .8])
        
        avg_topoET = cell(1,2);
        avg_topoET{1} = avg_shift{1}; avg_topoET{2} = avg_shiftALC{1};
        lgnd = {'wo/ alcohol', 'alcohol', 'Group difference'};
        fignum = 224;
        tit = 'Topographical distribution of 600-800ms potential on shift-trials (ET)';
        plot_topoplots(avg_topoET, fignum, lgnd, tit, [.6 .8])
        
        
        %% Topoplots of group differences with stats
        avg_shift_diff = cell (1,2); % empty cell for differences in both groups
        
        cfg_math = []; % subtract wo/ condition from ALC condition to use in script
        cfg_math.operation = 'subtract';
        cfg_math.parameter = 'avg';
        for n = 1:2
            avg_shift_diff{n} = arrayfun(@(q) ft_math(cfg_math, avg_shiftALC{n}{q}, avg_shift{n}{q}), ...
                1:numel(avg_shift{n}), 'Un', 0);
        end
        
        %% All memory trials at wo/ alcohol condition (21:25)
        if ~exist(fullfile(ROOTDIR, 'data', 'avg_memory_erp.mat'), 'file')
            estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [21:25], 'memory')
        end
        load(fullfile(ROOTDIR, 'data', 'avg_memory_erp.mat'))
        
        %% ET-patients vs control subjects (all) memory trials
        avg_all = avg; avg_memory=cell(1,2);
        for k = 1:2; avg_memory{k} = avg_all(idx_group{k}); end                            % separate average in groups
        ch              = {'Fz', 'Cz', 'Pz'};                                       % channels of interest
        toi = [-.15 1]; bsl = [-.15 0];
        
        lgnd = {'ET-patients', 'CTRL-subjects'};
        fignum = 11;
        nCol = 3;
        mcp = 'cluster_ft';
        tit = 'ET-patients vs. CTRL-subj on memory trials';
        plot_ERPcomparisions(avg_memory, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)
        
        
        %% All memory trials at alchol condition (121:125)
        if ~exist(fullfile(ROOTDIR, 'data', 'avg_memoryALC_erp.mat'), 'file')
            estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [121:125], 'memoryALC')
        end
        load(fullfile(ROOTDIR, 'data', 'avg_memoryALC_erp.mat'))
        
        %% ET-patients vs control subjects memory trials with alcohol
        avg_all = avg; avg_memoryALC=cell(1,2);
        for k = 1:2; avg_memoryALC{k} = avg_all(idx_group{k}); end                  % separate average in groups
        ch              = {'Fz', 'Cz', 'Pz'};                                       % channels of interest
        toi = [-.15 1]; bsl = [-.15 0];
        
        lgnd = {'ET-patients', 'CTRL-subjects'};
        fignum = 111;
        nCol = 3;
        mcp = 'cluster_ft';
        tit = 'ET-patients vs. CTRL-subj on memory trials (ALC)';
        plot_ERPcomparisions(avg_memoryALC, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)
        
        %% Topolots for the p300 answer for both groups and all memory trials
        lgnd = {'ET-patients', 'CTRL-subjects', 'Group difference'};
        fignum = 20;
        tit = 'Topographical distribution of p300 on memory-trials';
        plot_topoplots(avg_memory, fignum, lgnd, tit, [.2 .4])
        
        %% Differences at specific time points between groups (boxplots)
        load(fullfile(ROOTDIR, 'data', 'avg_memory_erp.mat'))
        avg_memory = avg;
        
        load(fullfile(ROOTDIR, 'data', 'avg_memoryALC_erp.mat'))
        avg_memoryALC = avg;
        
        %%
        plot_erp_amplitudes_boxplot(avg_memory, avg_memoryALC, all_subj, ...
            idx_group, 'POz', [.6 .8])
        
        %% Early memory trials at wo/ alcohol condition (21:23)
        if ~exist(fullfile(ROOTDIR, 'data', 'avg_earlymemory_erp.mat'), 'file')
            estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [21:23], 'earlymemory')
        end
        load(fullfile(ROOTDIR, 'data', 'avg_earlymemory_erp.mat'))
        
        %% ET-patients vs control subjects (early) memory trials (21:23)
        avg_all = avg; avg_earlymemory=cell(1,2);
        for k = 1:2; avg_earlymemory{k} = avg_all(idx_group{k}); end                            % separate average in groups
        ch              = {'Fz', 'Cz', 'Pz'};                         % channels of interest
        toi = [-.15 1]; bsl = [-.15 0];
        
        lgnd = {'ET-patients', 'CTRL-subjects'};
        fignum = 12;
        nCol = 3;
        mcp = 'cluster_ft';
        tit = 'ET-patients vs. CTRL-subj on early memory trials';
        plot_ERPcomparisions(avg_earlymemory, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)
        
        %% Early memory trials at alcohol condition (121:123)
        if ~exist(fullfile(ROOTDIR, 'data', 'avg_earlymemoryALC_erp.mat'), 'file')
            estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [121:123], 'earlymemoryALC')
        end
        load(fullfile(ROOTDIR, 'data', 'avg_earlymemoryALC_erp.mat'))
        
        %% ET-patients vs control subjects (early) memory trials with alcohol (21:23)
        avg_all = avg; avg_earlymemoryALC=cell(1,2);
        for k = 1:2; avg_earlymemoryALC{k} = avg_all(idx_group{k}); end                            % separate average in groups
        ch              = {'Fz', 'Cz', 'POO1'};                         % channels of interest
        toi = [-.15 1]; bsl = [-.15 0];
        
        lgnd = {'ET-patients', 'CTRL-subjects'};
        fignum = 112;
        nCol = 3;
        mcp = 'cluster_ft';
        tit = 'ET-patients vs. CTRL-subj on early memory trials (ALC)';
        plot_ERPcomparisions(avg_earlymemoryALC, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)
        
        
        %% Early memory trials (target_locked, 3:5)
        if ~exist(fullfile(ROOTDIR, 'data', 'avg_targetEarly_erp.mat'), 'file')
            estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [3:5], 'targetEarly')
        end
        load(fullfile(ROOTDIR, 'data', 'avg_targetEarly_erp.mat'))
        
        %% ET-patients vs control subjects (early) memory trials (21:23)
        avg_all = avg; avg_targetEarly=cell(1,2);
        for k = 1:2; avg_targetEarly{k} = avg_all(idx_group{k}); end                            % separate average in groups
        ch              = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz'};                         % channels of interest
        toi = [-.15 1]; bsl = [-.15 0];
        
        lgnd = {'ET-patients', 'CTRL-subjects'};
        fignum = 112;
        nCol = 3;
        mcp = 'cluster_ft';
        tit = 'ET-patients vs. CTRL-subj on early memory trials (target_locked)';
        plot_ERPcomparisions(avg_targetEarly, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)
        
        
        %% Late memory trials (24:25)
        if ~exist(fullfile(ROOTDIR, 'data', 'avg_latememory_erp.mat'), 'file')
            estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [24:25], 'latememory')
        end
        load(fullfile(ROOTDIR, 'data', 'avg_latememory_erp.mat'))
        
        %% ET-patients vs control subjects (late) memory trials (24:25)
        avg_all = avg; avg_latememory=cell(1,2);
        for k = 1:2; avg_latememory{k} = avg_all(idx_group{k}); end                            % separate average in groups
        ch              = {'Fz', 'Cz', 'Pz'};                         % channels of interest
        toi = [-.15 1]; bsl = [-.15 0];
        
        lgnd = {'ET-patients', 'CTRL-subjects'};
        fignum = 13;
        nCol = 3;
        mcp = 'cluster_ft';
        tit = 'ET-patients vs. CTRL-subj on late memory trials';
        plot_ERPcomparisions(avg_latememory, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)
        
        
        %% Differences at specific time points between groups (boxplots)
        load(fullfile(ROOTDIR, 'data', 'avg_shift_erp.mat'))
        avg_shift = avg;
        
        load(fullfile(ROOTDIR, 'data', 'avg_shiftALC_erp.mat'))
        avg_shiftALC = avg;
        
        %%
        plot_erp_amplitudes_boxplot(avg_shift, avg_shiftALC, all_subj, ...
            idx_group, 'POz', [.6 .8])
        
        
        %% Shift scatter plots
        % for g = 1:2 % loop through groups
        %     p = progressbar( numel(all_subj), 'percent' );                           % JSB routine for progress bars
        %     wrong_runs{g} = []; % % pre-allocates space for the error estimation
        %     for k = 1:numel(idx_group{g}) % subj. per group for individual results
        %         p.update( k )
        %         [tmp, ~] = select_runs(sprintf("S%s", ...
        %             num2str(all_subj(idx_group{g}(k)))), ROOTDIR);                   % obtains the error counts from the events-files
        %         wrong_runs{g} = [wrong_runs{g}; tmp];
        %     end
        % end
        % p.stop();
        %
        % avg2test = avg_shift;
        % chtemp = find(strcmp(avg2test{1}{1}.label, 'Pz'));
        % toi = dsearchn(avg2test{1}{1}.time.', [.45 .65].');
        %
        % lcp = nan(max(cell2mat(arrayfun(@(q) numel(avg2test{q}), 1:2, 'Un', 0))), 2);
        % for g = 1:2
        %     lcp(:,g) = cell2mat(arrayfun(@(q) squeeze(mean(avg2test{g}{q}.avg(1,chtemp,toi(1):toi(2)))), 1:numel(avg2test{g}), 'Un', 0)).';
        % end
        %
        % % Scatter plot
        % figure(997);
        % scatter(lcp(:,1), wrong_runs{1}(:,1))
        % [h, p] = corrcoef(lcp(:,1), wrong_runs{1}(:,1))
        %
        % % Scatter plot
        % figure(996);
        % scatter(lcp(:,2), wrong_runs{2}(:,1))
        % [h, p] = corrcoef(lcp(:,2), wrong_runs{2}(:,1))
        
        %% Topolots for the p300 answer for both groups and shift trials
        lgnd = {'ET-patients', 'CTRL-subjects', 'Group difference'};
        fignum = 26;
        tit = 'Topographical distribution of p300 on shift-trials';
        plot_topoplots(avg_shift, fignum, lgnd, tit, [.1 .3])
        
        
        
        %%%%%%%%%%%%%%%%%%%%ALCOHOL CONDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %% Shift trials comparison conditions (110,120)
        
        avg1 = load(fullfile(ROOTDIR, 'data', 'avg_shift_erp.mat'));
        avg2 = load(fullfile(ROOTDIR, 'data', 'avg_shiftALC_erp.mat'));
        
        %% ET-patients WO vs ALC condition in shift trials
        avg_all = avg; avg_shiftCONDshift=cell(1,2);
        for k = 1:2
            if k == 1
                avg_shiftCONDshift{k} = avg1.avg(idx_group{1});
            else
                avg_shiftCONDshift{k} = avg2.avg(idx_group{1});
            end
        end                            % separate average in groups
        ch              = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz'};                         % channels of interest
        toi = [-.15 1]; bsl = [-.15 0];
        
        lgnd = {'wo/ alcohol', 'w/ alcohol'};
        fignum = 999;
        nCol = 3;
        mcp = 'cluster_ft';
        tit = 'ET-patients without alcohol vs. alcohol condition on shift trials';
        plot_ERPcomparisions(avg_shiftCONDshift, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)
        
        
        %% CTRL-subjects WO vs ALC condition in shift trials
        avg_all = avg; avg_shiftCONDshift=cell(1,2);
        for k = 1:2
            if k == 1
                avg_shiftCONDshift{k} = avg1.avg(idx_group{2});
            else
                avg_shiftCONDshift{k} = avg2.avg(idx_group{2});
            end
        end                            % separate average in groups
        ch              = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz'};                         % channels of interest
        toi = [-.15 1]; bsl = [-.15 0];
        
        lgnd = {'wo/ alcohol', 'w/ alcohol'};
        fignum = 998;
        nCol = 3;
        mcp = 'cluster_ft';
        tit = 'CTRL-subjects without alcohol vs. alcohol condition on shift trials';
        plot_ERPcomparisions(avg_shiftCONDshift, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)
        
        %% Memory trials comparison conditions (21:25 and 121:125)
        
        avg1 = load(fullfile(ROOTDIR, 'data', 'avg_memory_erp.mat'));
        avg2 = load(fullfile(ROOTDIR, 'data', 'avg_memoryALC_erp.mat'));
        
        %% ET-patients WO vs ALC condition in memory trials
        avg_all = avg; avg_memoryCOND=cell(1,2);
        for k = 1:2
            if k == 1
                avg_memoryCOND{k} = avg1.avg(idx_group{1});
            else
                avg_memoryCOND{k} = avg2.avg(idx_group{1});
            end
        end                            % separate average in groups
        ch  = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz'};                         % channels of interest
        toi = [-.15 1]; bsl = [-.15 0];
        
        lgnd = {'wo/ alcohol', 'w/ alcohol'};
        fignum = 997;
        nCol = 3;
        mcp = 'cluster_ft';
        tit = 'ET-patients without alcohol vs. alcohol condition on memory trials';
        plot_ERPcomparisions(avg_memoryCOND, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)
        
        
        %% CTRL-subjects WO vs ALC condition in shift trials
        avg_all = avg; avg_memoryCOND=cell(1,2);
        for k = 1:2
            if k == 1
                avg_memoryCOND{k} = avg1.avg(idx_group{2});
            else
                avg_memoryCOND{k} = avg2.avg(idx_group{2});
            end
        end                            % separate average in groups
        ch = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz'};                         % channels of interest
        toi = [-.15 1]; bsl = [-.15 0];
        
        lgnd = {'wo/ alcohol', 'w/ alcohol'};
        fignum = 996;
        nCol = 3;
        mcp = 'cluster_ft';
        tit = 'CTRL-subjects without alcohol vs. alcohol condition on memory trials';
        plot_ERPcomparisions(avg_memoryCOND, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)
        
        %% Shift trials comparison conditions (10,20 and 110,120)
        
        avg1 = load(fullfile(ROOTDIR, 'data', 'avg_shift_erp.mat'));
        avg2 = load(fullfile(ROOTDIR, 'data', 'avg_shiftALC_erp.mat'));
        
end




