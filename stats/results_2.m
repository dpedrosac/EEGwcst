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
excl_subj = [ 16, 21, 45 ];                                                     % excluded subjects due to technical reasons
subj1proc = subj1(~ismember(subj1, excl_subj));
subj2proc = subj2(~ismember(subj2, excl_subj));

idx_group       = {find(ismember(all_subj, subj1proc)), ...                     % get 'group indices'
    find(ismember(all_subj, subj2proc))};

%% All shift trials at /wo condition (10,20)
if ~exist(fullfile(ROOTDIR, 'data', 'avg_shift_erp.mat'), 'file')
    estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [10, 20], 'shift')
end
load(fullfile(ROOTDIR, 'data', 'avg_shift_erp.mat'))

%% ET-patients vs control subjects shift trials
avg_allshift1 = avg; avg_shift=cell(1,2);
for k = 1:2; avg_shift{k} = avg_allshift1(idx_group{k}); end                % separate average in groups
% ch              = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz', 'POz'};                % channels of interest
ch              = {'Fz', 'Cz', 'Pz'};                                       % channels of interest
toi = [-.15 1]; bsl = [-.15 0];

lgnd = {'ET-patients', 'CTRL-subjects'};
fignum = 14;
nCol = 3;
mcp = 'cluster_ft';
tit = 'ET-patients vs. CTRL-subj on shift trials';
plot_ERPcomparisions(avg_shift, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)

%% All shift trials at alcohol condition (110,120)
if ~exist(fullfile(ROOTDIR, 'data', 'avg_shiftALC_erp.mat'), 'file')
    estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [110, 120], 'shiftALC')
end
load(fullfile(ROOTDIR, 'data', 'avg_shiftALC_erp.mat'))

%% ET-patients vs control subjects shift trials
avg_allshift2 = avg; avg_shiftALC=cell(1,2);
for k = 1:2; avg_shiftALC{k} = avg_allshift2(idx_group{k}); end             % separate average in groups
%ch              = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz', 'POz'};                 % channels of interest
ch              = {'Fz', 'Cz', 'Pz'};                                       % channels of interest
toi = [-.15 1]; bsl = [-.15 0];

lgnd = {'ET-patients', 'CTRL-subjects'};
fignum = 114;
nCol = 3;
mcp = 'cluster_ft';
tit = 'ET-patients vs. CTRL-subj on shift trials (ALC)';
plot_ERPcomparisions(avg_shiftALC, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)

%% Estimate differences at specific time points and between conditions,
%   that is repeated measures ANOVA

plot_erp_amplitudes_boxplot(avg_allshift1, avg_allshift2, all_subj, ...
idx_group, 'Pz', [.6 .8], 'p600-800_shift')

%% Insert boxplots here

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
ch              = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz'};                         % channels of interest
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
ch              = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz'};                         % channels of interest
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
ch              = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz'};                         % channels of interest
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






