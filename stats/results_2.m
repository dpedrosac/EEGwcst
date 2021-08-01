function results_2(control, patient, subj1, subj2, ROOTDIR, wdir)

%   Generate ERP between both groups with several available options

%   Copyright (C) July 2021
%   D. Pedrosa University Hospital of Gießen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

%% All memory trials (21:25)
cd(ROOTDIR)
print_legend
if ~exist(fullfile(ROOTDIR, 'data', 'avg_memory_erp.mat'), 'file')
   estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [21:25], 'memory') 
end
load(fullfile(ROOTDIR, 'data', 'avg_memory_erp.mat'))

%% ET-patients vs control subjects (all) memory trials
all_subj = [subj1, subj2];
avg_all = avg; avg_memory=cell(1,2);
idx_group       = {find(ismember(all_subj, subj1)), ...                     % get 'group indices'
    find(ismember(all_subj, subj2))};
for k = 1:2; avg_memory{k} = avg_all(idx_group{k}); end                            % separate average in groups
ch              = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz'};                         % channels of interest
toi = [-.15 1]; bsl = [-.15 0];

lgnd = {'ET-patients', 'CTRL-subjects'};
fignum = 11;
nCol = 3;
mcp = 'cluster_ft';
tit = 'ET-patients vs. CTRL-subj on memory trials';
plot_ERPcomparisions(avg_memory, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)

%% Early memory trials (21:23)
if ~exist(fullfile(ROOTDIR, 'data', 'avg_earlymemory_erp.mat'), 'file')
   estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [21:23], 'earlymemory') 
end
load(fullfile(ROOTDIR, 'data', 'avg_earlymemory_erp.mat'))

%% ET-patients vs control subjects (early) memory trials (21:23)
avg_all = avg; avg_earlymemory=cell(1,2);
idx_group       = {find(ismember(all_subj, subj1)), ...                     % get 'group indices'
    find(ismember(all_subj, subj2))};
for k = 1:2; avg_earlymemory{k} = avg_all(idx_group{k}); end                            % separate average in groups
ch              = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz'};                         % channels of interest
toi = [-.15 1]; bsl = [-.15 0];

lgnd = {'ET-patients', 'CTRL-subjects'};
fignum = 12;
nCol = 3;
mcp = 'cluster_ft';
tit = 'ET-patients vs. CTRL-subj on early memory trials';
plot_ERPcomparisions(avg_earlymemory, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)

%% Late memory trials (24:25)
if ~exist(fullfile(ROOTDIR, 'data', 'avg_latememory_erp.mat'), 'file')
   estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [24:25], 'latememory') 
end
load(fullfile(ROOTDIR, 'data', 'avg_latememory_erp.mat'))

%% ET-patients vs control subjects (late) memory trials (24:25)
avg_all = avg; avg_latememory=cell(1,2);
idx_group       = {find(ismember(all_subj, subj1)), ...                     % get 'group indices'
    find(ismember(all_subj, subj2))};
for k = 1:2; avg_latememory{k} = avg_all(idx_group{k}); end                            % separate average in groups
ch              = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz'};                         % channels of interest
toi = [-.15 1]; bsl = [-.15 0];

lgnd = {'ET-patients', 'CTRL-subjects'};
fignum = 13;
nCol = 3;
mcp = 'cluster_ft';
tit = 'ET-patients vs. CTRL-subj on late memory trials';
plot_ERPcomparisions(avg_latememory, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)


%% Topolots for the p300 answer for both groups and all memory trials
lgnd = {'ET-patients', 'CTRL-subjects', 'Group difference'};
fignum = 20;
tit = 'Topographical distribution of p300 on memory-trials';
plot_topoplots(avg_memory, fignum, lgnd, tit, [.2 .4])

%% All memory trials (21:25)
if ~exist(fullfile(ROOTDIR, 'data', 'avg_shift_erp.mat'), 'file')
   estimate_erp_averages_generic([subj1, subj2], ROOTDIR, wdir, [10, 20], 'shift') 
end
load(fullfile(ROOTDIR, 'data', 'avg_shift_erp.mat'))

%% ET-patients vs control subjects shift trials
all_subj = [subj1, subj2];
avg_all = avg; avg_shift=cell(1,2);
idx_group       = {find(ismember(all_subj, subj1)), ...                     % get 'group indices'
    find(ismember(all_subj, subj2))};
for k = 1:2; avg_shift{k} = avg_all(idx_group{k}); end                            % separate average in groups
ch              = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz', 'POz'};                         % channels of interest
toi = [-.15 1]; bsl = [-.15 0];

lgnd = {'ET-patients', 'CTRL-subjects'};
fignum = 14;
nCol = 3;
mcp = 'cluster_ft';
tit = 'ET-patients vs. CTRL-subj on shift trials';
plot_ERPcomparisions(avg_shift, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)

%% Topolots for the p300 answer for both groups and shift trials
lgnd = {'ET-patients', 'CTRL-subjects', 'Group difference'};
fignum = 26;
tit = 'Topographical distribution of p300 on shift-trials';
plot_topoplots(avg_shift, fignum, lgnd, tit, [.1 .3])

