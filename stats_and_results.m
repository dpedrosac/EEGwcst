unction stats_and_results

%% General settings and indices to use later in the analyses
wdir    = 'C:\Users\David\projekte\';
subj1   = [1,2,27,28,29,30,35,37,38,39,41,46,47,48,49,50];                    % CTRL-subjects
subj2   = [3,4,5,7,8,10,14,18,19,20,22,23,24,42,44,45];                       % ET-Patients
subj    = [subj1, subj2];
load(strcat(wdir, 'wcst\patdat.mat'));                                      %load metadata

subj1cll = []; subj2cll = [];                                               % clear variables from workspace
for m = 1:numel(subj1); subj1cll{m} = strcat('s', num2str(subj1(m))); end   %#ok<*AGROW> % create a cell array to compare the codes in control/patient with the data of the subjects
idx_ctrl = find(ismember({control.code}, subj1cll));
for m = 1:numel(subj2); subj2cll{m} = strcat('s', num2str(subj2(m))); end
idx_et = find(ismember({patient.code}, subj2cll));

%% First results come from the general data and two ANOVAS for the effects
%of alcohol on the amounr of two distinct errors and the effect of alcohol
% on two different repsonse times

npat = numel(idx_et);                                                       % number of patients of all included subjects
general_results(control(idx_ctrl),patient(idx_et),npat,subj1, subj2, wdir)  
%   displays general results for the data and the amount of errors in the 
%   groups and prepares data for ANOVAs in R and plots interaction effects

%% sensor level analysis
wdir                = 'C:\Users\David\projekte\';
conds.name          = {'rightWO', 'wrongWO2'};
conds.sb_title      = {'positive feedback', 'negative feedback'};
conds.ch            = {'Fz', 'Cz', 'Pz'};
conds.mcp.m         = 'cluster_ft';                                         % how to adjust for Type-I errors; two possible options: 'none' and 'mcp. the latter has many additional options to set
conds.mcp.method    = 'montecarlo';
conds.mcp.correctm  = 'cluster';
conds.mcp.numrand   = 5000;
conds.mcp.statistic = 'indepsamplesT';
conds.mcp.clusterthreshold = 'nonparametric_individual';
conds.mcp.alpha = .05;
conds.mcp.clusteralpha = .05;
conds.mcp.correcttail = 'alpha';
conds.mcp.clusterstatistic = 'maxsum';
conds.cmp.wcm_weight = .5;
conds.mcp.dimord = 'chan_time';
conds.mcp.connectivity = 0;

plotERP(subj, conds, wdir, subj1, subj2, 21)

%% sensor level analysis

wdir            = 'C:\Users\David\projekte\';
conds.name      = {'wrongWO', 'wrongALC'};
conds.sb_title  = {'wo\ alcoholtrials', 'w\ alcohol'};
conds.ch        = {'Fz', 'Cz', 'Pz'};
conds.mcp       = 'cluster_ft';                                                 % how to adjust for Type-I error correction

plotERP(subj, conds, wdir, subj1, subj2, 92)


