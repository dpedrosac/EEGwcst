function stats_and_results

%   This function runs all analyses for demographics, reaction times and
%   different errors between the two groups (CTRL-subjects and ET-patients)

%   Copyright (C) January 2018, modified June 2021
%   D. Pedrosa, Urs Kleinholdermann University Hospital of Gießen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

%% General settings and indices to use later in the analyses
[~, ROOTDIR] = EEGwcst_defaults; 
subj1   = [3,4,5,7,8,10,13,14,16,18,19,20,21, 22,23,24,42,44,45];           % ET-Patients
subj2   = [1,2,27,28,29,30,32,35,36,37,38,39,41,43,46,47,48,49,50];         % CTRL-subjects

load(fullfile(ROOTDIR, 'data', 'patdat.mat'));                   %#ok<LOAD> % load metadata

subj1cll = []; subj2cll = [];                                               % clear variables from workspace
for m = 1:numel(subj1); subj1cll{m} = strcat('s', num2str(subj1(m))); end
idx_et = find(ismember({patient.code}, subj1cll));
for m = 1:numel(subj2); subj2cll{m} = strcat('s', num2str(subj2(m))); end   %#ok<*AGROW> % create a cell array to compare the codes in control/patient with the data of the subjects
idx_ctrl = find(ismember({control.code}, subj2cll));


%% First results come from the general data, the reaction times and the 
% number of total errors

results_1(control(idx_ctrl),patient(idx_et),subj1, subj2, ROOTDIR)
