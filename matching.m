function matching(control, patient, save_dir)

%   This function creates a table, in which data is arranged in a way, that
%   a matching algorithm may be run in R. Specifically, both groups are
%   matched for gender and for age and eductation using the Mahalanobis
%   distance (for details, see MatchIt documentation)

%   Copyright (C) January 2018
%   D. Pedrosa, University Hospital of Gieﬂen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.


dattemp = [];                                                               % pre-allocate the space
idx = {[1:3, 5:9, 11:13, 15:17, 21:23], ...
    [1:7, 9:14, 16:22]};                                                    % indices of patients who need a 'matching partner'
Xs = []; dt = []; sub = [0, -2]; code_all = {};

for g = 1:2 % loop through both grous (ET-patients and ctrl-subj.)
    if g == 1; dattemp{g} = patient(idx{g}); else dattemp{g} = control(idx{g}); end % in order to make the code more efficient, a temporal data structure is defined
    for k = 1:numel(dattemp{g}) % loop through all participants in a group
        dt(k,1) = g+sub(g);                                                 % treatment, that is either (1) = ET or (0) = ctrl
        dt(k,2) = dattemp{g}(k).age;                                        % age
        dt(k,3) = dattemp{g}(k).gender;                                     % gender, that is either (1) = women, (0) = men
        dt(k,4) = dattemp{g}(k).demtect;                                    % total demtect score (not needed)
        try dt(k,5) = dattemp{g}(k).edu; catch; dt(k,5) = nan; end          % education years
        code{k,1} = dattemp{g}(k).code;                                     % patient code, to make teh identification easier
    end
    Xs = [Xs; dt]; code_all = cat(1,code_all, code); clear code             % (Xs) is the matrix with all values together
end
%%

tble = table( code_all, Xs(:,1), Xs(:,2), Xs(:,3), Xs(:,4), Xs(:,5), ...    % create a table that may be read in R
    'VariableNames', {'code', 'treatment', 'age', 'gender', 'demtect', 'edu'} );
writetable(tble, 'match_it.txt', 'Delimiter', '\t');

%%
match_matrix = [26,22,23,21,32,31,25,30,35,28,29,37,24,18,33,36,19];        % this matrix results from the algorithm in R, which aims at pairing the ET-patients and the control subjects

match_matrix = [21,22,33,26,34, 31, 25, 36, 24, 28, 29, 37, 23, 19, 35, 30, 18];

for k = 1:numel(control); control(k).matching = []; end
for k = 1:numel(patient); patient(k).matching = []; end

for m = 1:length(match_matrix)
    idx = find(strcmp({patient(:).code}, code_all{m}));                     % provides an index for the patient being analysed
    tmp = strsplit(code_all{match_matrix(m)}, 's');                         % remove the 'S' part of the code
    patient(idx).matching = code_all{match_matrix(m)};                      % sets the matching partner into the metadata file (patient)
    
    idx2 = find(strcmp({control(:).code}, code_all{match_matrix(m)}));      % provides the index for the control-subject corresponding to the first ET-patient
    control(idx2).matching = code_all{m};                                   % sets the matching partner into the metadata file (control)
end

save(strcat(save_dir, '\', 'patdat.mat'), 'control', 'patient', '-v7.3')