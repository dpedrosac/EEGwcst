function estimate_averages(all_subj, wdir)

%   This function loads all data and estimates the metrics needed, saving
%   them as structures (avg*x*). As a new feature, there is now a
%   conversion with ft_struct2single(x)

%   Copyright (C) April 2017, modified Mai 2021
%   D. Pedrosa, University Hospital of Gießen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

toi = {             % Different trial number available; for further details cf: ./paradigm/WCST-Tremorparadigma.sce
    [21,22,23], ... % Early trials (wo/ alcohol), formerly [2,3]
    [24,25], ...    % Late trials (wo/ alcohol), formerly [6,7]
    [102:103], ...  % Early trials (w/ alcohol)
    [106:107], ...  % Late trials (w/ alcohol) 
    [21:25], ...    % Right trials (wo/ alcohol)
    [121:125], ...  % Right trials (w/ alcohol)
    [10, 20], ...   % Shift trials(wo/ alcohol)
    [110, 120]};    % Shift trials(w/ alcohol)              %#ok<*NBRAK>

%%
iter = 0; clear avg*;                                          %#ok<*NBRAK>

for dId = 1:numel(all_subj)
    fprintf('\nthe subject being processed is: %s \n', ...
        strcat('S', num2str(all_subj(dId))));
    filename = fullfile(wdir, ['S', num2str(all_subj(dId))], ...
        ['data_final_S', num2str(all_subj(dId)),'.mat']); clear data_final
    load(filename); %#ok<LOAD>
    iter = iter +1;
    for t = 1:numel(toi)
        cfg = [];
        cfg.trials = find(ismember(data_final.trialinfo, toi{t}));
        data_temp = ft_selectdata(cfg, data_final);
        
        cfg = [];
        cfg.keeptrials = 'yes';
        data_temp = ft_timelockanalysis(cfg, data_temp);
        data_temp = ft_struct2single(data_temp); %#ok<NASGU>
        avg{t}{iter} = data_temp;        
        clear data_temp
    end
    clear data_final*;
end

for i = 1:numel(avg)
    eval(['avg',num2str(i),'=avg{',num2str(i),'};']);
end

save(fullfile(wdir,'avgs.mat'),'-regexp','avg[1-9]'); 