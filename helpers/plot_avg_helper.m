

for k = idx_group{1}
    close all;
    figure; hold on;
    fprintf('\nSubject being processed is no: %d (S%d)\n', k,  all_subj(k))
    data_temp = avg{k};
    mean_data = squeeze(nanmean(avg{k}.trial(:,112,:),1));
    if isstruct(data_temp)
        plot(avg{k}.time, squeeze(avg{k}.trial(:,112,:)).');
        plot(avg{k}.time, mean_data, '--k', 'LineWidth', 5)
        xlim([-.2 1])
        keyboard;
    else
        continue
    end
end
close all;

subj1   = [3,4,5,6,7,8,10,13,14,16,18,19,20,21,22,23,24,42,44,45];          % ET-Patients
subj2   = [41, 38, 50, 9, 27,47,36,32,43,48,28,35,37,29,49,2,30, 46,39,1]; %48, 49,        % CTRL-subjects

%%
% 1 -> 1 huge artifact in the ALC condition
% 30 -> ALC condition
% 21 -> blink artifacts
% 49 -> not many channels, blinks!

%%% Re estimate bad trtuals for subject
num = 49;
filename = fullfile(wdir, 'data_final', sprintf('data_final_erp_S%d.mat', num));
load(filename)

cfg = [];
cfg.bpfiler = 'yes';
cfg.bpfreq = [.1 30];
data_temp = ft_preprocessing(cfg, data_final);

cfg = [];
cfg.keeptrials = 'yes';

data_temp = ft_timelockanalysis(cfg, data_temp);
data_temp.avg = nanmean(data_temp.trial, 1);
data_temp = ft_struct2single(data_temp);                                % to save space, data is converted to "single"

figure; hold on;
trlnum = find(data_temp.trialinfo >=2 & data_temp.trialinfo <=7);
%trlnum = find(data_temp.trialinfo ==10 | data_temp.trialinfo ==20);
subplot(1,2,1); hold on;
plot(data_temp.time, squeeze(data_temp.trial(trlnum,112,:)).');
plot(data_temp.time, squeeze(nanmean(data_temp.trial(trlnum,112,:),1)), '--k', 'LineWidth', 5)
xlim([-.2 1])           
ylim([-1 1].*25)

subplot(1,2,2); hold on;
trlnum = find(data_temp.trialinfo >=102 & data_temp.trialinfo <=107);
%trlnum = find(data_temp.trialinfo ==110 | data_temp.trialinfo ==120);
plot(data_temp.time, squeeze(data_temp.trial(trlnum,112,:)).');
plot(data_temp.time, squeeze(nanmean(data_temp.trial(trlnum,112,:),1)), '--k', 'LineWidth', 5)
xlim([-.2 1])           
ylim([-1 1].*25)





