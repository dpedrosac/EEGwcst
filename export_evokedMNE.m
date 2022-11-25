
%%
conds = {'shift', 'shiftALC'}; % 'mem-cuelocked', 'mem-cuelockedALC', '
for c = 1:numel(conds)
    preprocess_cue(conds{c})
end

function preprocess_cue(cond)

load(sprintf('/media/storage/skripte/lambda/data/avg_%s_erp.mat', cond));
avg_temp = avg;
baseline_correct = 1;
reverseStr = '';
switch baseline_correct
    case (1) % perform baseline correction
        cfg = [];
        cfg.baseline = [-0.15 0];
        cfg.baselinetype = 'absolute';
        cfg.feedback = 'no';
        fprintf('\n')
        for k = 1:numel(avg)
            % Display the progress
            percentDone = 100 * k / numel(avg);
            msg = sprintf('Percent done: %3.1f', percentDone); %Don't forget this semicolon
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            avg_temp{k} = ft_timelockbaseline(cfg, avg_temp{k});
        end
end

cfg = [];
cfg.keepindividual = 'yes';
dataET = ft_timelockgrandaverage(cfg, avg_temp{idx_group{1}});
dataCTRL = ft_timelockgrandaverage(cfg, avg_temp{idx_group{2}});

save(sprintf('/media/storage/skripte/lambda/data/dataMNE-ET_%s_erp.mat', cond), ...
    'dataET', '-v7.3');
save(sprintf('/media/storage/skripte/lambda/data/dataMNE-CTRL_%s_erp.mat', cond), ...
    'dataCTRL', '-v7.3');

% fiff_file  = 'ctf_raw.fif';
% fieldtrip2fiff(fiff_file, avg_temp{1})

%% Define neighbours
cfg = [];
cfg.method = 'distance';
cfg.feedback = 'no';                                                       % can be set to 'yes' to see results
neighbours = ft_prepare_neighbours(cfg, avg_temp{1});

%% Start statistics
cfg                 = [];
cfg.channel          = 'all';
cfg.latency          = [0 1.2];
cfg.avgovergchan     = 'no';
cfg.method           = 'ft_statistics_montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'nonparametric_individual';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = cfg.tail;
cfg.alpha            = 0.05;
cfg.correcttail      = 'alpha';
cfg.computeprob      = 'yes';
cfg.numrandomization = 500;
cfg.neighbours       = neighbours;

design = zeros(1,size(dataET.individual,1) + size(dataCTRL.individual,1));
design(1,1:size(dataET.individual,1)) = 1;
design(1,(size(dataET.individual,1)+1):(size(dataET.individual,1)+size(dataCTRL.individual,1))) = 2;

cfg.design = design;
cfg.ivar   = 1;

stat3 = ft_timelockstatistics(cfg, dataET, dataCTRL);

%% Estimate average
cfg = [];
cfg.keepindividual = 'no';
avgET = ft_timelockgrandaverage(cfg, avg_temp{idx_group{1}});
avgCTRL = ft_timelockgrandaverage(cfg, avg_temp{idx_group{2}});

% Then take the difference of the averages using ft_math
cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffectETvsCTRL = ft_math(cfg, avgET, avgCTRL);

%% Start plotting
timestep      = 0.05; % timestep between time windows for each subplot (in seconds)
sampling_rate = 200; % Data has a temporal resolution of 300 Hz
sample_count  = length(stat.time);
% number of temporal samples in the statistics object
j = [0:timestep:1]; % Temporal endpoints (in seconds) of the ERP average computed ineach subplot
m = [1:timestep*sampling_rate:sample_count]; % temporal endpoints in M/EEG samples

for k = 1:20
    subplot(4,5,k);
    cfg = [];
    cfg.xlim = [j(k) j(k+1)];   % time interval of the subplot
    %cfg.zlim = [-2.5e-13 2.5e-13];
    % If a channel is in a to-be-plotted cluster, then
    % the element of pos_int with an index equal to that channel
    % number will be set to 1 (otherwise 0).
    
    % Next, check which channels are in the clusters over the
    % entire time interval of interest.
    %    pos_int = zeros(numel(raweffectFICvsFC.label),1);
    %    neg_int = zeros(numel(raweffectFICvsFC.label),1);
    %    pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
    %    neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);
    
    cfg.highlight   = 'on';
    % Get the index of the to-be-highlighted channel
    %    cfg.highlightchannel = find(pos_int | neg_int);
    cfg.comment     = 'xlim';
    cfg.commentpos  = 'title';
    cfg.layout      = 'EEG1005.lay';
    cfg.interactive = 'no';
    cfg.figure      = 'gca'; % plots in the current axes, here in a subplot
    ft_topoplotER(cfg, raweffectETvsCTRL);
end

%%
figure
cfg = [];
ft_multiplotER(cfg, dataET, dataCTRL);
end