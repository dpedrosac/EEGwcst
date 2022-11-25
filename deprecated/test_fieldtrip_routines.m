
cfg = [];
cfg.method = 'distance'; %'template' or 'distance' or 'triangulation'
cfg.neighbourdist = .07; % maximum distance between neighbouring sensors...
                        % (only for 'distance')
cfg.feedback = 'no'; % show a neighbour plot
neighbours = ft_prepare_neighbours(cfg,data_final.elec);


idx_all = [idx_group{1}, idx_group{2}];
grandavg = cell(1,36);
for k = 1:numel(idx_all) 
    
    load(fullfile(wdir, 'data_final', sprintf('data_final_erp_S%d.mat', all_subj(idx_all(k)))));
    
    cfg = [];
    cfg.detrend = 'yes';
    cfg.demean = 'yes';
    cfg.baseline = [-.2 0];
    dattemp = ft_preprocessing(cfg, data_final);
    
    cfg = [];
    cfg.latency = [-.2 1];
    cfg.trials = find(data_final.trialinfo>=21 & data_final.trialinfo<=25);    
    dattemp = ft_timelockanalysis(cfg, dattemp);
    
    cfg = [];
    grandavg{k} = ft_timelockgrandaverage(cfg, dattemp);
end 
      
cfg = [];
grandavg1 = ft_timelockgrandaverage(cfg, grandavg{1:16});

cfg = [];
cfg.baseline = [-.1 0];
grandavg1 = ft_timelockbaseline(cfg, grandavg1);

cfg = [];
grandavg2 = ft_timelockgrandaverage(cfg, grandavg{17:36});

cfg = [];
cfg.baseline = [-.1 0];
grandavg2 = ft_timelockbaseline(cfg, grandavg2);


cfg                  = [];
cfg.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
cfg.statistic        = 'indepsamplesT'; % use the independent samples T-statistic as a measure to
                                   % evaluate the effect at the sample level
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;       % alpha level of the sample-specific test statistic that
                                   % will be used for thresholding
cfg.clusterstatistic = 'maxsum';   % test statistic that will be evaluated under the
                                   % permutation distribution.
cfg.minnbchan        = 0;          % minimum number of neighborhood channels that is
                                   % required for a selected sample to be included
                                   % in the clustering algorithm (default=0).
cfg.neighbours     = neighbours; % see below
cfg.tail             = 0;          % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      = 0;
cfg.alpha            = 0.05;      % alpha level of the permutation test
cfg.numrandomization = 500;        % number of draws from the permutation distribution
cfg.neighbours       = neighbours;

cfg.design           = [ones(1,numel(idx_group{1})), ones(1,numel(idx_group{2}))*2]; % design matrix
cfg.ivar             = 1; % number or list with indices indicating the independent variable(s)
[stat] = ft_timelockstatistics(cfg, grandavg{1:16}, grandavg{17:36});


templateLayout = 'EEG1005.lay'; % one of the template layouts included
                                % in FieldTrip
cfg = [];
cfg.layout = which(templateLayout); 
layout = ft_prepare_layout(cfg);

% Plotting stat output
cfg = [];
cfg.layout = layout;
%cfg.parameter = 'stat';
%cfg.maskparameter = 'mask';
cfg.graphcolor = 'r';
cfg.showcomment = 'no';
cfg.showscale = 'no';
figure;
ft_multiplotER(cfg,grandavg1, grandavg2);



