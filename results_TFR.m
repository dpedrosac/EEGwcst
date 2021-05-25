wdir = 'C:\Users\David\projekte\';                                           % (wdir) defines the working drive for the script
loaddir = [wdir 'wcst\'];

dist = 30;                                      % distance neighbouring information is taken from (see define_neighbours)
neighbours = define_neighbours(data_merged, loaddir, dist, 0);                 % runs the function in which definition of neighbours is
ch_2plot = 'Pz';
idx_ch = neighbours(find(strcmp({neighbours(:).label}, ch_2plot))).neighblabel; %#ok<FNDSB>
trlsnum = [10, 20];

cfg = [];
cfg.trials = find(ismember(data_merged.trialinfo, trlsnum));
cfg.channel = idx_ch;
data_temp = ft_selectdata(cfg,data_final);

cfg = [];
cfg.output = 'pow';
cfg.channel = 'all';
cfg.method = 'mtmconvol';
cfg.foi = 8:.5:90;

cfg.t_ftimwin = 7./cfg.foi;
cfg.tapsmofrq = .3*cfg.foi;
cfg.toi = 'all';%-.5:.05:1.5;
TFRmult = ft_freqanalysis(cfg, data_temp);


%%
figure;
cfg = [];
cfg.baseline = [-.4 0];
cfg.baselinetype = 'db';

cfg.xlim = [-.4 1.0];
ft_singleplotTFR(cfg, TFR_avg);