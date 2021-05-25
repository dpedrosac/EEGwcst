subj = subj1;


for k = 1:numel(subj)
   load(strcat(wdir, 'wcst\','S', num2str(subj(k)), '\TFRwrongWO_S', num2str(subj(k)), '.mat'))
   cfg = [];
   TFRall{k} = ft_freqdescriptives(cfg, TFR_avg);
end

%%
cfg = [];
TFRga = ft_freqgrandaverage(cfg, TFRall{:});

%%
figure(345);
cfg = [];
cfg.baseline = [-.1 0];
cfg.baselinetype = 'db';
cfg.zlim = 'maxabs';
cfg.xlim = [-.2 1.2];
cfg.colormap =  othercolor('David1');  
ft_singleplotTFR(cfg, TFRga);