function run_groupstatisticsTFR(subj, wdir, outdir_tfr)

%   This function runs statistcis for the TFR representation of both
%   groups

%   ## Version 1.0 # first version

%   Copyright (C) December 2021
%   D. Pedrosa
%   University Hospital of Gießen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

load(fullfile(wdir, "patdat.mat"))
seq = 'subject'; iter = 0;                                                  % get metadata correctly along with next line
structTFR = cell(1, size(cat(2,subj{:}),2));
for k = 1:2 % load data for both groups
    tolom = subj{k}; % selects whether (1) control subjects or (2) patients are analysed
    if k == 1; temp = control; else; temp = patient; end
    codes = {temp(tolom).code};
    
    for p = 1:numel(codes) % loop through participants
        iter = iter + 1;
        filename_tfr_pre = fullfile(outdir_tfr, sprintf('tfr_%s', ...       % first part of filename which is completed later
            upper(codes{p})));
        load(sprintf('%s_memoryWO.mat', filename_tfr_pre))
        
        cfg_BL = [];
        cfg_BL.baselinetype= 'relchange';%'db';
        cfg_BL.baseline    = [-0.3 0];
        structTFR{iter} = tfreq; %ft_freqbaseline(cfg_BL, tfreq)
    end
end

%% Method = analytic

cfg = [];
cfg.method = 'analytic'; %% do a mass-univariate test
cfg.latency = [0 1.5];
cfg.alpha = 0.05; %% critical value around ± 2.09
cfg.statistic = 'indepsamplesT'; % use a dependent samples t-test
%cfg.correctm = 'holm';
cfg.design(1, :) = [1:numel(subj{1}) 1:numel(subj{2})];
cfg.design(2, :) = [ones(1, numel(subj{1})) 2 * ones(1, numel(subj{2}))];
%cfg.uvar = 1; % first row of cfg.design, containing the (u)nits (subjects)
cfg.ivar = 2; % second row of cfg.design, the (i)ndependent events (3&15)
stat = ft_freqstatistics(cfg, structTFR{:});

cfg = [];
cfg.parameter     = 'stat';
cfg.maskparameter = 'mask';
cfg.maskalpha     = .15; % opacity value for non-significant parts
%cfg.maskstyle     = 'outline';
cfgl.layout       = 'EEG1005.lay';
cfg.layout        = ft_prepare_layout(cfgl); clear cfgl;                    % prepares the layout according to the definition before
% cfg.colormap      = othercolor('David1');                                   % defines the colormap to use, see also omikron/othercolor
cfg.xlim          = [-.2 1.3];
cfg.ylim          = [12 85];
cfg.gridscale     = 40;                                                     % number of interpolation points, the higher the better the resolution
figure, ft_multiplotTFR(cfg,stat2)

%% Same Analysis with cluster randomised statistics

% Cluster based permutation
cfg = [];
cfg.channel          = 'EEG';
cfg.frequency        = [12 50];
cfg.latency          = [0 1.2];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 1;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, structTFR{1});

cfg.design(1, :) = [1:numel(subj{1}) 1:numel(subj{2})];
cfg.design(2, :) = [ones(1, numel(subj{1})) 2 * ones(1, numel(subj{2}))];

% design = zeros(1,size(freqFIC_planar_cmb.powspctrm,1) + size(freqFC_planar_cmb.powspctrm,1));
% design(1,1:size(structTFR{1}.powspctrm,1)) = 1;
% design(1,(size(freqFIC_planar_cmb.powspctrm,1)+1):(size(freqFIC_planar_cmb.powspctrm,1)+...
%     size(freqFC_planar_cmb.powspctrm,1))) = 2;

cfg.ivar             = 2;

[stat2] = ft_freqstatistics(cfg, structTFR{:});

%%
cfg = [];
effectCTRL = ft_freqgrandaverage(cfg, structTFR{1:20});
effectET = ft_freqgrandaverage(cfg, structTFR{21:end});

cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'subtract';%'(x1-x2) / (x1+x2)';
effectDiff = ft_math(cfg, effectCTRL, effectET);


cfg = [];
cfg.frequency = [stat2.freq(1) stat2.freq(end)];
cfg.latency   = [stat2.time(1) stat2.time(end)];

this_tfr = ft_selectdata(cfg, effectDiff);
this_tfr.mask = stat2.mask;

cfg = [];
cfg.parameter     = 'powspctrm';%'stat';
cfg.maskparameter = 'mask';
cfg.maskalpha     = .15; % opacity value for non-significant parts
cfg.maskstyle     = 'outline';
cfgl.layout       = 'EEG1005.lay';
cfg.layout        = ft_prepare_layout(cfgl); clear cfgl;                    % prepares the layout according to the definition before
% cfg.colormap      = othercolor('David1');                                   % defines the colormap to use, see also omikron/othercolor
%cfg.xlim          = [-.3 1.5];
%cfg.ylim          = [15 30];
%cfg.gridscale     = 40;                                                     % number of interpolation points, the higher the better the resolution
figure, ft_multiplotTFR(cfg,this_tfr)

stat2.raweffect = effectDiff.powspctrm;
cfg = [];
cfg.alpha  = 0.025;
cfg.parameter = 'stat';
cfg.zlim   = [-1e-27 1e-27];
cfgl.layout       = 'EEG1005.lay';
cfg.layout        = ft_prepare_layout(cfgl); clear cfgl;                    % prepares the layout according to the definition before
ft_clusterplot(cfg, stat2);


%% plot all conditions
cfg = [];
effectCTRL = ft_freqgrandaverage(cfg, structTFR{1:20});
effectET = ft_freqgrandaverage(cfg, structTFR{21:41});

cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = '(x1-x2)/(x1+x2)';
effectDiff = ft_math(cfg, effectET, effectCTRL);

cfgc = [];
cfgc.masknans           = 'yes';
cfgl.layout             = 'EEG1005.lay';
cfgc.layout             = ft_prepare_layout(cfgl); clear cfgl;          % prepares the layout according to the definition before
cfgc.colormap           = othercolor('David1');                         % defines the colormap to use, see also omikron/othercolor
cfgc.gridscale          = 40;                                           % number of interpolation points, the higher the better the resolution
cfgc.interplimits       = 'head';                                       % interpolates data to the edges of the head
cfgc.xlim               = [-0.1 : 1.2];                            % as only one frequency given, this value is 1. In cases of frequency vector this enables topoplots at different FOI
cfgc.ylim               = [12 90];
%cfgc.zlim               = [-2 2];
cfgc.layout.pos(:,1) = cfgc.layout.pos(:,1)*1.07;                   % the next two lines shift the values so that
cfgc.layout.pos(:,2) = cfgc.layout.pos(:,2)*1.2;                    % the interpolation of results is minimised
cfgc.baseline = [-.4 0];
cfgc.baselinetype = 'relchange';

ft_multiplotTFR(cfgc, effectET)

ft_multiplotTFR(cfgc, effectCTRL)

effectDiff.mask = stat.mask;
cfgc.mask = 'mask';
ft_multiplotTFR(cfgc, effectDiff)
%
%
%    build configuration
cfg = [];
cfg.method = 'analytic'; %% do a mass-univariate test
cfg.alpha = 0.05; %% critical value around ± 2.09
cfg.statistic = 'indepsamplesT'; % use a dependent samples t-test
cfg.design(1, :) = [1:numel(subj{1}) 1:numel(subj{2})];
cfg.design(2, :) = [ones(1, numel(subj{1})) 2 * ones(1, numel(subj{2}))];
%cfg.uvar = 1; % first row of cfg.design, containing the (u)nits (subjects)
cfg.ivar = 2; % second row of cfg.design, the (i)ndependent events (3&15)
stat = ft_freqstatistics(cfg, structTFR{:});

sig_tiles = find(stat.mask); % find significant time-frequency tiles
[chan freq time] = ind2sub(size(stat.mask),sig_tiles); % transform linear indices to subscript to extract significant channels, timepoints and frequencies
chan = unique(chan);

cfg               = [];
cfg.parameter     = 'stat';
cfg.maskparameter = 'mask';
cfg.maskalpha     = .15; % opacity value for non-significant parts
cfgl.layout             = 'EEG1005.lay';
cfg.layout             = ft_prepare_layout(cfgl); clear cfgl;          % prepares the layout according to the definition before
cfg.colormap           = othercolor('David1');                         % defines the colormap to use, see also omikron/othercolor
cfg.gridscale          = 40;                                           % number of interpolation points, the higher the better the resolution

figure, ft_multiplotTFR(cfg,stat)
%
cfg = [];
cfg.frequency = [stat.freq(1) stat.freq(end)];
cfg.latency   = [stat.time(1) stat.time(end)];

this_tfr = ft_selectdata(cfg, effectDiff);
this_tfr.mask = stat.mask;

cfg = [];
cfg.layout = 'elec1005.lay';
cfg.parameter = 'powspctrm';
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';

ft_multiplotTFR(cfg, this_tfr);
c = colorbar('location', 'southoutside');
c.Label.String = 'Power ratio (right over left)';
title(['Correction: analytical' ]);
%
%
%
%


end


function debugging


cfg = [];
%cfg.channel = 'CPP3h';
cfg.baseline = [-.3 0];
cfg.baselinetype='relchange';
cfg.xlim = [-.3 1.3];
%cfg.ylim = [1 75];
cfg.zlim = [-.5 .5];

iter = 0;
for k = 1:39
    if k < 21; temp = control(subj{1}); %continue; 
    elseif k == 21
        iter = 0;
        temp=patient(subj{2}); 
    else
        temp=patient(subj{2});
    end
    iter =iter+1;
    %if strcmpi(temp(iter).code, 'S43'); keyboard; end
    cfg.title = sprintf('Subject is %s', upper(temp(iter).code));
    figure;
    ft_singleplotTFR(cfg, structTFR{k});
end


end