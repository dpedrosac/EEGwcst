function plot_topoplots(avg, fignum, lgnd, tit, toi)



%% topoplots

% Formulae and general settings
bsl             = [-.1 0];
use_median      = 0;
% Start plottint the results:
figure(fignum); clf;                                                        % creates a figure
set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.1,0.1,0.6,0.6]) ;

cfg = []; cfgl = []; cfgc = [];
cfgl.layout             = 'EEG1005.lay';
cfgc.layout             = ft_prepare_layout(cfgl); clear cfgl;              % prepares the layout according to the definition before
cfgc.colormap           = othercolor('David1');                     % defines the colormap to use, see also omikron/othercolor
% cfgc.colorbar           = 'EastOutside';                            % position of the colorbar
cfgc.gridscale          = 40;                                       % number of interpolation points, the higher the better the resolution
cfgc.interplimits       = 'head';                                   % interpolates data to the edges of the head
cfgc.xlim               = toi;                                        % as only one frequency given, this value is 1. In cases of frequency vector this enables topoplots at different FOI
cfgc.baseline           = bsl;
cfgc.baselinetype       = 'relchange';
cfgc.zlim               = [-3 3];
cfgc.highlight          = 'on';                                     % with this options, positive clusters can be highlighted
cfgc.hotkeys            = 'yes';                                    % when enabled, manual change of colorbar is possible (arrowkeys)
cfgc.comment            = 'no';                                     % removes the comments that appear at lower leftcorner in FT topoplots
cfgc.layout.pos(:,1) = cfgc.layout.pos(:,1)*1.07;                   % the next two lines shift the values so that
cfgc.layout.pos(:,2) = cfgc.layout.pos(:,2)*1.2;                    % the interpolation of results is minimised

avg_all = cell(1,3);                                                        % Pre allocate space
cfg.method = 'across';

for n = 1:3 % loop through both groups
    if use_median == 1 && n < 3
        mdn = arrayfun(@(q) nanmedian(avg{n}{q}.trial,1), ...
            1:numel(avg{n}), 'Un', 0);
        for k = 1:numel(mdn); avg{n}{k}.mdn= mdn{k}; end
    end
    
    if n == 3 % plot difference
        cfg_math = [];
        cfg_math.operation = 'subtract';
        cfg_math.parameter = 'avg';
        avg_all{n} = ft_math(cfg_math, avg_all{2}, avg_all{1});
    else
        avg_all{n} = ft_timelockgrandaverage(cfg, avg{n}{:});
    end
    
    cfgc.comment = lgnd{n};
    cfgc.commentpos = 'rightbottom';
    subplot(1,3,n);
    ft_topoplotER(cfgc,avg_all{n})
end
sgtitle(tit)


%% Nonparametric tests

% Define neighbours
if ~exist([pwd, '\neighbours.mat'], 'file')
    define_neighbours(avg{1}{1}, pwd, 1, 0)
end
load([pwd, '\neighbours.mat'])


cfg = [];
cfg.channel     = 'EEG';
cfg.neighbours  = neighbours; % defined as above
cfg.latency     = toi;
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_indepsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'cluster';
cfg.correcttail = 'prob';
cfg.numrandomization = 5000;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;

Nsub = numel(avg{1});
% cfg.design(1,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
%cfg.uvar                = 1; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg, avg{1}{:}, avg{2}{:});

%% Plot only significant results of comparison

if (isfield(stat, 'posclusters') || isfield(stat, 'negclusters'))
    figure(fignum+1); clf;                                                  % creates a figure
    set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
        'Position', [0.1,0.1,0.6,0.6]) ;
    
    cfg = [];
    cfg.style               = 'blank';
    cfg.xlim                = toi;                                                    % as only one frequency given, this value is 1. In cases of frequency vector this enables topoplots at different FOI
    cfg.layout              = 'eeg1005.lay';
    cfg.highlight           = 'on';
    cfg.highlightchannel    = find(stat.mask);
    ft_topoplotER(cfg, avg_all{3})
    title('Nonparametric: significant with cluster-based multiple comparison correction') 
end