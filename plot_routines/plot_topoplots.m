function plot_topoplots(avg_topoplot, fignum, lgnd, tit, toi)

%   This function offers three distinct topoplots:
%   i) group 1, ii) group 2 and difference between groups with
%   significance (if present)

%   Copyright (C) July 2021
%   D. Pedrosa, University Hospital of Gießen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

% Formulae and general settings
bsl             = [-.15 0];
use_median      = 0;

% Define settings for all topoplots
cfg = []; cfgl = []; cfgc = []; % entire-cfg, cfg-elec and cfg-topoplot
cfgl.layout             = 'EEG1005.lay';
cfgc.layout             = ft_prepare_layout(cfgl); clear cfgl;              % prepares the layout according to the definition before
cfgc.colormap           = othercolor('David1');                     % defines the colormap to use, see also omikron/othercolor
% cfgc.colorbar           = 'EastOutside';                            % position of the colorbar
cfgc.gridscale          = 40;                                       % number of interpolation points, the higher the better the resolution
cfgc.interplimits       = 'head';                                   % interpolates data to the edges of the head
cfgc.xlim               = toi;                                        % as only one frequency given, this value is 1. In cases of frequency vector this enables topoplots at different FOI
cfgc.baseline           = bsl;
cfgc.baselinetype       = 'absolute'; %'relchange';
cfgc.zlim               = [-3 3];
cfgc.highlight          = 'on';                                     % with this options, positive clusters can be highlighted
cfgc.hotkeys            = 'yes';                                    % when enabled, manual change of colorbar is possible (arrowkeys)
cfgc.comment            = 'no';                                     % removes the comments that appear at lower leftcorner in FT topoplots
cfgc.layout.pos(:,1) = cfgc.layout.pos(:,1)*1.07;                   % the next two lines shift the values so that
cfgc.layout.pos(:,2) = cfgc.layout.pos(:,2)*1.2;                    % the interpolation of results is minimised

% Estimate values for groups and difference between both in a loop of three
avg_all = cell(1,3);                                                        % Pre allocate space
cfg.method = 'across';

%% Display topoplots for both groups and difference between them
for n = 1:3 % loop through both groups and estimate difference at the end
    figure(fignum+n-1); clf;                                                       % creates a figure with fignum as defined in the input of the function
    set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
        'Position', [0.1,0.1,0.6,0.6]) ;
    title_to_plot = sprintf('%s, %s', tit, lgnd{n});
    
    if use_median == 1 && n < 3
        mdn = arrayfun(@(q) nanmedian(avg_topoplot{n}{q}.trial,1), ...
            1:numel(avg_topoplot{n}), 'Un', 0);
        for k = 1:numel(mdn); avg_topoplot{n}{k}.mdn= mdn{k}; end
    end
    
    if n == 3 % plot difference
        cfg_math = [];
        cfg_math.operation = 'subtract';
        cfg_math.parameter = 'avg';
        avg_all{n} = ft_math(cfg_math, avg_all{2}, avg_all{1});
    else
        avg_all{n} = ft_timelockgrandaverage(cfg, avg_topoplot{n}{:});
    end
    
    cfgc.comment = lgnd{n};
    cfgc.commentpos = 'rightbottom';
    ft_topoplotER(cfgc,avg_all{n})
    sgtitle(title_to_plot)
end


%% Nonparametric tests between groups
% Define neighbours
if ~exist([pwd, '\neighbours.mat'], 'file')
    define_neighbours(avg_topoplot{1}{1}, pwd, 1, 0)
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

Nsub1 = numel(avg_topoplot{1});
Nsub2 = numel(avg_topoplot{2});

cfg.design(1,1:[Nsub1+Nsub2])  = [ones(1,Nsub1) 2*ones(1,Nsub2)];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable

stat = ft_timelockstatistics(cfg, avg_topoplot{1}{:}, avg_topoplot{2}{:});

%% Highlight significant results of comparison
if (isfield(stat, 'posclusters') || isfield(stat, 'negclusters'))
    figure(fignum+3); clf;                                                  % creates a figure
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

%% Plot ERP over time

cfgc = [];
cfgl.layout             = 'EEG1005.lay';
cfgc.layout             = ft_prepare_layout(cfgl); clear cfgl;              % prepares the layout according to the definition before
cfgc.colormap           = othercolor('David1');                     % defines the colormap to use, see also omikron/othercolor
cfgc.gridscale          = 40;                                       % number of interpolation points, the higher the better the resolution
cfgc.interplimits       = 'head';                                   % interpolates data to the edges of the head
cfgc.xlim               = [0.1 : 0.1 : 1.0];                                        % as only one frequency given, this value is 1. In cases of frequency vector this enables topoplots at different FOI
cfgc.baseline           = bsl;
cfgc.baselinetype       = 'absolute'; %'relchange';
cfgc.zlim               = [-3 3];
cfgc.layout.pos(:,1) = cfgc.layout.pos(:,1)*1.07;                   % the next two lines shift the values so that
cfgc.layout.pos(:,2) = cfgc.layout.pos(:,2)*1.2;                    % the interpolation of results is minimised

for k = 1:2
    figure(fignum+3+k); clf;                                                       % creates a figure with fignum as defined in the input of the function
    set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
        'Position', [0.1,0.1,0.6,0.6]) ;
    ft_topoplotER(cfgc, avg_all{k});
    sgtitle(sprintf('%s over time, %s', tit, lgnd{k}))
end
