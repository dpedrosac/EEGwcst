function plotERP(subj, conds, wdir, subj1, subj2, fignum)

%   This function plots the results for theh ERP and the topoplots for the
%   ERPs as the second part of the analysis

%   Copyright (C) March 2018
%   D. Pedrosa, University Hospital of Gieﬂen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

%%  plot only Fz and Pz and for the conditions the errors on top and the
%   late responses at the bottom, providing a 2x2 plot

%% General settings and formulae
p           = figure_params_gen;                                            % load general parameters for plots
p.colors    = fliplr(p.colors([3,4]));
leg         = {'CTRL-subj', 'ET-patients'};                                 % legend for all plots
fac         = 1;                                                            % factor which determines if data is reversed in order to make polarity "the right way round"
toi_time    = [-.2 1.2];                                                    % time of interest in the available data
bsl_time    = [-.2 0];                                                      % baseline definition

% formulae used throughout the script
fx_bslsub_all   = ...
    @(x,z) bsxfun(@minus,x,nanmean(nanmean(x(:,z(1):z(2)),1)));             % formula for baseline correction
fx_plots        = {@(x) fac*nanmean(x), ...                                 % plots of mean and SEM
    @(x) fac*(nanmean(x) - 1.96*nanstd(x)./sqrt(size(x,1))), ...
    @(x) fac*(nanmean(x) + 1.96*nanstd(x)./sqrt(size(x,1)))};

%% extract data from saved individual files:
fprintf('\n extracting data from individual subjects ERP \n');
for f = 1:numel(conds.name) % loops through all filenames and allocates space
    data_plot{f} = [];                                                      % creates an empty matrix in which data may be loaded
    subj_idx{f} = [];                                                       % pre-allocates space for subject index
end

pb = progressbar( numel(subj), 'percent' ); % JSB routine for progress bars
for k = 1:numel(subj) % loop through all available subjects
    pb.update( k )
    suffix = strcat(upper('s'), num2str(subj(k)));
    filename = arrayfun(@(x) fullfile(wdir, 'wcst', suffix, ...             % creates a set of filenames to later load the data
        strcat('ERP_', conds.name{x}, '_', suffix, '.mat')), ...
        1:numel(conds.name), 'Un', 0);
    for m = 1:numel(filename) % loop through the conditions to plot later
        load(filename{m});
        if k == 1 % only do this once for the first subject
            toi = dsearchn(ERP_avg.time', toi_time');                       % time of interest for ERP estimation
            time_vector = ERP_avg.time(toi(1):toi(2));                      % create a time_vector to make things easier
            bsl = dsearchn(ERP_avg.time', bsl_time');                       % indices for baseline period
            for chs = 1:numel(conds.ch) % loop through all channels of interest
                ch_idxs(1,chs) = ...
                    find(ismember(ERP_avg.label, conds.ch{chs}));           % find channel index for all channels listed in (conds.ch)
            end
        end
        
        eval(['avg' num2str(m), '{k} = ERP_avg;']);                         % save the data to one gigangtic structure (necessary for computation of grandaverage later);
        data_plot{m} = [data_plot{m};ERP_avg.trial(:,ch_idxs,:)];           % concatenaes the data to a single matrix of three dimensions (trials x channel x time)
        subj_idx{m} = [subj_idx{m}; ones(size(ERP_avg.trial,1),1)*subj(k)]; % returns an index from which it may be seen data corresponds to which subject
    end
end
pb.stop; clear pb

%% This line is only intended for debugging and to plot the individual
% and the means in order to get a hunch for possible problems (no further
% details are provided!)

flag_check = 0;
if flag_check == 1                                                          % additional sanity check which enables to plot all data generating ERPs to check for outliers
    idx_g = {1:16, 17:32};
    for  m = 1:numel(filename) % loop through conditions
        figure(1000+m); iter = 0;
        for chs = 1:numel(ch_idxs) % loop through channels to plot
            eval(['dat = arrayfun(@(q) squeeze(avg', num2str(m), ...
                '{q}.trial(:,' num2str(ch_idxs(chs)), ...
                ',:)), 1:32, ''Un'', 0);'])
            for g = 1:2 % loop through groups
                iter = iter + 1;
                dtemp = cat(1, dat{idx_g{g}}); %#ok<*USENS>
                dtemp  = fx_bslsub_all(dtemp, bsl);
                subplot(numel(ch_idxs), numel(filename), iter)
                
                plot(dtemp.'); hold on;
                plot(nanmean(dtemp,1).', 'LIneWidth', 5)
                ylim([-1 1].*40)
            end
        end
    end
end

%%  Start plotting data; first, a plot with the two groups, the conds.ch
%   channels and the conds.name trials is created

% - Compute #rows/cols, dimensions, and positions of lower-left corners.
nCol = numel(conds.name);  nRow = numel(conds.ch) ;                         % number of columns and number of rows;
rowH = 0.7 / nRow ;  colW = 0.4 / nCol ;                                    % defines the row-height and the column-width for the subplots
colX = 0.06 + linspace( 0, 0.96, nCol+1 ) ;  colX = colX(1:end-1);          % array defining the settings for the subplots
rowY = 0.1 + linspace( 0.9, 0, nRow+1 ) ;  rowY = rowY(2:end) ;

figure(fignum); clf;                                                        % create one figure for all data (ch x cond design)
set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.1,0.1,0.6,0.6] ) ;
iter = 0;                                                                   % needed later to get the subplots right
for dId = 1:numel(conds.ch) % loop through channels that will be plotted
    for cond = 1:numel(conds.name) % loop through conditions (late_responses vs. wrong_trials)
        iter = iter + 1;
        rowId = ceil( iter / nCol ); colId = iter - (rowId - 1) * nCol ;    % general settings for subplot position and axes of subplots
        axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;        % in these
        
        for g = 1:2 % loop through both groups (1) CTRL,; (2) ET;
            idx_g = ...
                find(ismember(subj_idx{cond}, eval(['subj', num2str(g)]))); % idx for the groups, needed to extract the data from the 3D (data_plot) array
            dat_all{g} = squeeze(data_plot{cond}(idx_g, dId,:));            % concatenates data into a  two-dimensional matrix (trials x time)
            dat_all{g} = fx_bslsub_all(dat_all{g}, bsl);                    % performs baseline subtraction
            
            indiv = 1;                                                      % Different possibility to plot data; instead of taking all data,
            switch indiv                                                    % averages per subject are subjected to plotting routine
                case (0)
                    dat_all2{g} = dat_all{g};
                    cfg = [];                                               % start at this point the configuration for the fieldtrip statistics routine
                case (1) % averages the subjects before plotting
                    npart = unique(eval(['subj', num2str(g)]));             % list all subjects in this group
                    dat_all2{g} = [];                                       % pre-allocate space
                    for l = 1:numel(npart) % loop trough all subjects in the group
                        idx_p = ...
                            find(ismember(subj_idx{cond}, npart(l))); % idx for the groups, needed to extract the data from the 3D (data_plot) array
                        dtemp = ...
                            nanmean(squeeze(data_plot{cond}(idx_p, dId,:)));% concatenates (averaged) data into a  two-dimensional matrix (trials x time)
                        dat_all2{g} = [dat_all2{g}; dtemp];
                        ch2test = [];
                        clear dtemp
                    end
                    dat_all2{g} = fx_bslsub_all(dat_all2{g}, bsl);          % performs baseline subtraction
                    cfg = [];                                               % start at this point the configuration for the fieldtrip statistics routine
            end
            
            clear tmp_data_avg data_temp                                    % clear data to avoid conflicts
            for fx = 1:numel(fx_plots) % loop through different metrics used for mean and SEM plotting (see fx_plots above for details)
                tmp_data_avg(fx,:) = ...
                    fx_plots{fx}(dat_all2{g}(:,toi(1):toi(2))); %#ok<*AGROW>
            end
            m(g) = plot(time_vector, tmp_data_avg(1,:), ...                 % plot average ERP result with a line color depending on the group
                'Color', p.colors{g}); hold on; set(gca, 'XMinorTick', 'on')
            
            fillx = [time_vector, fliplr(time_vector)];                     % the next two lines create the "filling" of the +/- SEM area below and above the mean
            filly = [tmp_data_avg(2,:), fliplr(tmp_data_avg(3,:))];
            f(g) = fill(fillx, filly, [p.colors{g}], ...                    % plots the (shaded) SEM area on top of the mean line
                'EdgeColor', 'none', 'FaceAlpha', .2); mygca(iter) = gca;   % saves the settings for the plot, to modify later
            set(gca, 'FontName', p.ftname, 'FontSize', p.ftsize(2), ...     % changes the axes at this point already (!! some redundancy here, please change to only one modification of the axes using mygca)
                'XMinorTick', 'on')
            
            if g == 2 % when second group is plotted, run stats on th ERP data
                switch conds.mcp.m
                    case 'none' % if no multiple-error correction is desired, only a t-test is run
                        [h,~] = ttest2(dat_all{1}(:,toi(1):toi(2)), ...
                            dat_all{2}(:,toi(1):toi(2)), 'Alpha', .01);
                    case 'cluster_ft' % applying the FT routines, a cluster-based permutation test is applied
                        cfg.neighbours          = [];
                        cfg.method              = conds.mcp.method;
                        cfg.correctm            = conds.mcp.correctm;
                        cfg.numrandomization    = conds.mcp.numrand;
                        cfg.statistic           = conds.mcp.statistic ;
                        cfg.clusterthreshold    = conds.mcp.clusterthreshold;
                        cfg.alpha               = conds.mcp.alpha;
                        cfg.clusteralpha        = conds.mcp.clusteralpha;
                        cfg.correcttail         = conds.mcp.correcttail;
                        cfg.clusterstatistic    = conds.mcp.clusterstatistic;
                        cfg.dimord              = conds.mcp.dimord;
                        cfg.connectivity        = conds.mcp.connectivity;
                        cfg.wcm_weights         = conds.cmp.wcm_weight;

                        
                        cfg.parameter = 'avg';
                        switch indiv
                            case(0)
                                cfg.ivar = 1;
                                cfg.design = [ones(1,size(dat_all{1},1)), ...       % design matrix to put into fieldtrip routine
                                    2*ones(1,size(dat_all{2},1))];
                                cfg.dim = [1 size(dat_all{1}(:,toi(1):toi(2)),2)];  % dimensions of data
                                mask{iter} = ft_statistics_montecarlo(cfg, ...
                                    [dat_all{1}(:,toi(1):toi(2)).', ...
                                    dat_all{2}(:,toi(1):toi(2)).'], cfg.design);
                            case(1)
                                cfgl = []; 
                                cfgl.latency = [toi_time(1) toi_time(2)];
                                cfgl.channel = conds.ch(dId);

                                fx_selectdata = @(x) arrayfun(@(q) ft_selectdata(cfgl, x{q}), 1:numel(x), 'Un', 0);
                                if cond == 1
                                    dat_stats1 = avg1(1:numel(subj1));
                                    dat_stats2 = avg1(numel(subj1)+1:end);                                    
                                else
                                    dat_stats1 = avg2(1:numel(subj1));
                                    dat_stats2 = avg2(numel(subj1)+1:end);                                
                                end
                                dat_stats1 = fx_selectdata(dat_stats1);
                                dat_stats2 = fx_selectdata(dat_stats2);
                                
                                cfg.ivar = 1;
                                cfg.design = [ones(1,numel(dat_stats1)), ...    % design matrix to put into fieldtrip routine
                                    2*ones(1,numel(dat_stats2))];
                                
                                mask{iter} = ...
                                    ft_timelockstatistics(cfg, dat_stats1{:}, dat_stats2{:});
                                clear dats_stats*
                        end
                        h = mask{iter}.mask.';
                end
                
                if ~isempty(h) && any(h==1)                                 % create the indices for the significance statement (grey shaded)
                    idx_h1 = find(diff(h) == 1);
                    idx_h2 = find(diff(h) == -1);
                    idx_h1 = idx_h1(idx_h1>toi(1)& idx_h1<toi(2));
                    idx_h2 = idx_h2(idx_h2>toi(1)& idx_h2<toi(2));
                    
                    if size(idx_h2,2) < size(idx_h1,2) && ~isempty(idx_h2)  % there are several problems with the way
                        idx_h2(end+1) = length(h)-1;
                    elseif size(idx_h1,2) == 1 && isempty(idx_h2)
                        idx_h2 = length(h)-1;
                    elseif size(idx_h2,2) > size(idx_h1,2)
                        idx_h1 = [toi(1), idx_h1];
                    end
                    
                    for sig_h = 1:numel(idx_h1) % loops through the amount of significant areas (shaded)
                        fillx =[time_vector(idx_h1(sig_h):idx_h2(sig_h)),...
                            fliplr(time_vector(idx_h1(sig_h):idx_h2(sig_h)))];
                        filly = [ones(1,size(fillx,2)/2).*6 ...
                            ones(1,size(fillx,2)/2).*-.5];
                        
                        try                                                 % kind of pointless try...catch argument, because code never breaks
                            s(g) = fill(fillx, filly, [p.greys{1}], ...     % but without there is a bug that cannot be understood ...
                                'EdgeColor', 'none', 'FaceAlpha', .1);
                        catch
                            keyboard;
                        end
                    end
                end
            else
                yl_text=get(gca, 'Ylim'); xl_text=[1 1].*max(time_vector);  % positions for markers for channel name
                text(xl_text(2),yl_text(2),sprintf('%s', conds.ch{dId}), ...% plots the channel names to distinguish them in the plot
                    'FontName', p.ftname, 'FontSize', p.ftsize(1), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', ...
                    'BackgroundColor', ones(1,3), 'EdgeColor', p.greys{1})
            end
        end
        plot([time_vector(1) time_vector(end)], [0 0], 'k');
        grid on; box off; plot([0 0],[-1 1]*6, 'k--');
        
        if iter == numel(conds.ch)*2
            legend([m(1), m(end)], leg, 'Location', 'SouthEast')
            yl = cell2mat(get(mygca, 'Ylim'));
            ylnew = [min(yl(:,1)) max(yl(:,2))];
            set(mygca, 'Ylim', .8.*ylnew, ...
                'Xlim', [time_vector(1) time_vector(end)])
        end
        
        ylabel('µs', 'FontName', p.ftname, 'FontSize', p.ftsize(2));
        
        if rowId == 1
            axes( 'Position', [colX(colId), .95, colW, rowH] ) ;
            set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
            text( 0.5, 0, sprintf('%s', conds.sb_title{cond}), 'FontName', p.ftname, 'FontSize', p.ftsize(2), 'FontWeight', 'Bold', ...
                'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
        elseif rowId == numel(conds.ch)
            xlabel('time [in s.]', 'FontName', p.ftname, 'FontSize', p.ftsize(2));
        end
        
    end
end

%% Second plot of differences between conditions
% - Compute #rows/cols, dimensions, and positions of lower-left corners.
nCol = numel(conds.name)-1;  nRow = numel(conds.ch)+1 ;
rowH = 0.7 / nRow ;  colW = 0.32 / nCol ;
colX = 0.06 + linspace( 0, 0.96, nCol+1 ) ;  colX = colX(1:end-1) ;
rowY = 0.1 + linspace( 0.9, 0, nRow+1 ) ;  rowY = rowY(2:end) ;

figure(fignum-1); clf;                                                      % creates the figure according to (fignum)
set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.1,0.1,0.6,0.6] ) ;
iter = 0;

for dId = 1:numel(conds.ch) % loop through channels that will be plotted
    iter = iter + 1;
    rowId = ceil( iter / nCol ); colId = iter - (rowId - 1) * nCol;
    axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] );
    
    for g = 1:2 % loop through both groups (1) CTRL,; (2) ET;
        if g == 1; dattemp = subj1; else dattemp = subj2; end
        dat_all = nan(length(dattemp), size(data_plot{1},3));
        clear dat_subj
        for sj = 1:length(dattemp) % loop through subject in groups
            idx_g = arrayfun(@(q) find(ismember(subj_idx{q}, ...
                dattemp(sj))), 1:numel(subj_idx), 'Un', 0); % idx for subject to be computed
            dat_subj = arrayfun(@(q) squeeze(data_plot{q}(idx_g{q}, dId,:)), ...
                1:numel(subj_idx), 'Un', 0); % concatenates data into a  two-dimensional matrix (trials x time)
            dat_subj = arrayfun(@(q) fx_bslsub_all(dat_subj{q}, bsl), ...
                1:numel(subj_idx), 'Un', 0);                    % performs baseline subtraction
            dat_all(sj,:) = bsxfun(@minus, nanmean(dat_subj{2},1), nanmean(dat_subj{1},1));
        end
        
        clear tmp_data_avg data_temp
        for fx = 1:numel(fx_plots) % loop through different metrics
            tmp_data_avg(fx,:) = ...
                fac * fx_plots{fx}(dat_all(:,toi(1):toi(2))); %#ok<*AGROW>
        end
        m(g) = plot(time_vector, tmp_data_avg(1,:), ...
            'Color', p.colors{g}); hold on;
        set(gca, 'XMinorTick', 'on')
        
        fillx = [time_vector, fliplr(time_vector)];
        filly = [tmp_data_avg(2,:), fliplr(tmp_data_avg(3,:))];
        f(g) = fill(fillx, filly, [p.colors{g}], 'EdgeColor', 'none', 'FaceAlpha', .2);
        mygca(iter) = gca;
        set(gca, 'FontName', p.ftname, 'FontSize', p.ftsize(2), 'XMinorTick', 'on')
        
        if g == 2
            switch conds.mcp.m
                case 'none'
                    [h,~] = ttest2(dat_all(:,toi(1):toi(2)), dat_all(:,toi(1):toi(2)), 'Alpha', .01);
                case 'cluster_ft'
                    cfg = [];
                    cfg.method = 'montecarlo';
                    cfg.correctm = 'cluster';
                    cfg.neighbours = [];
                    cfg.numrandomization = 100;
                    cfg.statistic = 'indepsamplesT';
                    cfg.clusterthreshold = 'nonparametric_individual';
                    cfg.alpha = .05;
                    cfg.clusteralpha = .05;
                    cfg.correcttail = 'alpha';
                    cfg.clusterstatistic = 'maxsum';
                    cfg.parameter = 'trial';
                    cfg.ivar = 1;
                    cfg.design = [ones(1,size(dat_all,1)), 2*ones(1,size(dat_all,1))];
                    cfg.dimord = 'chan_time';
                    cfg.dim = [1 size(dat_all(:,toi(1):toi(2)),2)];
                    cfg.neighbours = [];
                    cfg.connectivity = 0;
                    mask{dId} = ft_statistics_montecarlo(cfg, [dat_all(:,toi(1):toi(2)).', dat_all(:,toi(1):toi(2)).'], cfg.design);
                    h = mask{dId}.mask.';
            end
            
            if ~isempty(h) && any(h==1)
                idx_h1 = find(diff(h) == 1);
                idx_h2 = find(diff(h) == -1);
                if size(idx_h2,2) < size(idx_h1,2) && ~isempty(idx_h2)
                    idx_h2(end+1) = length(h)-1;
                elseif size(idx_h1,2) == 1 && isempty(idx_h2)
                    idx_h2 = length(h)-1;
                end
                
                for sig_h = 1:numel(idx_h1)
                    fillx = [time_vector(idx_h1(sig_h):idx_h2(sig_h)), fliplr(time_vector(idx_h1(sig_h):idx_h2(sig_h)))];
                    filly = [ones(1,size(fillx,2)/2).*50, ones(1,size(fillx,2)/2).*(-.5)];
                    s(g) = fill(fillx, filly, [p.greys{1}], 'EdgeColor', 'none', 'FaceAlpha', .1);
                end
            end
        else
            yl_text = get(gca, 'Ylim'); xl_text = get(gca, 'Xlim');
            text(xl_text(2),yl_text(2),sprintf('%s', conds.ch{dId}), ...
                'FontName', p.ftname, 'FontSize', p.ftsize(1), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'BackgroundColor', ones(1,3), ...
                'EdgeColor', p.greys{1})
        end
    end
    plot([time_vector(1) time_vector(end)], [0 0], 'k');
    grid on; box off
    plot([0 0],[-1 1]*6, 'k--')
    
    if iter == numel(conds.ch)*2
        legend([m(1), m(end)], leg, 'Location', 'SouthEast')
        yl = cell2mat(get(mygca, 'Ylim'));
        ylnew = [min(yl(:,1)) max(yl(:,2))];
        set(mygca, 'Ylim', .75.*ylnew, ...
            'Xlim', [time_vector(1) time_vector(end)])
    end
    
    if rowId == 1
        axes( 'Position', [colX(colId), .95, colW, rowH] ) ;
        set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
        text( 0.5, 0, sprintf('%s', conds.sb_title{cond}), 'FontName', p.ftname, 'FontSize', p.ftsize(2), 'FontWeight', 'Bold', ...
            'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
    elseif rowId == 2
        xlabel('time [in s.]', 'FontName', p.ftname, 'FontSize', p.ftsize(1));
    end
end


%% topoplots
%formulae and general settings
% fx_bslsub_all   = @(x,z) bsxfun(@minus, x, nanmean(x(:,z(1):z(2)),2));
idx_group       = {find(ismember(subj, subj1)), find(ismember(subj, subj2))};
toi             = [.6 .8];
bsl             = [-.15 0];

% Average the timelockstatistics results
% avgs = {avg2, avg6};

cfg = []; cfgl = []; cfgc = [];
cfgl.layout             = 'EEG1005.lay';
cfgc.layout             = ft_prepare_layout(cfgl);                  % prepares the layout according to the definition before
cfgc.colormap           = othercolor('David1');                     % defines the colormap to use, see also omikron/othercolor
cfgc.colorbar           = 'EastOutside';                            % position of the colorbar
cfgc.channel            = {'EEG'};
cfgc.gridscale          = 40;                                       % number of interpolation points, the higher the better the resolution
cfgc.parameter          = 'avg';
cfgc.interplimits       = 'head';                                   % interpolates data to the edges of the head
cfgc.xlim               = toi;                                        % as only one frequency given, this value is 1. In cases of frequency vector this enables topoplots at different FOI
cfgc.baseline           = bsl;
cfgc.baselinetype       = 'db';
% cfgc.ylim               = 1;                                        % as only one frequency given, this value is 1. In cases of frequency vector this enables topoplots at different FOI
cfgc.zlim               = [-1 1].*10;
cfgc.highlight          = 'on';                                     % with this options, positive clusters can be highlighted
cfgc.hotkeys            = 'yes';                                    % when enabled, manual change of colorbar is possible (arrowkeys)
cfgc.comment            = 'no';                                     % removes the comments that appear at lower leftcorner in FT topoplots
cfgc.layout.pos(:,1) = cfgc.layout.pos(:,1)*1.07;                   % the next two lines shift the values so that
cfgc.layout.pos(:,2) = cfgc.layout.pos(:,2)*1.2;                    % the interpolation of results is minimised

figure
cfg.keepindividual = 'yes';
for n = 1:2 % loop throughj both groups
    iter = iter + 1;
    %     rowId = ceil( iter / nCol ); colId = iter - (rowId - 1) * nCol ;
    %     axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
    
    avg_all{n} = eval('ft_timelockgrandaverage(cfg, avg2{idx_group{n}})');
    avg_all{n}.avg = squeeze(nanmean(avg_all{n}.individual,1));
    avg_all{n}.mdn = squeeze(nanmedian(avg_all{n}.individual,1));
    subplot(1,2,n)
    ft_topoplotER(cfgc,avg_all{n})
end

keyboard;