subj1 = [1, 2, 27, 28, 29, 30, 35, 37, 38, 39, 41, 46, 47, 48, 49, 50];     % CTRL-subjects
subj2 = [3, 4, 5, 7, 8, 10, 14, 18, 19, 20, 22, 23, 24, 42, 44, 45];        % ET-Patients
subj = [subj1, subj2];

%% general settings
wdir = 'C:\Users\David\projekte\wcst\';
metric = {'ERP'};                                                             % two different options available (1) ERP and (2) TFR

for k = 1:numel(metric)
    filename = strcat(metric{k}, '_results.mat');
    
    try
        load(fullfile(wdir, 'results', filename))
    catch
        fprintf('\n%s-results not available, starting the estimation\n', metric{k})
        results_estimation(wdir, subj, metric{k});
    end
end

%%
iter = 0; clear AVG*; avg1 = {}; avg2 = {};
for dId = 1:numel(subj)
    fprintf('\nthe subject being processed is: %s \n', strcat('S', num2str(subj(dId))));
    iter = iter +1;
    filename = strcat(wdir, 'S', num2str(subj(dId)), '\data_final_S', num2str(subj(dId)), '.mat');
    clear data_final
    load(filename);
    trlsnum = 2:3;
    cfg = [];
    cfg.trials = find(ismember(data_final.trialinfo, trlsnum));
    data_early = ft_selectdata(cfg,data_final);
    
    trlsnum = 6:7;
    cfg.trials = find(ismember(data_final.trialinfo, trlsnum));
    data_late = ft_selectdata(cfg,data_final);
    
    trlsnum = [10,20];
    cfg.trials = find(ismember(data_final.trialinfo, trlsnum));
    data_wrong = ft_selectdata(cfg,data_final);
    
    trlsnum = [110,120];
    cfg.trials = find(ismember(data_final.trialinfo, trlsnum));
    data_wrong_alc = ft_selectdata(cfg,data_final);
    
    trlsnum = [21:25];
    cfg.trials = find(ismember(data_final.trialinfo, trlsnum));
    data_right = ft_selectdata(cfg,data_final);
    
    trlsnum = [121:125];
    cfg.trials = find(ismember(data_final.trialinfo, trlsnum));
    data_right_alc = ft_selectdata(cfg,data_final);
    
    trlsnum = 102:103;
    cfg.trials = find(ismember(data_final.trialinfo, trlsnum));
    data_early_alc = ft_selectdata(cfg,data_final);
    
    trlsnum = 106:107;
    cfg.trials = find(ismember(data_final.trialinfo, trlsnum));
    data_late_alc = ft_selectdata(cfg,data_final);
    
    cfg = [];
    cfg.keeptrials = 'yes';
    avgEARLY = ft_timelockanalysis(cfg, data_early);
    avgLATE = ft_timelockanalysis(cfg, data_late);
    avgEARLYa = ft_timelockanalysis(cfg, data_early_alc);
    avgLATEa = ft_timelockanalysis(cfg, data_late_alc);
    avgWRONG = ft_timelockanalysis(cfg, data_wrong);
    avgRIGHT = ft_timelockanalysis(cfg, data_right);
    avgWRONGa = ft_timelockanalysis(cfg, data_wrong_alc);
    avgRIGHTa = ft_timelockanalysis(cfg, data_right_alc);
    
    avg1{iter} = avgEARLY;
    avg2{iter} = avgLATE;
    avg3{iter} = avgEARLYa;
    avg4{iter} = avgLATEa;
    avg5{iter} = avgRIGHT;
    avg6{iter} = avgWRONG;
    avg7{iter} = avgRIGHTa;
    avg8{iter} = avgWRONGa;
    
    clear data_early data_late avgEARLY avgLATE;
end
%%
idx_ctrl = find(ismember(subj, subj1));
idx_et = find(ismember(subj, subj2));
cond = 2;
cfg = [];
% cfg.keepindividual = 'yes';
switch cond
    case (1)
        avg_all1 = ft_timelockgrandaverage(cfg, avg1{:});
        avg_all2 = ft_timelockgrandaverage(cfg, avg2{:});
        avg_all3 = ft_timelockgrandaverage(cfg, avg3{:});
        avg_all4 = ft_timelockgrandaverage(cfg, avg4{:});
        
    case (2)
        avg_all1 = ft_timelockgrandaverage(cfg, avg1{idx_ctrl});
        avg_all2 = ft_timelockgrandaverage(cfg, avg2{idx_ctrl});
        avg_all3 = ft_timelockgrandaverage(cfg, avg1{idx_et});
        avg_all4 = ft_timelockgrandaverage(cfg, avg2{idx_et});
end

%% Replicate Barcelos results
p = figure_params_gen;                                                      % load general parameters for plots
% data_total = {avg_all1, avg_all2, avg_all3, avg_all4};                      % all data
ch = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz', 'POz'};                                      % channels of interest
rows = 2;
sb = reshape(1:numel(ch), [rows, ceil(numel(ch)/rows)]);%[1,2,3,5,7,8,9];
fac = 1;
indiv = 0;
fx_bslsub = @(x,y,z) bsxfun(@minus, x(y,:), nanmean(x(y,z(1):z(2))));
fx_bslsub_all = @(x,z) bsxfun(@minus, x, nanmean(nanmean(x(:,z(1):z(2)),1)));
idx_group = {find(ismember(subj, subj1)), find(ismember(subj, subj2))};
idx_cond = {'avg1_work', 'avg2_work'};
mcp = 'cluster_ft';

opt = 'early-late';
switch opt
    case 'group'
        fc = 0;
        avg1_work = avg1; avg2_work = avg2;
        tit = {'early responses', 'late responses'};
        leg = {'CTRL-subj', 'ET-patients'};
    case 'early-late'
        fc = 10;
        avg1_work = {avg1{idx_group{1}}, avg2{idx_group{1}}};
        avg2_work = {avg1{idx_group{2}}, avg2{idx_group{2}}};
        leg = {'early responses', 'late responses'};
        tit = {'CTRL-subj', 'ET-patients'};
    case 'condition'
        fc = 20;
        avg1_work = {avg5{idx_group{1}}, avg6{idx_group{1}}};
        avg2_work = {avg5{idx_group{2}}, avg6{idx_group{2}}};
        leg = {'right answer', 'wrong answer'};
        tit = {'CTRL-subj', 'ET-patients'};
    case 'condition_group'
        fc = 40;
        avg1_work = avg5; avg2_work = avg6;
        tit = {'right answer', 'wrong answer'};
        leg = {'CTRL-subj', 'ET-patients'};
    case 'alc_group'
        fc = 60;
        avg1_work = {avg6{idx_group{1}}, avg8{idx_group{1}}};
        avg2_work = {avg6{idx_group{2}}, avg8{idx_group{2}}};
        leg = {'wo/ alcohol', 'with alcohol'};
        tit = {'CTRL-subj', 'ET-patients'};
end

% - Compute #rows/cols, dimensions, and positions of lower-left corners.
nCol = 3 ;  nRow = ceil( numel(ch) / nCol ) ;
rowH = 0.7 / nRow ;  colW = 0.65 / nCol ;
colX = 0.06 + linspace( 0, 0.96, nCol+1 ) ;  colX = colX(1:end-1) ;
rowY = 0.1 + linspace( 0.9, 0, nRow+1 ) ;  rowY = rowY(2:end) ;
toi = [-.2 1]; bsl = [-.25 0];
toi = dsearchn(avg1_work{1}.time', toi');                                    % time of interest for ERP estimation
time_vector = avg1{1}.time(toi(1):toi(2));                    % create a time_vector to make things easier

bsl = dsearchn(avg1_work{1}.time', bsl');                                    % indices for baseline period

for fig = 1:2 % loops through early and late trials
    figure(100-fig-fc); clf;
    set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
        'Position', [0.1,0.1,0.6,0.6] ) ;
    avg_dat = eval(idx_cond{fig});
    for dId = 1:numel(ch) % loop through the channels of interest
        rowId = ceil( dId / nCol ); colId = dId - (rowId - 1) * nCol ;
        axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
        
        chtemp = find(strcmp(avg1_work{1}.label, ch{dId}));
        dat_all = cell(1,2);                                                    % pre-allocate space
        
        for g = 1:2 % loop through both groups (1) CTRL,; (2) ET;
            data_temp = arrayfun(@(x) squeeze(avg_dat{x}.trial(:,chtemp,:)), idx_group{g}, 'Un', 0);
            dat_all{g} = cat(1, data_temp{:});
            dat_all{g} = fx_bslsub_all(dat_all{g}, bsl); clear data_temp
            
            clear tmp_data_avg;
            fx_plots = {@(x) fac*nanmean(x), ...
                @(x) fac*(nanmean(x) - 1.96*nanstd(x)./sqrt(size(x,1))), ...
                @(x) fac*(nanmean(x) + 1.96*nanstd(x)./sqrt(size(x,1)))};
            for fx = 1:numel(fx_plots) % loop through different metrics
                tmp_data_avg(fx,:) = fac * fx_plots{fx}(dat_all{g}(:,toi(1):toi(2)));
            end
            m(g) = plot(time_vector, tmp_data_avg(1,:), ...
                'Color', p.colors{g+2}); hold on;
            fillx = [time_vector, fliplr(time_vector)];
            filly = [tmp_data_avg(2,:), fliplr(tmp_data_avg(3,:))];
            f(g) = fill(fillx, filly, [p.colors{g+2}], 'EdgeColor', 'none', 'FaceAlpha', .2);
            mygca(dId) = gca;
            set(gca, 'FontName', p.ftname, 'FontSize', p.ftsize(2))
            if g == 2
                switch mcp
                    case 'none'
                        [h,ps] = ttest2(dat_all{1}(:,toi(1):toi(2)), dat_all{2}(:,toi(1):toi(2)), 'Alpha', .01);
                    case 'cluster_ft'
                        cfg = [];
                        cfg.method = 'montecarlo';
                        cfg.correctm = 'cluster';
                        cfg.neighbours = [];
                        cfg.numrandomization = 1000;
                        cfg.statistic = 'indepsamplesT';
                        cfg.clusterthreshold = 'nonparametric_individual';
                        cfg.alpha = .05;
                        cfg.clusteralpha = .05;
                        cfg.correcttail = 'alpha';
                        cfg.clusterstatistic = 'maxsum';
                        cfg.parameter = 'trial';
                        cfg.ivar = 1;
                        cfg.design = [ones(1,size(dat_all{1},1)), 2*ones(1,size(dat_all{2},1))];
                        cfg.dimord = 'chan_time';
                        cfg.dim = [1 size(dat_all{1}(:,toi(1):toi(2)),2)];
                        cfg.neighbours = [];
                        cfg.connectivity = 0;
                        mask{dId} = ft_statistics_montecarlo(cfg, [dat_all{1}(:,toi(1):toi(2)).', dat_all{2}(:,toi(1):toi(2)).'], cfg.design);
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
                        filly = [ones(1,size(fillx,2)/2).*5 ones(1,size(fillx,2)/2).*-5];
                        s(g) = fill(fillx, filly, [p.greys{1}], 'EdgeColor', 'none', 'FaceAlpha', .1);
                    end
                end
            else
                yl_text = get(gca, 'Ylim'); xl_text = get(gca, 'Xlim');
                text(xl_text(2),yl_text(2),sprintf('%s', ch{dId}), ...
                    'FontName', p.ftname, 'FontSize', p.ftsize(1), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', 'BackgroundColor', ones(1,3), ...
                    'EdgeColor', p.greys{1})
            end
        end
        
        %         switch indiv
        %             case (1)
        %                 colors = {[0, 0, 120]./255, [120, 0, 0]./255};
        %                 for m = 1:numel(avg1)
        %                     if ismember(m,idx_ctrl); c = 1; else; c = 2; end
        %                     plot(dat2.time(toi(1):toi(2)), nanmean(squeeze(avg1{m}.trial(:, chtemp, toi(1):toi(2))),1), 'Color', [colors{c}, .2])
        %                 end
        %         end
        
        plot([time_vector(1) time_vector(end)], [0 0], 'k');
        grid on; box off
        plot([0 0],[-1 1]*6, 'k--')
        
        if dId == numel(ch)
            legend([m(1), m(end)], leg, 'Location', 'SouthEast')
            yl = cell2mat(get(mygca, 'Ylim'));
            ylnew = [min(yl(:,1)) max(yl(:,2))];
            set(mygca, 'Ylim', .75.*ylnew, ...
                'Xlim', [time_vector(1) time_vector(end)])
        end
        
        if rowId == 2
            xlabel('time [in s.]', 'FontName', p.ftname, 'FontSize', p.ftsize(1));
        end
    end
    % - Build title axes and title.
    axes( 'Position', [0, 0.95, 1, 0.05] ) ;
    set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
    text( 0.5, 0, tit{fig}, 'FontName', p.ftname, 'FontSize', p.ftsize(2), 'FontWeight', 'Bold', ...
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
end
%%  plot only Fz and Pz and for the conditions the errors on top and the
%   late responses at the bottom, providing a 2x2 plot

p = figure_params_gen;                                                      % load general parameters for plots

% formulae used in the script:
fx_bslsub = @(x,y,z) bsxfun(@minus, x(y,:), nanmean(x(y,z(1):z(2))));       % baseline subtraction
fx_bslsub_all = @(x,z) bsxfun(@minus, x, nanmean(nanmean(x(:,z(1):z(2)),1)));
fx_plots = {@(x) fac*nanmean(x), ...                                         % plots of mean and SEM
    @(x) fac*(nanmean(x) - 1.96*nanstd(x)./sqrt(size(x,1))), ...
    @(x) fac*(nanmean(x) + 1.96*nanstd(x)./sqrt(size(x,1)))};
idx_group = {find(ismember(subj, subj1)), find(ismember(subj, subj2))};

% settings for the plots
ch          = {'Fz', 'Pz'};                                                 % channels of interest
rows        = 2;                                                            % rows to plot
fac         = 1;                                                            % factor which determines if data is reversed in order to make polarity "the right way round"
indiv       = 0;                                                            % plot individual data (not implemented at the end)
mcp         = 'cluster_ft';                                                 % how to adjust for Type-I error correction
fc          = 10;                                                           % this is subtracted from figure number in order to make it always the same plot
avg1_work   = avg2; avg2_work = avg6;
idx_cond    = {'avg1_work', 'avg2_work'};
tit         = 'Event-related potentials (ERP)';                             % title
leg         = {'CTRL-subj', 'ET-patients'};                                 % legend
sb_title    = {'late responses', 'wrong_trials'};
toi         = [-.2 1];
bsl         = [-.25 -0.05];

% - Compute #rows/cols, dimensions, and positions of lower-left corners.
nCol = 2 ;  nRow = 2 ;
rowH = 0.7 / nRow ;  colW = 0.65 / nCol ;
colX = 0.06 + linspace( 0, 0.96, nCol+1 ) ;  colX = colX(1:end-1) ;
rowY = 0.1 + linspace( 0.9, 0, nRow+1 ) ;  rowY = rowY(2:end) ;

toi = dsearchn(avg1_work{1}.time', toi');                                   % time of interest for ERP estimation
time_vector = avg1{1}.time(toi(1):toi(2));                                  % create a time_vector to make things easier
bsl = dsearchn(avg1_work{1}.time', bsl');                                   % indices for baseline period

figure(100-fc); clf;                                                    % create one figure for all data (2x2 design)
set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.1,0.1,0.6,0.6] ) ;
iter = 0;
for dId = 1:numel(ch) % loop through channels that will be plotted
    for cond = 1:2 % loop through conditions (late_responses vs. wrong_trials)
        iter = iter + 1;
        clear avg_dat; avg_dat = eval(idx_cond{cond});
        rowId = ceil( iter / nCol ); colId = iter - (rowId - 1) * nCol ;
        axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
        chtemp = find(strcmp(avg1_work{1}.label, ch{dId}));                 % find channel index
        dat_all = cell(1,2);                                                % pre-allocate space
        
        for g = 1:2 % loop through both groups (1) CTRL,; (2) ET;
            data_temp = ...                                                 % the next few lines extract only the data for groupm (g)
                arrayfun(@(x) squeeze(avg_dat{x}.trial(:,chtemp,:)), ...
                idx_group{g}, 'Un', 0);
            dat_all{g} = cat(1, data_temp{:});                              % concatenates data to a two-dimensional matrix (trials x time)
            dat_all{g} = fx_bslsub_all(dat_all{g}, bsl);                    % performs baseline subtraction
            
            clear tmp_data_avg data_temp
            for fx = 1:numel(fx_plots) % loop through different metrics
                tmp_data_avg(fx,:) = fac * fx_plots{fx}(dat_all{g}(:,toi(1):toi(2)));
            end
            m(g) = plot(time_vector, tmp_data_avg(1,:), ...
                'Color', p.colors{g+2}); hold on;
            fillx = [time_vector, fliplr(time_vector)];
            filly = [tmp_data_avg(2,:), fliplr(tmp_data_avg(3,:))];
            f(g) = fill(fillx, filly, [p.colors{g+2}], 'EdgeColor', 'none', 'FaceAlpha', .2);
            mygca(iter) = gca;
            set(gca, 'FontName', p.ftname, 'FontSize', p.ftsize(2))
            
            if g == 2
                switch mcp
                    case 'none'
                        [h,ps] = ttest2(dat_all{1}(:,toi(1):toi(2)), dat_all{2}(:,toi(1):toi(2)), 'Alpha', .01);
                    case 'cluster_ft'
                        cfg = [];
                        cfg.method = 'montecarlo';
                        cfg.correctm = 'cluster';
                        cfg.neighbours = [];
                        cfg.numrandomization = 1000;
                        cfg.statistic = 'indepsamplesT';
                        cfg.clusterthreshold = 'nonparametric_individual';
                        cfg.alpha = .05;
                        cfg.clusteralpha = .05;
                        cfg.correcttail = 'alpha';
                        cfg.clusterstatistic = 'maxsum';
                        cfg.parameter = 'trial';
                        cfg.ivar = 1;
                        cfg.design = [ones(1,size(dat_all{1},1)), 2*ones(1,size(dat_all{2},1))];
                        cfg.dimord = 'chan_time';
                        cfg.dim = [1 size(dat_all{1}(:,toi(1):toi(2)),2)];
                        cfg.neighbours = [];
                        cfg.connectivity = 0;
                        mask{dId} = ft_statistics_montecarlo(cfg, [dat_all{1}(:,toi(1):toi(2)).', dat_all{2}(:,toi(1):toi(2)).'], cfg.design);
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
                        filly = [ones(1,size(fillx,2)/2).*5 ones(1,size(fillx,2)/2).*-5];
                        s(g) = fill(fillx, filly, [p.greys{1}], 'EdgeColor', 'none', 'FaceAlpha', .1);
                    end
                end
            else
                yl_text = get(gca, 'Ylim'); xl_text = get(gca, 'Xlim');
                text(xl_text(2),yl_text(2),sprintf('%s', ch{dId}), ...
                    'FontName', p.ftname, 'FontSize', p.ftsize(1), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', 'BackgroundColor', ones(1,3), ...
                    'EdgeColor', p.greys{1})
            end
        end
        plot([time_vector(1) time_vector(end)], [0 0], 'k');
        grid on; box off
        plot([0 0],[-1 1]*6, 'k--')
        
        if iter == numel(ch)*2
            legend([m(1), m(end)], leg, 'Location', 'SouthEast')
            yl = cell2mat(get(mygca, 'Ylim'));
            ylnew = [min(yl(:,1)) max(yl(:,2))];
            set(mygca, 'Ylim', .75.*ylnew, ...
                'Xlim', [time_vector(1) time_vector(end)])
        end
        
        if rowId == 2
            xlabel('time [in s.]', 'FontName', p.ftname, 'FontSize', p.ftsize(1));
        end
    end
    % - Build title axes and title.
    axes( 'Position', [0, 0.95, 1, 0.05] ) ;
    set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
    text( 0.5, 0, tit, 'FontName', p.ftname, 'FontSize', p.ftsize(2), 'FontWeight', 'Bold', ...
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
end

%% topoplots
figure(2);
%formulae and general settings
fx_bslsub_all   = @(x,z) bsxfun(@minus, x, nanmean(x(:,z(1):z(2)),2));
idx_group       = {find(ismember(subj, subj1)), find(ismember(subj, subj2))};
toi             = [.3 .7];
bsl             = [-.25 -0.05];

% Average the timelockstatistics results
avgs = {avg2, avg6};

cfg = []; cfgl = []; cfgc = [];
cfgl.layout             = 'EEG1005.lay';
cfgc.layout             = ft_prepare_layout(cfgl);                  % prepares the layout according to the definition before
cfgc.colormap           = othercolor('David1');                     % defines the colormap to use, see also omikron/othercolor
% cfgc.colorbar           = 'EastOutside';                            % position of the colorbar
cfgc.gridscale          = 40;                                       % number of interpolation points, the higher the better the resolution
cfgc.parameter          = 'mdn';
cfgc.interplimits       = 'head';                                   % interpolates data to the edges of the head
cfgc.xlim               = toi;                                        % as only one frequency given, this value is 1. In cases of frequency vector this enables topoplots at different FOI
cfgc.baseline           = bsl;
cfgc.baselinetype       = 'absolute';
% cfgc.ylim               = 1;                                        % as only one frequency given, this value is 1. In cases of frequency vector this enables topoplots at different FOI
% cfgc.zlim               = [-3 3];
cfgc.highlight          = 'on';                                     % with this options, positive clusters can be highlighted
cfgc.hotkeys            = 'yes';                                    % when enabled, manual change of colorbar is possible (arrowkeys)
cfgc.comment            = 'no';                                     % removes the comments that appear at lower leftcorner in FT topoplots
cfgc.layout.pos(:,1) = cfgc.layout.pos(:,1)*1.07;                   % the next two lines shift the values so that
cfgc.layout.pos(:,2) = cfgc.layout.pos(:,2)*1.2;                    % the interpolation of results is minimised

cfg.keepindividual = 'yes';
for n = 1:2 % loop throughj both groups
    avg_all{n} = eval('ft_timelockgrandaverage(cfg, avg2_work{idx_group{n}})');
    avg_all{n}.avg = squeeze(nanmean(avg_all{n}.individual,1));
    avg_all{n}.mdn = squeeze(nanmedian(avg_all{n}.individual,1));
    
    subplot(1,2,n);
    ft_topoplotER(cfgc,avg_all{n})
end






%%
figure
cfg = [];
cfg.showlabels = 'yes';
cfg.fontsize = 6;
cfg.layout = ft_prepare_layout(avg1{1});
% cfg.ylim = [-3e-13 3e-13];
cfg.xlim = [-0.2 1.4];
cfg.baseline = [-0.3 0];
ft_multiplotER(cfg, avg_all2, avg_all4);

%%
cfg = [];
cfg.xlim = [-0.2:.1:1.2];
cfg.colorbar = 'yes';
cfg.layout = ft_prepare_layout(avg1{1});
cfg.baseline = [-.2 0];
ft_topoplotER(cfg,avg_all4);


%%
timepoint=0.2;
cfg = [];
cfg.zlim='maxmin';
cfg.xlim=[.25 .55];
cfg.layout = 'EEG1005';

figure;
subplot(1,2,1)
ft_topoplotER(cfg,avg_all1)
title ('set shift')
subplot(1,2,2)
ft_topoplotER(cfg,avg_all2)
title ('memory')


%% ft_timelockstatistics

idx_ctrl = find(ismember(subj, subj1));
idx_et = find(ismember(subj, subj2));

avg_stats_ctrl = avg1(idx_ctrl);
avg_stats_et = avg1(idx_et);
fx_sum = @(x) nansum(cell2mat(arrayfun(@(q) size(x{q}.trial,1), 1:numel(x), 'Un', 0)));

cfg = [];
cfg.channel = 'Fz';
cfg.method = 'montecarlo';
cfg.correctm = 'cluster';
cfg.neighbours = [];
cfg.numrandomization = 500;
cfg.statistic = 'indepsamplesT';
cfg.alpha = .05;
cfg.correcttail = 'alpha';
cfg.parameter = 'trial';
cfg.ivar = 1;
cfg.design = [ones(1,fx_sum(avg_stats_et)), 2*ones(1,fx_sum(avg_stats_ctrl))];
stats_ERP = ft_timelockstatistics(cfg, avg_stats_ctrl{:}, avg_stats_et{:});

% look for significant channels in ERP stats
[chans time] = find(stats_ERP.mask)