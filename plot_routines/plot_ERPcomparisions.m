function plot_ERPcomparisions(avg, lgnd, fignum, ch, toi, bsl, nCol, mcp, tit)

%   This function plots differences between two averages; several
%   potential options may be provided

%   Copyright (C) Mai 2021
%   D. Pedrosa, University Hospital of Gießen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

if nargin < 9
    tit = '';
elseif nargin < 8
    mcp = 'cluster_ft';
elseif nargin < 7
    nCol = 3;
elseif nargin < 6
    nCol = 3; bsl = [-.2 0];
elseif nargin < 5
    nCol = 3; bsl = [-.2 0]; toi = [-.2 1];
elseif nargin < 4
    nCol = 3; bsl = [-.2 0]; toi = [-.2 1];
    ch = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz', 'Oz'};
elseif nargin < 3
    nCol = 3; bsl = [-.2 0]; toi = [-.2 1];
    ch = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz', 'Oz'};
    fignum = 1;
elseif nargin < 2
    nCol = 3; bsl = [-.2 0]; toi = [-.2 1];
    ch = {'Fz', 'FCz', 'Cz', 'CPz', 'Pz', 'Oz'};
    fignum = 1;
    lgnd = {'condition1', 'condition2'};
elseif nargin < 1
    fprintf('Please enter some data')
    return
end

if numel(avg)~=2
    fprintf('Two structures with timelockaverages needed!')
    return
end

% Formulae used throughout script
fx_bslsub       = @(x,y,z) bsxfun(@minus, x(y,:), nanmean(x(y,z(1):z(2))));
fx_bslsub_all   = ...
    @(x,z) bsxfun(@minus, x, nanmean(nanmean(x(:,z(1):z(2)),1)));
fac             = 1;                                                        % factor to create SEM for plots
fx_plots = {@(x) fac*nanmean(x), ...
    @(x) fac*(nanmean(x) - 1.96*nanstd(x)./sqrt(size(x,1))), ...
    @(x) fac*(nanmean(x) + 1.96*nanstd(x)./sqrt(size(x,1)))};
p = figure_params_gen;

% Create metainfo needed to plot results
nRow = ceil( numel(ch) / nCol ) ;
rowH = 0.7 / nRow ;  colW = 0.65 / nCol ;
colX = 0.06 + linspace( 0, 0.96, nCol+1 ) ;  colX = colX(1:end-1) ;
rowY = 0.1 + linspace( 0.9, 0, nRow+1 ) ;  rowY = rowY(2:end) ;
toi = dsearchn(avg{1}{1}.time', toi');                                      % time of interest for ERP estimation
time_vector = avg{1}{1}.time(toi(1):toi(2));                                       % create a time_vector to make things easier
bsl = dsearchn(avg{1}{1}.time', bsl');                                      % indices for baseline period

% Start plotting results
figure(fignum); clf;                                                        % creates a figure
set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.1,0.1,0.6,0.6] ) ;
all_labels = avg{1}{1}.label;
for dId = 1:numel(ch) % loop through the channels of interest
    rowId = ceil( dId / nCol ); colId = dId - (rowId - 1) * nCol ;
    axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
    chtemp = find(strcmp(all_labels, ch{dId}));
    dat_all = cell(1,2);                                                    % pre-allocate space
    
    for g = 1:2 % loop through both conditions (avg{1}/{2}) cf. lgnd;
        data_temp = arrayfun(@(x) squeeze(avg{g}{x}.trial(:,chtemp,:)), ...
            1:numel(avg{g}), 'Un', 0);                         %#ok<FNDSB>  % extract channel of interest
        dat_all{g} = cat(1, data_temp{:});
        dat_all{g} = fx_bslsub_all(dat_all{g}, bsl); clear data_temp % subtract baseline from all trials
        
        tmp_data_avg = nan(numel(fx_plots), ...
            length(toi(1):toi(2)));                           % pre-allocate space
        for fx = 1:numel(fx_plots) % loop through different metrics
            tmp_data_avg(fx,:) = ...
                fac * fx_plots{fx}(dat_all{g}(:,toi(1):toi(2)));
        end
        
        % Line plots and fill out in between to show Confidence interval
        m(g) = plot(time_vector, tmp_data_avg(1,:), ...
            'Color', p.colors{g+2}); hold on;
        fillx = [time_vector, fliplr(time_vector)];
        filly = [tmp_data_avg(2,:), fliplr(tmp_data_avg(3,:))];
        f(g) = fill(fillx, filly, [p.colors{g+2}], ...
            'EdgeColor', 'none', 'FaceAlpha', .2);
        mygca(dId) = gca;
        set(gca, 'FontName', p.ftname, 'FontSize', p.ftsize(2))
        set(mygca(dId), 'ydir', 'reverse')
        % Start with statistics toi compare conditions
        if g == 2
            [h, mask, ps] = typeI_correction(dat_all, toi, mcp);

            switch mcp
                case 'none_old'                                                 % simple t-test no multiple comparison correction
                    [h,ps] = ttest2(dat_all{1}(:,toi(1):toi(2)), ...
                        dat_all{2}(:,toi(1):toi(2)), 'Alpha', .01);
                case 'cluster_ft_old'
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
                    mask{dId} = ft_statistics_montecarlo(cfg, ...
                        [dat_all{1}(:,toi(1):toi(2)).', ...
                        dat_all{2}(:,toi(1):toi(2)).'], cfg.design);
                    h = mask{dId}.mask.';
            end
            
            % Highlight significant differences (if present=
            if ~isempty(h) && any(h==1)
                idx_h1 = find(diff(h) == 1);
                idx_h2 = find(diff(h) == -1);
                if size(idx_h2,2) < size(idx_h1,2) && ~isempty(idx_h2)
                    idx_h2(end+1) = length(h)-1;
                elseif size(idx_h1,2) == 1 && isempty(idx_h2)
                    idx_h2 = length(h)-1;
                end
                
                for sig_h = 1:numel(idx_h1)
                    fillx = [time_vector(idx_h1(sig_h):idx_h2(sig_h)), ...
                        fliplr(time_vector(idx_h1(sig_h):idx_h2(sig_h)))];
                    filly = [ones(1,size(fillx,2)/2).*5 ones(1,size(fillx,2)/2).*-5];
                    s(g) = fill(fillx, filly, [p.greys{1}], ...
                        'EdgeColor', 'none', 'FaceAlpha', .1);
                end
            end
        else
            yl_text = get(gca, 'Ylim'); xl_text = get(gca, 'Xlim');
            text(xl_text(2),yl_text(2),sprintf('%s', ch{dId}), ...
                'FontName', p.ftname, 'FontSize', p.ftsize(1), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'BackgroundColor', ones(1,3), ...
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
        legend([m(1), m(end)], lgnd, 'Location', 'SouthEast')
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
text( 0.5, 0, tit, 'FontName', p.ftname, 'FontSize', p.ftsize(2), ...
    'FontWeight', 'Bold', ...
    'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
end

function [h, mask, ps] = typeI_correction(dat, toi, mcp)

%   Subfunction to run statistical tests for comparisons above

%   Copyright (C) Mai 2021
%   D. Pedrosa, Urs Kleinholdermann University Hospital of Gießen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.


switch mcp
    case 'none'                                                             % simple t-test no multiple comparison correction
        [h,ps] = ttest2(dat{1}(:,toi(1):toi(2)), ...
            dat{2}(:,toi(1):toi(2)), 'Alpha', .01);
        mask = [];
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
        cfg.design = [ones(1,size(dat{1},1)), 2*ones(1,size(dat{2},1))];
        cfg.dimord = 'chan_time';
        cfg.dim = [1 size(dat{1}(:,toi(1):toi(2)),2)];
        cfg.neighbours = [];
        cfg.connectivity = 0;
        mask = ft_statistics_montecarlo(cfg, ...
            [dat{1}(:,toi(1):toi(2)).', dat{2}(:,toi(1):toi(2)).'], ...
            cfg.design);
        h = mask.mask.';
        ps = [];
end
end