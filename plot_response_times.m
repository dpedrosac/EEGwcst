% script to plot the data from the time dependent reaction times
% function plot_rt(rt_time, figure_no)
figure(figure_no); hold on; par = figure_params_gen;
whisker = .015; width = 3;
fx_form = @(x) nanmean(x,1);

% nCol = 1; nRow = 1;                                                         % number of columns/rows in the final plot
% rowH = 0.8/nRow; colW = .8/nCol ;
% colX = 0.10 + linspace( 0, 0.9, nCol+1 ); colX = colX(1:end-1);
% rowY = 0.10 + linspace( 0.9, 0, nRow+1 ); rowY = rowY(2:end);
% rowId = 1; colId = 1;% plot at second row and sixth column
% axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] );
colors = {[204 59 45]./255, [82 114 153]./255, };                           % specific red and blue for the plots according to color scheme before
loc = {'ET-patients', 'control subjects'};
offset = [-.2, .2];                                                                % offset to create overlapping bars
iter = 0;

for g = 1:2
    rt_times_both{1,g}{1}(rt_times_both{1,g}{1} == 0) = nan;                                     % replaces zeros by nan, so that they may be plotted
end
rt_times = {rt_times_both{1,1}{2}, rt_times_both{1,2}{2}};
for bin = 1:size(rt_times{g},2) % every available bin is plotted
    for g = 2:-1:1 % groups, that is ET-patients and control subjects
        if any(isfinite(rt_times{g}(:,bin)))
            iter = iter +1;
            means{g,bin} = fx_form(rt_times{g}(:,bin));
            sems{g,bin} = ( (nanstd(rt_times{g}(:,bin))) ./ sqrt(sum(~isnan(rt_times{g}(:,bin)))) ) * 1.96;
            
            h(g) = bar(iter+offset(g), means{g,bin}, ...
                'FaceColor', par.colors{g}, 'FaceAlpha', .84); hold on;
            set(h(g), 'BarWidth', .8, 'FaceColor', par.colors{g});          % change bar-width and colors
            set(gca, 'FontName', par.ftname, 'FontSize', par.ftsize(2));
            plot([iter+offset(g),iter+offset(g)], [means{g,bin}-sems{g,bin}, means{g,bin}+sems{g,bin}], '-k', 'LineWidth', width);
            plot([iter+offset(g)-whisker iter+offset(g)+whisker], means{g,bin}-sems{g,bin}*[1 1], '-k', 'Linewidth', width);
            plot([iter+offset(g)-whisker iter+offset(g)+whisker], means{g,bin}+sems{g,bin}*[1 1], '-k', 'Linewidth', width);
            
            %         if c == 2
            %             [~,ps(iter),~,stats{iter}] = ...
            %                 ttest2(data_jitter{1}{2}(:,bin), data_jitter{2}{2}(:,bin));% run statistics on both groups of burst counts
            %         end
        else
            continue
        end
%         legend(loc, 'FontSize', par.ftsize(1));
    end
end
box off; grid off; xlim([0 iter+1]);
set(gca, 'Xtick', 1.5:2:iter, 'XminorTick', 'off', 'Box', 'off', ...
        'LineWidth', par.lnsize(2), 'XTickLabel', {'Resp1','Resp2','Resp3','Resp4','Resp5'});
    set(gca, 'FontName', par.ftname, 'FontSize', par.ftsize(2));