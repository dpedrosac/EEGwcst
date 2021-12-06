function plot_boxplots_erps(data_wo, data_alcohol, index, toi, legend)

%   This function displays boxplots of erp amplitudes at different 
%   conditions (wo/ vs. with alcohol) at a time of interest (toi). Besides,
%   statistics are run at using R for ANOVA

%   Copyright (C) October 2021
%   D. Pedrosa, University Hospital of Gießen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

figure(num_figure); hold on; par = figure_params_gen; kind = [1,3;2,4];
fx_tf = @(x) (x);                                                           % transformation of data if necessary, e.g. to make normality of residuals in ANOVA more probable
fx_vec = @(x) x(:);                                                         % formula to vectorise data
fx_metric = @(x) nanmean(x);
ylbl = 'Error count';
error_bars = 0;
whisker = .015; width = .75;
offset = [-.013, .02, .018, -.015];
maxval = cellfun(@nanmean, data_raw, 'un', 0);
ylims = [];
for k = 1:2
    semsval = bsxfun(@rdivide, nanstd(data_raw{k}), sqrt(size(data_raw{k},1)));
    ylims(k,:) = maxval{k}+semsval;
end
ylims= nanmax(nanmax(ylims));

switch plot_type
    case (1) % scatter
        l = {'without alcohol', 'after alcohol intake', ''};
        for ttype = 1:2 % loop through different tremor types (PT vs. IT)
            scatter(fx_tf(data_raw{2}(:,kind(ttype,1))), ...
                fx_tf(data_raw{2}(:,kind(ttype,2))), par.symbols{ttype+1});
            Rtemp = corrcoef(fx_tf(data_raw{2}(:,kind(ttype,1))), ...
                fx_tf(data_raw{2}(:,kind(ttype,2))));
            R(ttype) = Rtemp (2,1);
            dat_range = fx_tf([min(fx_vec(data_raw{2}(:, kind))) ...         % the next two lines estimates minima/maxima for ylim later
                max(fx_vec(data_raw{2}(:, kind)))]);
            
            xlabel(sprintf('log_{10}(average postural tremor amplitude)'));
            if ttype == 1
                ylabel(sprintf('log_{10}(average intentional tremor amplitude)'));
            end
            set(gca, 'FontName', par.ftname, 'FontSize', par.ftsize(2));
            grid on; box off
            legend(l, 'Location', 'NorthEastOutside', 'FontName', par.ftname, ...
                'FontSize', par.ftsize(2));
        end
        % at this point, regression lines could be plotted to see if alcohol changes direction
        
    case (2) % ANOVA plot to see differences between groups and within condition
        plot_style = {'k', '--k'};                                          % for visualisation postural tremor is a solid and intentional tremor a dashed line
        tit = {'ET-patients', 'control subjects'};
        for n = 1:2 % groups, i.e. (1) ET-patients (2) control subjects
            iter = 0; subplot(1,2,n); hold on;% subplots for ctrl subjects and for ET-patients
            for l = 1:2 % type of error, that is (1) set-shifting or (2) memory error
                for c = 1:2 % condition, i.e (1) wo/ alcohol, (2) with alcohol
                    iter = iter +1;
                    scatter(c+offset(iter),fx_metric(fx_tf(data_raw{n}(:,iter))), ...
                        par.ftsize(3)+20, par.scatter{l});                        % this plots the mean of the data for every group/condition and tremor type
                    switch error_bars
                        case (1)
                            means = fx_metric(fx_tf(data_raw{n}(:,iter))) ;
                            sems = nanstd(fx_tf(data_raw{1}(:,iter)))/sqrt(numel(data_raw{1}(:,iter)));
                            plot([c+offset(iter) c+offset(iter)], ...
                                [means means+sems], plot_style{1}, 'Linewidth', width);
                            plot([c+offset(iter) c+offset(iter)], ...
                                [means means-sems], plot_style{1}, 'Linewidth', width);
                            plot([c-whisker+offset(iter) c+whisker+offset(iter)], means-sems*[1 1], plot_style{1}, 'Linewidth', width);
                            plot([c-whisker+offset(iter) c+whisker+offset(iter)], means+sems*[1 1], plot_style{1}, 'Linewidth', width);
                    end
                    if c == 2                                           % plot the line between the two conditions (wo/ alcohol and with alcohol)
                        switch error_bars
                            case (0)
                                p(iter) = plot([1 2], ...
                                    [fx_metric(fx_tf(data_raw{n}(:,iter-1))) ...
                                    fx_metric(fx_tf((data_raw{n}(:,iter))))], ...
                                    plot_style{l}, 'LineWidth', par.lnsize(2));
                                
                            case (1)
                                p(iter) = plot([1+offset(iter-1) 2+offset(iter)], ...
                                    [fx_metric(fx_tf(data_raw{n}(:,iter-1))) ...
                                    fx_metric(fx_tf((data_raw{n}(:,iter))))], ...
                                    plot_style{l}, 'LineWidth', par.lnsize(2));
                        end
                    end
                end
                set(gca, 'XTick', [1 2], 'XTickLabel', ...              % settings for the plot
                    {'wo/ alcohol', 'alcohol'}, 'FontName', par.ftname, ...
                    'FontSize', par.ftsize(2), 'TickLength', [0 0], ...
                    'Ylim', [0 1.05*ylims])
            end
            xlim([.5 2.5]); grid on; %ylim([3.5 ylims])                   % limits, grid
            title(sprintf(tit{n}));                                     % title for every subplot
            if n == 1
                ylabel(ylbl, 'FontName', par.ftname, ...
                    'FontSize', par.ftsize(2))                          % ylabel only for first subplot
            else
                legend([p(2) p(4)], 'set-shifting errors', ...
                    'memory errors')                                % legend, only for second subplot
            end
        end
end
