function rtimes_all(subj, subj1, subj2, wdir)

% This function


%%
% general settings
conds = {'WO', 'ALC'};                                                      % two different conditions available
plot_style = {'--k', 'k'};                                                  % for visualisation, ctrl subjects are shown as dashed and ET as solid line
for g = 1:2 % loop through both groups
    for k = 1:7; rt_temp{1,k} = {[], []}; end                                % üre-allocate space for later estimating the RTs
    
    if g == 1; sjts = subj1; else; sjts = subj2; end
    for s = 1:numel(sjts) % loop through subjects
        for c = 1:numel(conds) % loop through the available conditions
            filename = fullfile(wdir, 'wcst', strcat('S', num2str(sjts(s))), ...
                strcat('events_', conds{c}, '_', strcat('S', num2str(sjts(s))), '.mat'));
            load(filename);                                                 % load data for individual subject
            for k = 1:7 % repetitions within a block
                idx = find(strcmp({events.value}, strcat('S', {'  '}, ...   % find the index for the different trials witin a block
                    num2str(k))));
                idx50 = find(cell2mat(arrayfun(@(q) ...
                    strcmp(events(q).value,'S 50'), idx+2, 'Un', 0)));
                idx = setdiff(idx, idx(idx50));
                idx40 = find(cell2mat(arrayfun(@(q) ...
                    strcmp(events(q).value,'S 40'), idx+2, 'Un', 0)));
                idx = setdiff(idx, idx(idx40));
                tmp = [([events(idx+1).sample]./5000 - ...                  % time between disply of the trial S1, S2, ..., S7 until the reaction, that is S31:S34
                    [events(idx).sample]./5000).', ones(length(idx),1) *sjts(s)];
                
                if max(tmp(:,1)) >8 || min(tmp(:,1)) < .3
                    idx_outl = find(tmp(:,1) >8 | min(tmp(:,1) < .2));
                    idx_val = setdiff(1:size(tmp,1), idx_outl);
                    tmp = tmp(idx_val,:);
                end
                rt_temp{1,k}{c} = cat(1, [rt_temp{1,k}{c}], tmp);           % concatenate data to one structure per condition and trial
            end
        end
    end
    
    if g == 1 % CTRL subjects
        rt_ctrl = rt_temp;
    else
        rt_et = rt_temp;
    end
end
%% Plot results
figure(95); hold on;
p = figure_params_gen;                                                      % load general parameters for plots
fx_plots        = {@(x) nanmean(x), ...                                 % plots of mean and SEM
    @(x) (nanmean(x) - 1.96*nanstd(x)./sqrt(size(x,1))), ...
    @(x) (nanmean(x) + 1.96*nanstd(x)./sqrt(size(x,1)))};
fx_outliers = @(x) x(x>.2 & x < 10);

shft = [-.035, .035];
offset = .025;
% colors = {p.colors{4}, p.colors{1}};
leg = {'CTRL-subj.', 'ET-patients'};
for c = 1:2
    subplot(1,2,c); hold on; iter = 0;
    for k = 1:7 % loop through different trials
        for g = 1:2 % loop through both groups
            iter = iter +1;
            if g == 1; rt_all = rt_ctrl; else; rt_all = rt_et; end
            dat = fx_outliers(rt_all{k}{c}(:,1));
            ps(iter) = scatter(k+shft(g), fx_plots{1}(dat), 100, p.scatter{g});
            if k > 1
                ls(iter) = plot([k-1+shft(g) k+shft(g)],[fx_plots{1}(fx_outliers(rt_all{k-1}{c}(:,1))) ...
                    fx_plots{1}(fx_outliers(rt_all{k}{c}(:,1)))], plot_style{g});
                %                'Color', colors{g}, 'LineWidth', p.lnsize(3));
            end
            %plot bars
            plot([shft(g)+k shft(g)+k], [fx_plots{1}(dat) fx_plots{3}(dat)], ...
                '-k', 'LineWidth', p.lnsize(3));
            plot([shft(g)+k shft(g)+k], [fx_plots{2}(dat) fx_plots{1}(dat)], ...
                '-k', 'LineWidth', p.lnsize(3));
            % plot whisker
            plot([k+shft(g)-offset k+shft(g)+offset], [1 1].*fx_plots{2}(dat), ...
                '-k', 'LineWidth', p.lnsize(3))
            plot([k+shft(g)-offset k+shft(g)+offset], [1 1].*fx_plots{3}(dat), ...
                '-k', 'LineWidth', p.lnsize(3))
            
            %         if g == 2
            %             rt_all1 = rt_ctrl;
            %             rt_all2 = rt_et;
            %             [h, psa(1,k)] =  ttest2(fx_outliers(rt_all1{k}{1}(:,1)), fx_outliers(rt_all2{k}{1}(:,1)));
            %         end
            
            set(gca, 'FontName', p.ftname, 'FontSize', p.ftsize(2))
            ylim([1.0 3.5]); xlim([.5 7.5]);
            xlabel(sprintf('Trial within a block'), 'FontSize', p.ftsize(2));
            grid on;
            if c == 1
                ylabel(sprintf('Reaction time [in s.]'))
            end
            if g == 2
                
                %                sjts_meanc = unique(rt_ctrl{k}{1}(:,2));
                %                sjts_meane = unique(rt_et{k}{1}(:,2));
                %
                %                rt_tot = nan(numel(sjts_meanc),2);
                %                for sj = 1:numel(sjts_meane)
                %                    idx1 = find(rt_ctrl{k}{c}(:,2) == sjts_meanc(sj));
                %                    idx2 = find(rt_et{k}{c}(:,2) == sjts_meane(sj));
                %                    rt_tot(sj,:) = [nanmean(rt_ctrl{k}{c}(idx1,1)), ...
                %                        nanmean(rt_et{k}{c}(idx2,1))];
                %                end
                [h,px] = ttest2(rt_ctrl{k}{c}(:,1), rt_et{k}{c}(:,1));
                
                if px < .05/(k*2)
                    text(k,3.2,'*', ...
                        'FontName', p.ftname, 'FontSize', p.ftsize(1), ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'bottom', 'BackgroundColor', ones(1,3), ...
                        'EdgeColor', p.greys{1})
                end
            end
        end
    end
end
legend(ls([3,k+1]), leg)

%% Plot interaction results

figure; clf; clear rtall
plot_style = {'xk', 'ok'};                                                  % for visualisation, ctrl subjects are shown as dashed and ET as solid line

for g = 1:2
    if g == 1; rtall = rt_ctrl; else rtall = rt_et; end
    for k = 1:7
        subplot(1,7,k);
        for c = 1:2
            scatter([c c], [1 1].*nanmean(rtall{k}{c}(:,1)),plot_style{g}); hold on
        end
        xlim([.5 2.5]); ylim([1.2 3.5]);
    end
end


%% Plot differeces in RTs with and without alcohol

figure(94); hold on;
p = figure_params_gen;                                                      % load general parameters for plots
dat_diff = [];
shft = [-.08, .08]; iter = 0;
offset = .13;
colors = {p.colors{4}, p.colors{1}};
leg = {'CTRL-subj.', 'ET-patients'};
for k = 1:7 % loop through different trials
    for g = 1:2 % loop through both groups
        iter = iter +1;
        if g == 1; rt_all = rt_ctrl; else; rt_all = rt_et; end
        sub = unique(rt_all{1}{1}(:,2));
        for s = 1:numel(sub) % loop though different subjects
            idx1 = find(rt_all{k}{1}(:,2)==sub(s) & rt_all{k}{1}(:,1) > .4 & rt_all{k}{1}(:,1) < 8);
            idx2 = find(rt_all{k}{2}(:,2)==sub(s) & rt_all{k}{2}(:,1) > .4 & rt_all{k}{2}(:,1) < 8);
            dat_diff{g}(s,k) = nanmean(rt_all{k}{2}(idx2,1),1) - ...
                nanmean(rt_all{k}{1}(idx1,1),1);
        end
    end
end

% Start plotting

for g = 1:2 % loop through both groups
    for k = 1:7 % loop through the seven trials within blocks
        iter = iter +1;
        dat = dat_diff{g}(:,k);
        ps(iter) = plot([k+shft(g)-offset k+shft(g)+offset], ...
            [1 1].*fx_plots{1}(dat), 'Color', colors{g}, 'LineWidth', p.lnsize(3));
        plot([shft(g)+k shft(g)+k], [fx_plots{1}(dat) fx_plots{3}(dat)], '-k');
        plot([shft(g)+k shft(g)+k], [fx_plots{2}(dat) fx_plots{1}(dat)], '-k');
        set(gca, 'FontName', p.ftname, 'FontSize', p.ftsize(2))
        xlim([.5 7.5]);
        xlabel(sprintf('Trial within a block'), 'FontSize', p.ftsize(2));
        ylabel(sprintf('Difference in reaction times between the w/ alcohol \nand the wo/ alcohol condition [in s.]'))
        
    end
end
grid on;
% legend(ps([1,3]), leg)

%% Put data together as a table
dat = []; iter = 0;

for g = 1:2 % loop through the groups
    if g == 1; rtall = rt_ctrl; else; rtall = rt_et; end
    sjts = unique(rtall{1}{1,c}(:,2));
    for s = 1:numel(sjts)
        iter = iter +1;
        for c = 1:2 % loop through conditions
            for k = 1:7 % loop through trials
                idx = find(rtall{k}{1,c}(:,2) == sjts(s));
                tmp = [sjts(s), iter, g, nanmean(rtall{1,k}{1,c}(idx,1)), c, k];
                dat = [dat; tmp];
            end
        end
    end
end

tbl_anova = table(dat(:,1), dat(:,2),dat(:,3),dat(:,4),dat(:,5),dat(:,6));
tbl_anova.Properties.VariableNames = ...
    {'ID', 'subj', 'group', 'rt', 'cond', 'trial'};
writetable(tbl_anova, strcat(fullfile(wdir, 'wcst', 'results'), ...
    '\anova_rtimes.txt'), 'Delimiter', '\t');


%% Correlation between reaction time and PSP?
trunc = fullfile(wdir, 'wcst\');
file_start = 'ERP_rightWO_';
rt_tot = [];
sjt = subj;
for k = 1:numel(sjt)
    suffix = strcat('S', num2str(sjt(k)));
    filename = strcat(trunc, suffix, '\', file_start, suffix, '.mat');
    load(filename);
    tme = [.25 .45];
    toi = dsearchn(ERP_avg.time', tme');
    ch = find(strcmp(ERP_avg.label, 'Pz'));
    
    idx = find(dat(:,1) == sjt(k) & dat(:,5) == 1 & dat(:,6) >=6);
    tmp = [sjt(k), nanmean(dat(idx,4)), nanmean(nanmean(ERP_avg.trial(:,ch,toi(1):toi(2)),3))];
    rt_tot = [rt_tot; tmp];
end

figure
scatter(rt_tot(:,2), rt_tot(:,3))
