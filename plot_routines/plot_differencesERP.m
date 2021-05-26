function plot_differencesERP(avg)

%   This function creates boxplots between averages at a specific
%   time-of-interest;

%   Copyright (C) Mai 2021
%   D. Pedrosa, University Hospital of Gießen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

ch_sub = {'CPz'};
tit = sprintf('Differences of p300 at %s between groups', cell2mat(ch_sub));
bsl = [-.25 0]; toi = [.25 .45];
fignum = 30;
lgnd1 = {'CTRL-subjects', 'ET-patients'}; % change this to a sprintf including numel(avg)
lgnd2 = {'early', 'late'};

if nargin < 1
    fprintf('Please enter some data')
    return
end

% Formulae and general variables used throughout script
p                   = figure_params_gen;
toi                 = dsearchn(avg{1}{1}.time', toi');                      % time of interest for ERP estimation
bsl                 = dsearchn(avg{1}{1}.time', bsl');                      % indices for baseline period
dat_all             = reformat_data(avg, toi, bsl, ch_sub);                 % returns data cokumwise (see subfunction below)
linewidth_boxplot   = 1.5;
positions           = [1 1.15 1.35 1.5];

%%
% Start plotting results
figure(fignum); clf;                                                        % creates a figure
set( gcf, 'Color', 'White', 'Unit', 'Normalized', ...
    'Position', [0.1,0.1,0.6,0.6] ) ;
boxplot(dat_all, 'Positions', positions)
set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) ])
set(gca,'xticklabel',lgnd1) % Groups on left and right side
hOutliers = findobj(gca,'Tag','Outliers');
delete(hOutliers)

h = findobj(gca,'Tag','Box');
for j=1:numel(h)
    if mod(j,2)~=0; opt = [3, .2]; else; opt = [2, .7]; end
    patch(get(h(j),'XData'), get(h(j),'YData'), p.greys{opt(1)}, ...
        'LineWidth', linewidth_boxplot, 'FaceAlpha',opt(2), ...
        'EdgeColor', 'none');
end

h = findobj(gca,'tag','Median');
set(h,'Color',p.greys{3}, 'LineWidth', 2.5);

h = findobj(gca, 'Tag', 'Upper Whisker');
set(h, 'Color', p.greys{3}, ...
    'LineWidth', linewidth_boxplot, 'Linestyle', '-');
h = findobj(gca, 'Tag', 'Lower Whisker');
set(h, 'Color', p.greys{3}, ...
    'LineWidth', linewidth_boxplot, 'Linestyle', '-');

h = findobj(gca, 'Tag', 'Upper Adjacent Value');
for j = 1:numel(h)
    set(h(j), 'Color', p.greys{3}, ...
        'LineWidth', linewidth_boxplot, ...
        'Linestyle', '-', 'MarkerSize', 1.5);
    xdata = mean(get(h(j),'XData'));
    ywhisker_high_temp = get(h(j),'YData');
    ywhisker_high(j) = ywhisker_high_temp(1);
    xdata = [xdata-.015 xdata+.015];
    set(h(j), 'XData',xdata)
end
ywhisker_high = fliplr(ywhisker_high);

h = findobj(gca, 'Tag', 'Lower Adjacent Value');
for j = 1:numel(h)
    set(h(j), 'Color', p.colors{3}, 'LineWidth', linewidth_boxplot, 'Linestyle', '-', 'MarkerSize', 1.5);
    xdata = mean(get(h(j),'XData'));
    ywhisker_low_temp = get(h(j),'YData');
    ywhisker_low(j) = ywhisker_low_temp(1); %#ok<*AGROW>
    xdata = [xdata-.015 xdata+.015];
    set(h(j), 'XData',xdata)
end
ywhisker_low = fliplr(ywhisker_low);

hold on; grid on; box off

scatter_plot = 1;
if scatter_plot == 1
    for j = 1:numel(avg)
        clear ind
        sizemin_scatter = 10; sizemax_scatter = 20;
        scatter_color = p.colors{3};
        if j <3 
            scatter_color = p.colors{4}; 
        end
        upper_outlier = ywhisker_high(j);%prctile(all_data(:,index{group}(j)), 75) + 1.5*(prctile(all_data(:,index{group}(j)), 75) - prctile(all_data(:,index{group}(j)), 25));
        lower_outlier = ywhisker_low(j);%prctile(all_data(:,index{group}(j)), 25) - 1.5*(prctile(all_data(:,index{group}(j)), 75) - prctile(all_data(:,index{group}(j)), 25));
        
        ind = find(dat_all(:,j) > lower_outlier & ...
            dat_all(:,j) < upper_outlier);
        scatter_data = dat_all(ind,j); %#ok<FNDSB>
        scatter_data = sort(scatter_data);
        if mod(numel(scatter_data),2) ==1
            scatter_data = scatter_data(1:end-1,:);
        end
        
        vector_size = [linspace(sizemin_scatter,sizemax_scatter,numel(ceil(scatter_data))/2).';...
            linspace(sizemax_scatter,sizemin_scatter,(numel(scatter_data))/2).'];
        temp_pos = positions(j);
        x_data_rand =  temp_pos-.05 + (temp_pos+.05-[temp_pos-.05]).*rand(numel(scatter_data),1);
        vector_size_log = [logspace(-1,0,numel(scatter_data)/2),...
            logspace(0,-1,numel(scatter_data)/2)]';
        vector_sign = ones(numel(scatter_data),1);
        vector_sign(x_data_rand < temp_pos) = -1; % what is this for?
        x_data_rand = (abs(x_data_rand-temp_pos).*vector_sign.*vector_size_log)+temp_pos;
        scatter(x_data_rand, scatter_data, vector_size, ...
            scatter_color, 'filled', 'MarkerFaceAlpha',4/8);
        % scatter(j*ones(size(scatter_data,1),1), scatter_data, vector_size, scatter_color, 'filled');
        clear scatter_data
    end
end

% Run simple stats on meaningful combinations
combs = {[1,2], [1,3], [2,4], [3,4]};
offset = 1.05;
for c =1:numel(combs)
    dat_stat1 = dat_all(:,combs{c}(1));
    dat_stat2 = dat_all(:,combs{c}(2));
    [~,ps] = ranksum(dat_stat1, dat_stat2);
    if ps < .01
        x1 = positions(combs{c}(1)) + .01;
        x2 = positions(combs{c}(2)) - .01;
        plot([x1 x2], [max(ywhisker_high)*offset max(ywhisker_high)*offset], ...
            'Color', [84 84 84]./255)
        txt1 = '*';
        h = text(mean(positions(combs{c})),max(ywhisker_high)*offset*1.1,txt1);
        set(h, 'FontName', p.ftname, 'FontSize', p.ftsize(3));
        offset = offset*1.1;
    end
end

% change appearance at the end
c = get(gca, 'Children');
ws = warning('off', 'MATLAB:legend:IgnoringExtraEntries');
hleg1 = legend(c([9,12]), lgnd2, ...
    'FontName', p.ftname, 'FontSize', p.ftsize(1), ...
    'Location', 'NorthEast' );
set(gca, 'FontName', p.ftname, 'FontSize', p.ftsize(2));
title(tit, 'FontSize', p.ftsize(2))
ylabel('ERP [in a.u.]', 'FontSize', p.ftsize(2));

end

function dat_all = reformat_data(avg, toi_idx, bsl_idx, ch_sub)

%   Subfunction processing data the exact same way as ERP as estimated 
%   cf. plot_ERPcomparisons.m

%   Copyright (C) Mai 2021
%   D. Pedrosa, University Hospital of Gießen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

% Formulae used and general variables
fx_nanmean = @(x,y) nanmean(nanmean(x(:,y(1):y(2)),1));
all_labels = avg{1}{1}.label;
chtemp = strcmp(all_labels, ch_sub);

dat_all = nan(max(cell2mat(cellfun(@numel,avg,'uni',0))),numel(avg));       % pre-allocate space
clear tem*
for d = 1:numel(avg) % loop through all groups
    temp = arrayfun(@(x) squeeze(avg{d}{x}.trial(:,chtemp,:)), ...
        1:numel(avg{d}), 'Un', 0);
    dat_concat = cat(1, temp{:});
    mean_dat = nanmean(nanmean(dat_concat(:,bsl_idx(1):bsl_idx(2)),1));
    temp_corr = arrayfun(@(q) bsxfun(@minus, temp{q}, mean_dat), ...
        1:numel(temp), 'Un', 0);
    dat_all(:,d)= cell2mat(arrayfun(@(q) fx_nanmean(temp_corr{q}, ...
        toi_idx), 1:numel(temp_corr), 'Un', 0)).';
end

end
