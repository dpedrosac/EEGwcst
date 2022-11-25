function plot_erp_amplitudes_boxplot(data_wo, data_alc, all_subj, ...
    idx_group, channel, toi, file_appendix, ROOTDIR)

%   function used for plotting differences between groups and conditions

%   Version 1.3, 2021-10-13 DP

%   Copyright (C) October 2021
%   D. Pedrosa, University Hospital of Giessen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

if nargin < 8
    [~, ROOTDIR] = EEGwcst_defaults(1);
end

%% Define neighbours
cfg = [];
cfg.method = 'distance';
cfg.feedback = 'no';                                                       % can be set to 'yes' to see results
cfg.neighbourdist = .02;
neighbours = ft_prepare_neighbours(cfg, data_wo{1});

%% Extract data for both groups and all conditions
% TODO get data into long format with subject + condition rowwise and one colum for
% value and run stats: https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/#two-way-repeated-measures-anova
all_labels  = data_wo{1}.label;
ch_oi       = [ channel, neighbours(strcmp(all_labels, channel)).neighblabel ];
toi_idx     = dsearchn(data_wo{1}.time', toi');                             % time of interest for ERP estimation
bsl         = [-.1 0];
bsl_idx     = dsearchn(data_wo{1}.time', bsl');                             % indices for baseline period
iter        = 0;

data_all = nan(numel(idx_group) * sum(cell2mat(cellfun(@numel, ...
    idx_group, 'un', 0))), 4);      % pre-allocate space for iD x group x cond x value
for g = 1:numel(idx_group) % loop over all groups
    data_to_process = data_wo;
    for c = 1:2
        if c == 2; data_to_process = data_alc; end
        for k = 1:numel(idx_group{g}) % loop over all subjects within group
            iter = iter +1;
            data_all(iter, 1) = all_subj(idx_group{g}(k));                      % ID
            data_all(iter, 2) = g;                                              % group
            data_all(iter, 3) = c;                                              % condition
            data_all(iter, end) = ...
                extract_data_ftstruct(data_to_process{idx_group{g}(k)}, ...
                ch_oi, toi_idx, bsl_idx, all_labels);
        end
    end
end
%%
data_table = table(data_all(:,1), data_all(:,2), data_all(:,3), data_all(:,4),...
    'VariableNames', {'ID', 'group', 'condition', 'meanERP'});
writetable(data_table, fullfile(ROOTDIR, 'results', ...
    sprintf('ERP_resuts_per_group_%s_%s.csv', channel, file_appendix)))
data_wide = nan(20, 4); iter = 0; % the next few lines result in a subj x [g1c1 x g1c2 x g2c1 x g2c2] matrix
for g = 1:2
    for k = 1:2
        iter = iter + 1;
        idx_subj = find(data_all(:,2)==g & data_all(:,3)==k);
        dattemp = data_all(idx_subj,4); %#ok<FNDSB>
        data_wide(1:length(dattemp),iter) = dattemp;
    end
end

figure; boxplot(data_wide, 'symbol', '*'); ylim([-5 5]);

end

function value = ...
    extract_data_ftstruct(data_subj, ch_oi, toi_idx, bsl_idx, all_labels)

% subfunction extracting data from fieldtrip structure in order to
% create the boxplots at later stage

fx_bslsub_all   = ...
    @(x,z) bsxfun(@minus, x, mean(x(:,z(1):z(2)),2));
chtemp = find(ismember(all_labels, ch_oi));
data_temp = squeeze(data_subj.avg(:,chtemp,:));                             % extract channels of interest
if size(chtemp,1) == 1; data_temp = data_temp.'; end
dat_all = fx_bslsub_all(data_temp, bsl_idx); clear data_temp                % subtract baseline from all trials

value = mean(mean(dat_all(:, toi_idx(1):toi_idx(2))));                      % returns the average for the subject

end