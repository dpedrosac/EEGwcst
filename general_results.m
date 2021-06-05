function general_results(control, patient, npat, subj1, subj2, paths)

%   This function selects the subjects to be analyzed and provides the
%   demographics for both groups and also the results of the intakte on
%   alcohol on set-shift and memory errors

%   Copyright (C) March 2018
%   D. Pedrosa, University Hospital of Gieﬂen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

if nansum(~cellfun(@isempty, {control(:).matching})) ~= npat || ...         % if the number of matched pairs in the metadata, does not correspond to the number of
        nansum(~cellfun(@isempty, {control(:).matching})) ~= npat           % subjects that should be included, a script is run in which the matching is performed again
    matching(control, patient, paths.data_dir)
    load(fullfile(paths.data_dir, 'patdat.mat')) %#ok<*LOAD>
end

cd(paths.data_dir)
idx = {find(strncmp({patient(:).matching}, 's', 1)), ...
    find(strncmp({control(:).matching}, 's', 1))};                          % indices of patients to which there was a matching (n = 19 * 2, in total)
fx_print = @(x) sprintf('%.1f', x); fx_print2 = @(x) sprintf(' %.1f', x);   % formula needed to get the data into the right format

for g = 1:2 % loop through both groups (1) ET-patients, (2) control subjects
    if g == 1; dattmp = patient(idx{g}); else dattmp = control(idx{g}); end% in order to make the code more efficient, a temporal data structure is defined
    
    age{g} = strcat(fx_print(nanmean([dattmp(:).age])), ' ±',  ...         % summarises age
        fx_print2(nanstd([dattmp(:).age])));
    gndr{g} = strcat(num2str(sum(deal([dattmp(:).gender]))), ':', ...      % summarises gender
        num2str([numel(dattmp) - sum(deal([dattmp(:).gender]))]));
    dd{g} = strcat(fx_print(nanmean([dattmp(:).dd])), ' ±',  ...           % summarises diseas_duration
        fx_print2(nanstd([dattmp(:).dd])));
    for d = 1:numel(dd)                                                     % the next few lines replace the nan's in the control subjects' data with a term for the table 'n.a.'
        if any(ismember(dd{d}, 'NaN')); dd{d} = 'n.a.'; end
    end
    
    rt_all = rtimes_trials(subj1, subj2, paths);                       % uses a different script to extract the response times for the differentv trials;
    % rt_all consists of two cells (CTRL{1} vs. ET{2}) with four columns:
    %   - shift-wo reaction times
    %   - shift-alc reaction times
    %   - memory-wo reaction time
    %   - memory-alc reaction time
    
    errors{1} = {}; errors{2} = {};                                         % pre-allocates space for the error estimation    
    fprintf('estimating response times for group %s \n', num2str(g));
    p = progressbar( numel(dattmp), 'percent' );                            % JSB routine for progress bars
    for k = 1:numel(dattmp) % loop through patients, to get the individual values
        p.update( k )
        err_temp = rt_gen(dattmp(k).code, 'error', paths);             % obtains the error counts from the events-files
        errors = arrayfun(@(x) [errors{x}; err_temp{x}], 1:2, 'Un', 0);     % concatenates the errors for every group in a matrix
    end
    p.stop();
    
    errors = arrayfun(@(x) cat(2,errors{x}{:}), 1:2, 'Un', 0);              % concatenates all errors from both groups to one structure
    parti = unique({errors{1}(:).code}, 'stable');                          % creates a list of all available participants
    
    err_tbl{g} = [];                                           %#ok<*AGROW> % pre-allocates space for
    for k = 1:numel(parti) % loop through all participants
        idx_subj = parti{k};
        idx_all = arrayfun(@(x) find(strcmp({errors{x}(:).code}, idx_subj)), ...
            1:2, 'Un', 0);
        
        err_tbl{g}(k,1:2) = cell2mat(arrayfun(@(x) ...                      % the first two columns provide the anticipations for the WO and the ALC condition
            numel(find(strcmp({errors{x}(idx_all{x}).value}, 'S 40'))), ...
            1:2, 'Un', 0));
        
%         err_tbl{g}(k,3:4) = cell2mat(arrayfun(@(x) ...                      % the third and fourth columns provide the efficient error for the WO and the ALC condition
%             numel(find(strcmp({errors{x}(idx_all{x}).value}, 'S 20'))), ...
%             1:2, 'Un', 0));

        err_tbl{g}(k,3:4) = cell2mat(arrayfun(@(x) ...                      % the third and fourth columns provide the efficient error for the WO and the ALC condition
            numel(find(strcmp({errors{x}(idx_all{x}).value}, 'S 20') | ...
            strcmp({errors{x}(idx_all{x}).value}, 'S 10'))), ...
            1:2, 'Un', 0));

        
        err_tbl{g}(k,5:6) = cell2mat(arrayfun(@(x) ...                      % the fifth and sixth columns provide the set shift error for the WO and the ALC condition
            numel(find(strcmp({errors{x}(idx_all{x}).value}, 'S 50') & ...
            (strcmp({errors{x}(idx_all{x}).prev}, 'S  3') | ...
            strcmp({errors{x}(idx_all{x}).prev}, 'S  4') | ...
            strcmp({errors{x}(idx_all{x}).prev}, 'S  5')))), ...
            1:2, 'Un', 0));
        
        err_tbl{g}(k,7:8) = cell2mat(arrayfun(@(x) ...                      % the seventh and eighth columns provide the errors of holding the rule in memory for the WO and the ALC condition
            numel(find(strcmp({errors{x}(idx_all{x}).value}, 'S 50') & ...
            (strcmp({errors{x}(idx_all{x}).prev}, 'S  6') | ...
            strcmp({errors{x}(idx_all{x}).prev}, 'S  7')))), ...
            1:2, 'Un', 0));
    end
    
    rt_shift1{g} = strcat(fx_print(nanmean(rt_all{g}(:,1))), ' ±',  ...     % response time for shift trials (wo/alcohol)
        fx_print2(nanstd(rt_all{g}(:,1))));
    rt_shift2{g} = strcat(fx_print(nanmean(rt_all{g}(:,2))), ' ±',  ...     % response time for shft trials (w/alcohol)
        fx_print2(nanstd(rt_all{g}(:,2))));
    rt_rest1{g} = strcat(fx_print(nanmean(rt_all{g}(:,3))), ' ±',  ...      % response time for other trials (wo/alcohol)
        fx_print2(nanstd(rt_all{g}(:,3))));
    rt_rest2{g} = strcat(fx_print(nanmean(rt_all{g}(:,4))), ' ±',  ...      % response time for other trials (w/alcohol)
        fx_print2(nanstd(rt_all{g}(:,4))));
    
    ssh1{g} = strcat(fx_print(nanmean(err_tbl{g}(:,5))), ' ±',  ...          % set shifting errors (wo/ condition)
        fx_print2(nanstd(err_tbl{g}(:,5))));
    ssh2{g} = strcat(fx_print(nanmean(err_tbl{g}(:,6))), ' ±',  ...          % set shifting errors (alc/ consition)
        fx_print2(nanstd(err_tbl{g}(:,6))));
    
    mem1{g} = strcat(fx_print(nanmean(err_tbl{g}(:,7))), ' ±',  ...          % "memory errors" (wo/ condition)
        fx_print2(nanstd(err_tbl{g}(:,7))));
    mem2{g} = strcat(fx_print(nanmean(err_tbl{g}(:,8))), ' ±',  ...          % "memory errors" (alc/ condition)
        fx_print2(nanstd(err_tbl{g}(:,8))));
    data{g} = dattmp;                                                       % needed for statistical analyses
end

%% Estimate the significance levels using non-parametric tests
p_gndr = '-';
p_dd = '-';

test = 'unpaired';
switch test % switch to change between paired and unpaired t-tests
    case 'unpaired'
        fx_sig = @(x,y) ranksum(x,y); % Wilcoxon Rank sum test for unpaired samples
        srt = 1:length(data{1});
    case 'paired'
        fx_sig = @(x,y) signrank(x,y); % Signed rank test for paired samples
        for k = 1:numel(data{1}) % loop through all participants
            srt(k) = find(strcmp({data{2}.code}, data{1}(k).matching));
        end
end

p_age = fx_sig([data{1}(:).age].', [data{2}(srt).age].');
p_ssh = fx_sig(err_tbl{1}(:,5), err_tbl{2}(srt,5));                         % next two lines are not needed, because ANOVA is computed
p_mem = fx_sig(err_tbl{1}(:,7), err_tbl{2}(srt,7));

idx_pat = {1:npat, 1:npat};                                                 % indices if ET-patinets and CTRL-subjects; because balanced, a 1:numsubj is introduced
data_raw = arrayfun(@(x) [err_tbl{x}(idx_pat{x},5:8), ...
    reshape(1:length(data{x}(idx_pat{x})), ...
    [length(data{x}(idx_pat{x})),1])], 1:2, 'Un', 0);
plot_effect_alcohol(data_raw, 99, 2)
dat_plot{1} = [rt_all{2}(:,1:2), rt_all{1}(:,1:2)];
dat_plot{2} = [rt_all{2}(:,3:4), rt_all{1}(:,3:4)];
plot_rt_anova(dat_plot, 98, 1);

test_anova = 1;
switch test_anova
    case(1)
        % Prepare data for the ANOVA in Matlab and in R
        fx_tf = @(x) (x);
        num_v = 4;                                                          % number of different variables per subject, that is (1) pers_err_wo, (2) pers_err_alc, (3) mem_err_wo and (4) mem_err_wo
        dat_all = [data_raw{1}; data_raw{2}];                               % concatenates data to one matrix
        limit = length(data_raw{1});                                        % number of subjects per group
        Xs = nan(length(dat_all)*num_v,5);                                  % pre-allocates space for computations later
        
        iter = 0;
        for n = 1:num_v:length(Xs) % loops through all data
            iter = iter + 1;
            Xs(n:n+num_v-1,1) = fx_tf([dat_all(iter,1:num_v)]).';
            if iter <= limit
                Xs(n:n+num_v-1,2) = 1;
            else
                Xs(n:n+num_v-1,2) = 2;
            end
            Xs(n:n+num_v-1,3) = [1;1;2;2];
            Xs(n:n+num_v-1,4) = [1;2;1;2];
            Xs(n:n+num_v-1,5) = ones(num_v,1)*iter;%dat_all(iter,5);
        end
        
        tble = table( Xs(:,1), Xs(:,2), Xs(:,3), Xs(:,4),Xs(:,5), ...
            'VariableNames', {'error_count', 'group', 'error_type', 'condition', 'subj'} );
        writetable(tble, 'anova_data.txt', 'Delimiter', '\t');
end
anova_script(Xs,.05);                                                       % for comparison with R script, the matlab anova sctipt is introduced here

%% Run split-plot ANOVA to test for differences between groups and within subjects foro different conditions
% Create the ANOVA table from the available variables
fx_tf = @(x) (x);

ID = nan(npat*2,1); iter = 0;
for g = 1:2 % loop through both groups
    for k = 1:length(data{g})
        iter = iter + 1;
        idx_tmp = regexp(data{g}(k).code, 's', 'split');
        
        ID(iter,1) = str2num(idx_tmp{2}); clear idx_tmp
    end
end

subj = [1:length(ID)].';
group = [ones(length(data{1}),1);ones(length(data{2}),1)*2];
memserr_wo = fx_tf([err_tbl{1}(:,7); err_tbl{2}(:,7)]);                     % memory errors (WO)
memserr_alc = fx_tf([err_tbl{1}(:,8); err_tbl{2}(:,8)]);                    % memory errors (ALC)
sserr_wo = fx_tf([err_tbl{1}(:,5); err_tbl{2}(:,5)]);                       % set-shift error (WO)
sserr_alc = fx_tf([err_tbl{1}(:,6); err_tbl{2}(:,6)]);                      % set-shift error (ALC)

tbl_anova = table(ID, subj, group, memserr_wo, memserr_alc, sserr_wo, sserr_alc);
writetable(tbl_anova, 'anova_mixedmodel.txt', 'Delimiter', '\t');

%% Plot the table to be introduced in presentation/manuscript etc,
row_names = {'Age [in years]', 'Gender [?:?]', ...
            'Set-shifting error', 'Memory error'};
column_names = {sprintf('ET-patients \n n = %d', numel(patient(idx{1}))), ...
    sprintf('Control subjects \n n = %d', numel(control(idx{2})))};

fx_nx = @(x) num2cell(x);
T = cell2table([{age{1}; gndr{1}; dd{1}; rt_shift1{1}; rt_shift2{1}; rt_rest1{1}; rt_rest2{1}; ssh1{1}; ssh2{1}; mem1{1}; mem2{1}},...
    {age{2}; gndr{2}; dd{2}; rt_shift1{2}; rt_shift2{2};  rt_rest1{2}; rt_rest2{2}; ssh1{2}; ssh2{2}; mem1{2}; mem2{2}},...
    {fx_nx(p_age); p_gndr; p_dd; 'n.a.';  'n.a.';'n.a.'; 'n.a.'; p_ssh; p_ssh; p_mem; p_mem}], ...
    'VariableNames', {'ET', 'ctrl', 'sig'}, 'RowNames', {'Age', 'Gender', ...
    'disease_duration', 'rt_shift_wo', 'rt_shift_alc', 'rt_gen_wo', 'rt_gen_alc', 'set_shifting_wo', 'set_shifting_alc', 'memory_wo', 'memory_alc' })
writetable(T,'demographics.xls', 'WriteRowNames', true);