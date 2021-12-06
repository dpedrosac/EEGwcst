function results_1(control, patient, subj1, subj2, ROOTDIR)

%   Estimate the wrong runs, extract reaction times and generate table 1

%   Copyright (C) June 2021, modified July 2021
%   D. Pedrosa, Urs Kleinholdermann University Hospital of Gießen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

%%
cd(ROOTDIR)
idx = {[1:numel(patient)],[1:numel(control)]};

%% Estimate reaction times for all subjects
rt_all = rtimes_trials(subj1, subj2, ROOTDIR);                              % uses a different script to extract the response times for the differentv trials;
rtimes_per_block(subj1, subj2, ROOTDIR)

%%  ==============================================================   %%
% rt_all consists of two cells (CTRL{1} vs. ET{2}) with four columns:
%   - shift-wo reaction times
%   - shift-alc reaction times
%   - memory-wo reaction time
%   - memory-alc reaction time

data = cell(1,2);
for g = 1:2 % loop through groups (1) ET-patients, (2) control subjects
    if g == 1; mta_dat=patient(idx{g}); else mta_dat=control(idx{g}); end   % in order to make the code more efficient, a temporal data structure is defined
    
    %% Select wrong and mark complete trials
    % only complete trials should be analysed according to the definition
    % provided by Barceló et al. 2003, that is trials without
    % anticipations, without perseverations and without memory errors at
    % the end
    wrong_runs{g} = []; % % pre-allocates space for the error estimation
    complete_runs{g} = []; % % pre-allocates space for the error estimation
    pers_err{g} = [];
    fprintf("====\nestimating wrong runs/perseverations for group:%s \n====\n", num2str(g))

    p = progressbar( numel(mta_dat), 'percent' );                           % JSB routine for progress bars
    for k = 1:numel(mta_dat) % subj. per group for individual results
        p.update( k )
        [tmp, c, ~] = select_runs(mta_dat(k).code, ROOTDIR);                   % obtains the error counts from the events-files
        pe = select_perseveration(mta_dat(k).code, ROOTDIR);
        wrong_runs{g} = [wrong_runs{g}; tmp];
        complete_runs{g} = [complete_runs{g}; c];
        pers_err{g} = [pers_err{g}; pe];
    end
    p.stop();
    data{g} = mta_dat;
end

fprintf("====\nWrong runs marked and saved to ./data/header_and_events/events_xx.mat\n====\n")

create_tableOne(data, ...
    arrayfun(@(q) [[rt_all{q}(:,1); rt_all{q}(:,3)], ...
    [rt_all{q}(:,2); rt_all{q}(:,4)]], 1:2, 'Un', 0), ...
    wrong_runs, pers_err, ROOTDIR)
    
end

function table1 = create_tableOne(data, rtimes, incomplete_runs, pers_err, ROOTDIR)

%   This function replaces the code in a previous version which saved
%   data during the main function and moves it into subfunction, to tidy
%   up the entire script
for g = 1:numel(data)
    fx_table1 = @(x) sprintf('%.1f ± %.1f', nanmean(x), nanstd(x));
    age{g} = fx_table1([data{g}(:).age]);
    gndr{g} = strcat(num2str(sum(deal([data{g}(:).gender]))), ':', ...       % summarises gender
        num2str([numel(data{g}) - sum(deal([data{g}(:).gender]))]));
    dd{g} = fx_table1([data{g}(:).dd]);
    for d = 1:numel(dd)                                                     % the next few lines replace the nan's in the control subjects' data with a term for the table 'n.a.'
        if any(ismember(dd{d}, 'NaN')); dd{d} = 'n.a.'; end
    end
    rt_total_wo{g} = fx_table1(rtimes{g}(:,1)); %#ok<*AGROW>
    rt_total_alc{g} = fx_table1(rtimes{g}(:,2));

    completed_wo{g} = fx_table1(incomplete_runs{g}(:,1));
    completed_alc{g} = fx_table1(incomplete_runs{g}(:,2));

    perseveration_wo{g} = fx_table1(pers_err{g}(:,1));
    perseveration_alc{g} = fx_table1(pers_err{g}(:,2));

end

p_gndr = '-';
p_dd = '-';

test = 'unpaired';
switch test % switch to change between paired and unpaired t-tests
    case 'unpaired'
        fx_sig = @(x,y) ranksum(x,y); % Wilcoxon Rank sum test for unpaired samples
        srt = 1:length(data{1});
    case 'paired' % deprecated
        fx_sig = @(x,y) signrank(x,y); % Signed rank test for paired samples
        for k = 1:numel(data{1}) % loop through all participants
            srt(k) = find(strcmp({data{2}.code}, data{1}(k).matching));
        end
end

p_age = fx_sig([data{1}(:).age].', [data{2}(srt).age].');
p_rt_wo = fx_sig(rtimes{1}(:,1), rtimes{2}(:,2));
p_rt_alc = fx_sig(rtimes{1}(:,2), rtimes{2}(:,2));

p_complete_wo = fx_sig(incomplete_runs{1}(srt,1), incomplete_runs{2}(srt,1));
p_complete_alc = fx_sig(incomplete_runs{1}(srt,2), incomplete_runs{2}(srt,2));

p_pers_err_wo = fx_sig(pers_err{1}(srt,1), pers_err{2}(srt,1));
p_pers_err_alc = fx_sig(pers_err{1}(srt,2), pers_err{2}(srt,2));

%% Plot the table to be introduced in presentation/manuscript etc,
% row_names = {'Age [in years]', 'Gender [?:?]', ...
%     'Reaction time', 'Completed trials'};
% column_names = {sprintf('ET-patients \n n = %d', numel(data{1})), ...
%     sprintf('Control subjects \n n = %d', numel(data{2}))};

fx_nx = @(x) num2cell(x);
T = cell2table([{age{1}; gndr{1}; dd{1}; rt_total_wo{1}; rt_total_alc{1}; completed_wo{1}; completed_alc{1}; perseveration_wo{1}; perseveration_alc{1}},...
    {age{2}; gndr{2}; dd{2}; rt_total_wo{2}; rt_total_alc{2}; completed_wo{2}; completed_alc{2}; perseveration_wo{2}; perseveration_alc{2}},...
    {fx_nx(p_age); p_gndr; p_dd; p_rt_wo;  p_rt_alc; p_complete_wo; p_complete_alc; p_pers_err_wo; p_pers_err_alc}], ...
    'VariableNames', {'ET', 'ctrl', 'sig'}, 'RowNames', {'Age', 'Gender', ...
    'disease_duration', 'reaction_times_wo', 'reaction_times_alc', 'incomplete_wo', 'incompleted_alc', 'pers_err_wo', 'pers_err_alc'});
filename_save = fullfile(ROOTDIR, 'results', 'demographics.xls');
writetable(T,filename_save, 'WriteRowNames', true);
disp(T)

end
