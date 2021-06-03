function dat_processed = subselect_trials(ft_dat, event_def)

% General settings

if nargin < 2
    fprintf("\nPlease provide some trials to process!")
    return
elseif nargin < 1
    fprintf("\nNo data provided!")
    return
end
event_def = cell2mat(event_def);
dat_processed = cell(1, numel(ft_dat));
warning('off')
p = progressbar( numel(ft_dat), 'percent' );                                % JSB routine for progress bars
for subj = 1:numel(ft_dat)
    dat_temp = ft_dat{subj};
    p.update( subj )
    if isfield(dat_temp, 'avg'); dat_temp = rmfield(dat_temp, 'avg'); end
    idx_trials = find(ismember(dat_temp.trialinfo, event_def));
    cfg = [];
    cfg.trials = idx_trials;
    cfg.feedback = 'no';
    dat_processed{subj} = ft_selectdata(cfg, dat_temp);
    
    cfg = [];
    cfg.bpfiler = 'yes';
    cfg.bpfreq = [.1 30];
    dat_processed{subj} = ft_preprocessing(cfg, dat_processed{subj});
    
    dat_processed{subj}.avg = nanmean(dat_processed{subj}.trial, 1);
end
p.stop()
warning('on')