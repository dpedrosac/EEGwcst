all_subj = [subj1, subj2];
excl_subj = [ 16, 21, 45 ];                                                     % excluded subjects due to technical reasons
subj1proc = subj1(~ismember(subj1, excl_subj));
subj2proc = subj2(~ismember(subj2, excl_subj));

idx_group       = {find(ismember(all_subj, subj1proc)), ...                     % get 'group indices'
    find(ismember(all_subj, subj2proc))};


%% Load data if available otherwise run grandaverage for TFR per subj

load_dir = fullfile(ROOTDIR, 'data', 'tfr_data');                           % Folder to load data from
datEfERR = cell(2,1);
subjs = subj1proc; group = 'ET';
for g = 1:2 % loop over groups
    if g == 2; subjs = subj2proc; group = 'CTRL'; end
    datEfERR{g} = cell(numel(subjs),1); iter = 0;
    
    fprintf('Estimating grandaverages for %s: ...', group)
    
    filename_output = sprintf(fullfile(load_dir, ...
        sprintf('grandaverageEfErr_%s.mat', group)));
    if exist(filename_output, 'file')
        fprintf('\n... already finished; \n ...loading data for %s-group\n', group)
        datEfERR{g} = load(filename_output);
    else
        struct_data2read =  [];                                             % pre-allocate data
        for k = 1:numel(subjs) % loop over CTRL subjects
            iter = iter +1;
            fprintf('\n\t ... adding subj: S%d; (%d of %d) ', subjs(k), ...
                k, numel(subjs));
            filename_subj = sprintf('tfr_S%d_efficientErrorWO.mat', subjs(k));
            if exist(filename_subj, 'file')
                struct_data2read = fullfile(load_dir, filename_subj);
            else
                continue
            end
            
            cfg = [];
            cfg.keepindividual  = 'no';
            cfg.inputfile       = struct_data2read;
            % cfg.outputfile      = filename_output;
            datEfERR{g}{iter} = ft_freqgrandaverage(cfg);
        end
    end
    fprintf('\n\t ... done adding files!\n\n')
    datEfErr = datEfERR{g};
    save(filename_output, 'datEfErr', '-v7.3');
end


%% Show basic TFR plots

% 
% cfg = []; cfgl = [];
% cfgl.layout         ='EEG1005.lay';
% cfg.layout          = ft_prepare_layout(cfgl); clear cfgl;
% cfg.baseline        = [-.15 0];
% cfg.xlim            = [-.2 1];
% cfg.baselinetype    = 'db';
% % cfg.channel         = {'EEG', '-AF7', '-TPP9h'}
% figure; ft_topoplotTFR(cfg, dat{1}.grandavg)
% figure; ft_topoplotTFR(cfg, dat{2}.grandavg)


%% Load data if available otherwise run grandaverage for TFR per subj

datMEM = cell(2,1);
subjs = subj1proc; group = 'ET';
for g = 1:2 % loop over groups
    if g == 2; subjs = subj2proc; group = 'CTRL'; end
    datMEM{g} = cell(numel(subjs),1); iter = 0;
    
    fprintf('Estimating grandaverages for %s: ...', group)
    
    filename_output = sprintf(fullfile(load_dir, ...
        sprintf('grandaverageMem_%s.mat', group)));
    if exist(filename_output, 'file')
        fprintf('\n... already finished; \n ...loading data for %s-group\n', group)
        datMEM{g} = load(filename_output);
    else
        struct_data2read =  [];                                             % pre-allocate data
        for k = 1:numel(subjs) % loop over CTRL subjects
            iter = iter +1;
            fprintf('\n\t ... adding subj: S%d; (%d of %d) ', subjs(k), ...
                k, numel(subjs));
            filename_subj = sprintf('tfr_S%d_memoryWO.mat', subjs(k));
            if exist(filename_subj, 'file')
                struct_data2read = fullfile(load_dir, filename_subj);
            else
                continue
            end
            
            cfg = [];
            cfg.keepindividual  = 'no';
            cfg.inputfile       = struct_data2read;
            % cfg.outputfile      = filename_output;
            datMEM{g}{iter} = ft_freqgrandaverage(cfg);
        end
    end
    fprintf('\n\t ... done adding files!\n\n')
    datMem = datMEM{g};
    save(filename_output, 'datMem', '-v7.3');
end


%% Stats


% Cluster based permutation
cfg = [];
cfg.channel          = 'EEG';
cfg.frequency        = [12 50];
cfg.latency          = [0 1.2];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 1;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, dat{1}.grandavg);

cfg.design(1, :) = [1:size(dat{3}.grandavg.powspctrm,1) 1:size(dat{4}.grandavg.powspctrm,1)];
cfg.design(2, :) = [ones(1, size(dat{3}.grandavg.powspctrm,1)), 2 * ones(1, size(dat{4}.grandavg.powspctrm,1))];

% design = zeros(1,size(freqFIC_planar_cmb.powspctrm,1) + size(freqFC_planar_cmb.powspctrm,1));
% design(1,1:size(structTFR{1}.powspctrm,1)) = 1;
% design(1,(size(freqFIC_planar_cmb.powspctrm,1)+1):(size(freqFIC_planar_cmb.powspctrm,1)+...
%     size(freqFC_planar_cmb.powspctrm,1))) = 2;

cfg.ivar             = 2;

[stat2] = ft_freqstatistics(cfg, dat{3}.grandavg, dat{4}.grandavg);

%% Show basic TFR plots


cfg = []; cfgl = [];
cfgl.layout         ='EEG1005.lay';
cfg.layout          = ft_prepare_layout(cfgl); clear cfgl;
cfg.baseline        = [-.15 0];
cfg.xlim            = [-.2 1];
cfg.baselinetype    = 'db';
% cfg.channel         = {'EEG', '-AF7', '-TPP9h'}
figure; ft_topoplotTFR(cfg, dat{3}.grandavg)
figure; ft_topoplotTFR(cfg, dat{4}.grandavg)

%%
cfg = [];cfgl = [];
cfgl.layout         ='EEG1005.lay';
cfg.layout          = ft_prepare_layout(cfgl); clear cfgl;
cfg.parameter = 'powspctrm';
cfg.frequency = [2 45];
cfg.time = [-.2 1.1];
cfg.baseline        = [-.15 0];
cfg.xlim            = [-.2 1];
cfg.baselinetype    = 'db';
ft_multiplotTFR(cfg, dat{2}.grandavg)

figure;
ft_multiplotTFR(cfg, dat{3}.grandavg)
figure;
ft_multiplotTFR(cfg, dat{4}.grandavg)
