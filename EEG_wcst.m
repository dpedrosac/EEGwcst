
%%  Control file for analysis of beta at expectation. This file controls 
%   for all the differentsteps in order to get the preprocessing and 
%   the analysis of the source data

<<<<<<< HEAD
[wdir, ROOTDIR] = EEGwcst_defaults(0);
addpath('d:\skripte\othercolor')
=======
restoredefaultpath                                                          % removes confounding software to avoid interaction with fielsdtrip 
addpath('d:\skripte\fieldtrip'); ft_defaults                                        % set feidltrip defaults
addpath(genpath('d:\skripte\lambda')); addpath('d:\skripte\othercolor')                             % change folder to current folder and add a series of colorbars to plot data later
addpath('d:\skripte\lambda\multitaper\m\');                                            % add JSB routines for TFR analyses

%% general settings
wdir = 'd:\skripte\lambda\';                                           % (wdir) defines the working drive for the script
save_dir = strcat(wdir, 'data\');
cd(save_dir);
>>>>>>> beab179d0fc15269759e57e1d637962ca00d14c4

%%
% subj{1} = controls, subj{2} = patients
subj{1} = [1:7, 9:14, 16:22];                                               % all available control subjects for analyses
subj{2} = [1:3, 5:17, 21:23];

for adap = 3
    for n = 1:2
        if n == 1, type = 'c'; else, type = 'p'; end
        switch adap
            case (1)
                read_data(subj, wdir, type)

            case (2)
                clean_data(subj, wdir, type)
                
            case (3)
                preprocess_data(subj, wdir, type)
                
            case (4)
                sensor_analysis(subj, save_dir, type)
                
            case (5)
                sensor_stats(subj, wdir, type)
                
            case (6)
                prepare_source_analysis(subj, save_dir, type)
                
            case (7)
                source_analysis(subj, wdir, type)
                
            case (8)
                source_stats(subj, wdir, type)
                
            case (9)
                correlations_sources(subj, wdir, type)
                
            case (10)
                tremor_analyses(subj, wdir, type);
                
            case (99)
                stats_and_results(subj, wdir, type);
        end
        
    end
end