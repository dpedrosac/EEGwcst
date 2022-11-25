
%%  Control file for analysis of beta at expectation. This file controls 
%   for all the different steps in order to get the preprocessing and 
%   the analysis of the source data
[wdir, ROOTDIR] = EEGwcst_defaults(0);
addpath(fullfile(ROOTDIR, 'othercolor'))


%%
% subj{1} = controls, subj{2} = patients
subj{1} = [1:7, 9:14, 16:22]; %7                                              % all available control subjects for analyses
subj{2} = [1:23];% [1:7, 11:13, 15:23];%[1:3, 5:9, 11:17, 21:23]; %subj13, np=8 and subj16, np=10 must be excluded

for adap = 99
    for n = 1:2
        if n == 1, type = 'c'; else, type = 'p'; end
        switch adap
            case (1)
                read_data(subj, wdir, type)

            case (2)
                clean_data(subj, ROOTDIR, wdir, type)
                
            case (3)
                preprocess_data(subj, ROOTDIR, wdir, type)

            case (4)
                freqstats(subj, ROOTDIR, wdir, type)

            case (5)
                preprocess_sourceanalysis(subj, ROOTDIR, wdir, type)
            
            case (99)
                stats_and_results;
        end
        
    end
end