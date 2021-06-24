function filter_rawdata(filename_clean, filename_filtered, hpf, lpf, wdir)

%   This function filters the data according to the settings provided as a
%   step of the preprocessing routine

%   Copyright (C) June 2021
%   D. Pedrosa, University Hospital of Gießen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

load(filename_clean);                                            %#ok<LOAD> % this line loads the cleaned data into workspace (for details see clean_data.m)
data_clean.elec = ft_read_sens([wdir 'brainamp.sfp']);                      % this block assigns the standard electrode positions to the data
<<<<<<< HEAD
[data_clean.label,I] = sort(data_clean.label);
data_clean.trial{1,1} = data_clean.trial{1,1}(I,:);
=======
data_clean = sorted_data(data_clean, 1);                                    % function to sort EEG channels alphabetically
>>>>>>> beab179d0fc15269759e57e1d637962ca00d14c4

% Do the pre-processing in two steps:
%   (1) preprocess for the ERP estimation,
%   (2) preprocess for the TRF computation

cfg = [];
cfg.reref       = 'yes';                                                    % the next few lines are intended to create average-referenced data
cfg.refchannel  = 'EEG';                                                    % to average reference data, the refchannel is set to 'all'
cfg.detrend     = 'yes';
cfg.bpfilter    = 'yes';                                                    % band-pass filter between .1 and 30Hz
cfg.bpfilttype  = 'firws';                                                  % Filter-type: FIR (window-scaled), for visualisation set cfg.plotfiltresp to 'yes'
cfg.bpfreq      = [hpf(1) lpf];                                             % filter frequencies
cfg.plotfiltresp= 'no';
data_preproc{1} = ft_preprocessing(cfg, data_clean);

cfg = [];
cfg.reref       = 'yes';                                                    % the next few lines are intended to create average-referenced data
cfg.refchannel  = 'EEG';                                                    % to average reference data, the refchannel is set to 'all'
cfg.detrend     = 'yes';
cfg.hpfilter    = 'yes';                                                    % in a second step, data is filtered in order to make TFR estimations possible, therefore, only hpf is used (lpf may be applied later)
cfg.hpfilttype  = 'firws';                                                  % Filter-type: FIR (window-scaled), for visualisation set cfg.plotfiltresp to 'yes'
cfg.hpfreq      = hpf(2);                                                   % cutoff frequency
cfg.plotfiltresp= 'no';
data_preproc{2} = ft_preprocessing(cfg, data_clean);

save(strcat(outdir, filename_filtered), 'data_preproc', '-v7.3');           % saves the preprocessed EEG to a file
