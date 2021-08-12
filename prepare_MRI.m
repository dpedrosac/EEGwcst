function prepare_MRI(code_parti, filename_erpfinal, mri_file_raw, ...
    mri_file_processed)

%   This function is prepares the 'raw' MRI in order to create a source
%   model and a lead field in later stages

%   ## Version 1.0

%   Copyright (C) August 2021
%   D. Pedrosa
%   University Hospital of Gießen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

%% General settings
[wdir, ~]       = EEGwcst_defaults(0); cd(wdir)
flag_check      = 1;                                                        % defines whether results should be plotted (1) or not (0)

%% Start with ICA in order to remove blink artefacts
fprintf('\n\t loading MRI for subj:\t%s', code_parti)
mri = ft_read_mri(fullfile(wdir, 'raw_MRI', mri_file_raw));                                            % this line loads the raw MRIinto workspace
load(filename_erpfinal)                                         %#ok<LOAD>  % loads data into workspace
try data_final = sorted_data(data_final, 1); catch; end         %#ok<NODEF> % in order to make FT recognize 'Iz', it is necessary to rename it; besides, elec labels are sorted alphabetically for consistency

% Realign data
cfg = [];
cfg.method      = 'interactive';
cfg.coordsys    = 'ctf';
cfg.viewresult  = 'yes';
mri             = ft_volumerealign(cfg, mri);

% Reslice to isotropic voxel size of 1mm
cfg = [];
cfg.resolution  = 1;
cfg.dim         = ones(1,3).*256;
mri             = ft_volumereslice(cfg, mri);

% Write MRI file to disk
cfg = [];
cfg.filename    = mri_file_processed;
cfg.filetype    = 'nifti_gz';
cfg.parameter   = 'anatomy';
ft_volumewrite(cfg, mri);

%  Plot MRI after processing
if flag_check == 1
    ft_sourceplot([], mri)
    keyboard;
end

cfg = [];
cfg.output      = {'brain','skull', 'scalp'};
mri_segmented   = ft_volumesegment(cfg, mri);


end