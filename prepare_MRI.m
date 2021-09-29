function prepare_MRI(code_parti, filename_erpfinal, mri_file_raw, ...
    mri_file_processed)

%   This function is prepares the 'raw' MRI in order to create a source
%   model and a lead field in later stages; everything is aligned to the
%   MNI standard and segmentation is performed according to the default

%   ## Version 1.1

%   Copyright (C) August 2021, modified September 2021
%   D. Pedrosa
%   University Hospital of Gießen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

%% General settings
[wdir, ~]       = EEGwcst_defaults(0); cd(wdir)
flag_check      = 0;                                                        % defines whether results should be plotted (1) or not (0)

%% Preprocess MRI for every subject
fprintf('\n\t loading MRI for subj:\t%s\n', code_parti)
mri = ft_read_mri(fullfile(wdir, 'raw_MRI', mri_file_raw));                 % loads the raw MRI into workspace
load(filename_erpfinal)                                         %#ok<LOAD>  % loads ERP data into workspace
try data_final = sorted_data(data_final, 1); catch; end         %#ok<NODEF> % in order to make FT recognize 'Iz', it is necessary to rename it; besides, elec labels are sorted alphabetically for consistency

template = ft_read_mri(fullfile(wdir, 'templateMRI', ...
    'mni_icbm152_t1_tal_nlin_sym_09c.nii'));                                % template MRI

fsldir = '/opt/fsl/'; fsldirmpath = sprintf('%s/etc/matlab',fsldir);
setenv('FSLDIR', fsldir); setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
path(path, fsldirmpath); clear fsldir fsldirmpath;

% Realign data to default template (in acpc space)
cfg = [];
cfg.method              = 'fsl';
cfg.fsl.path            = '/opt/fsl/bin';
cfg.coordsys            = 'acpc';
cfg.fsl.interpmethod    = 'trilinear';
mri                     = ft_volumerealign(cfg, mri, template);
mri.coordsys = 'acpc';

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

%  Segment MRI to different compartments
cfg = [];
cfg.output      = {'gray','white','csf','skull','scalp'};
mri_segmented   = ft_volumesegment(cfg, mri);

if flag_check == 1 %% NOT WORKING!?!
    seg_i = ft_datatype_segmentation(mri_segmented,'segmentationstyle','indexed');
    
    cfg              = [];
    cfg.funparameter = 'seg';
    cfg.funcolormap  = lines(6); % distinct color per tissue
    cfg.location     = 'center';
    %cfg.atlas        = seg_i;    % the segmentation can also be used as atlas
    ft_sourceplot(cfg, seg_i);
    
end

filename_segmented = fullfile(wdir, 'mri_preprocessed', ...
    sprintf('segmentedMRI_%s.mat', code_parti));

save(filename_segmented, 'mri_segmented', '-v7.3');

end