function prepare_MRI(code_parti, filename_erpfinal, mri_file_raw)

%   This function is prepares the 'raw' MRI in order to create a source
%   model and a lead field in later stages; everything is aligned to the
%   MNI standard and segmentation is performed according to the default

%   ## Version 1.3 #    complete modification of the pipeline

%   Copyright (C) August 2021, modified September to November 2021
%   D. Pedrosa
%   University Hospital of Gießen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

%% Check https://github.com/meeg-cfin/nemolab/blob/master/basics/nemo_mriproc.m
% steps according to https://www.fieldtriptoolbox.org/workshop/oslo2019/forward_modeling/#head-models-component-1

%% General settings
[wdir, ~]       = EEGwcst_defaults(0); cd(wdir)
flag_check      = 0;                                                        % defines whether results should be plotted (1) or not (0)

%% Preprocess MRI for every subject
fprintf('\n\t loading MRI for subj:\t%s\n', code_parti)
try load(filename_erpfinal); catch; return; end                             % loads ERP data into workspace

version = 1;
switch version
    case 0 % old version using FLS scripts to realign data
        template = ft_read_mri(fullfile(wdir, 'templateMRI', ...
            'MRItrem_template0.nii.gz'));                                           % subject specific template that resulted from ANTs routine
        %    'mni_icbm152_t1_tal_nlin_sym_09c.nii'));                               % template MRI
        
        mri = ft_read_mri(fullfile(wdir, 'raw_MRI', mri_file_raw));         % loads the raw MRI into workspace
        fsldir = '/opt/fsl/'; fsldirmpath = sprintf('%s/etc/matlab', fsldir);
        setenv('FSLDIR', fsldir); setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
        path(path, fsldirmpath); clear fsldir fsldirmpath;
        
        % Realign data to default template (in lps space, i.e. DICOM space)
        cfg = [];
        cfg.method              = 'fsl';
        cfg.fsl.path            = '/opt/fsl/bin';
        cfg.coordsys            = 'lps';
        cfg.fsl.interpmethod    = 'trilinear';
        mri                     = ft_volumerealign(cfg, mri, template);
        mri.coordsys = 'lps';
        
        % Save data to disk
        cfg = [];
        cfg.filename    = mri_file_processed;
        cfg.filetype    = 'nifti_gz';
        cfg.parameter   = 'anatomy';
        ft_volumewrite(cfg, mri);
        
    case(1)
        filename_mri_warped = fullfile(wdir, 'mri_preprocessed', ...
            sprintf('MRIprocessed_%s_Warped.nii.gz', upper(code_parti)));
        try
            mri = ft_read_mri(filename_mri_warped);
        catch
            fprintf("\nMRI for subj: %s not existing, will use template instead\n", upper(code_parti))
            check_for_template(wdir)
            return % in case templates are necessary, this should be performed once, so that a return command is set
        end
end
mri.coordsys = 'lps';

%%  Segment MRI to different compartments
fprintf("\n\tStep 1: Segmenting MRI with different compartments for subj: %s\n", upper(code_parti))
filename_segmented = fullfile(wdir, 'mri_preprocessed', ...
    sprintf('segmentedMRI_%s.mat', code_parti));

if ~exist(filename_segmented, 'file')
    fprintf("\n\t ... Processing ... \t\t\n")
    cfg = [];
    cfg.output      = {'gray','white','csf','skull','scalp'};
    mri_segmented   = ft_volumesegment(cfg, mri);
    
    % Sanity check for adequate segmentation
    if flag_check == 1
        seg_i = ft_datatype_segmentation(mri_segmented, ...
            'segmentationstyle', 'indexed');
        
        cfg              = [];
        cfg.funparameter = 'tissue';
        cfg.funcolormap  = ...
            lines(numel(seg_i.tissuelabel)); % colors according to tissuelabel
        cfg.location     = 'center';
        ft_sourceplot(cfg, seg_i);
    end
    
    fprintf("\n\t ... Done! \t\t\n")
    save(filename_segmented, 'mri_segmented', '-v7.3');
else
    load(filename_segmented)
    fprintf("\n\tStep 1: Already completed for subj: %s\n", upper(code_parti))
end

%% Create mesh for important compartments
fprintf("\n\tStep 2: Create mesh for compartments (subj: %s)\n", upper(code_parti))
filename_mesh = fullfile(wdir, 'sourcemodel', ...
    sprintf('mesh_%s.mat', code_parti));
if ~exist(fullfile(wdir, "sourcemodel"), 'dir'); mkdir(fullfile(wdir, "sourcemodel")); end

if ~exist(filename_mesh, 'file')
    fprintf("\n\t ... Processing ... \t\t\n")
    
    compartments = {'brain', 'skull', 'scalp'};
    vertices     = [3000, 2000, 1000];
    mesh_brain = cell(3,1);
    cfg             = [];
    for k = 1:numel(compartments)
        cfg.method      = 'projectmesh';
        cfg.tissue      = compartments{k};
        cfg.numvertices = vertices(k);
        mesh_brain{k} = ft_prepare_mesh(cfg, mri_segmented);
        mesh_brain{k} = ft_convert_units(mesh_brain{k}, 'm'); % Use SI Units % Use SI Units
    end
    mesh_eeg = [mesh_brain{1,1} mesh_brain{2} mesh_brain{3}]; % concatenate results
    fprintf("\n\t ... Done! \t\t\n")
    save(filename_mesh, 'mesh_eeg', '-v7.3'); 
    clear mesh_brain compartments vertices
else
    load(filename_mesh)
    fprintf("\n\tStep 2: Already completed (subj: %s)\n", upper(code_parti))
end

if flag_check == 1
    figure % plot the created mesh for all tissue types of interest
    ft_plot_mesh(mesh_eeg(1), 'edgecolor', 'none', 'facecolor', 'r') % brain
    ft_plot_mesh(mesh_eeg(2), 'edgecolor', 'none', 'facecolor', 'g') % skull
    ft_plot_mesh(mesh_eeg(3), 'edgecolor', 'none', 'facecolor', 'b') % scalp
    alpha 0.3; view(132, 14)
end

%% Create headmodel (returns error on a Windows machine)
fprintf("\n\tStep 3: Create headmodel for subj: %s.\n", upper(code_parti))
filename_headmodel = fullfile(wdir, 'sourcemodel', ...
    sprintf('headmodel_%s.mat', code_parti));

if ispc
    fprintf("\nHeadmodel cannot be constructed with dipoli option on a window machine. Please reconsider other method or machine!")
    return
end

if ~exist(filename_headmodel, 'file')
    fprintf("\n\t ... Processing ... \t\t\n")
    
    cfg              = [];
    cfg.method       = 'dipoli'; % will not work with Windows
    cfg.conductivity = [1 1/80 1] * (1/3); % S/m
    headmodel_dipoli = ft_prepare_headmodel(cfg, mesh_eeg);
    headmodel_dipoli = ft_convert_units(headmodel_dipoli, 'm'); % Use SI Units
    fprintf("\n\t ... Done! \t\t\n")
    
    save(filename_headmodel, 'headmodel_dipoli', '-v7.3');
else
    load(filename_headmodel) %#ok<*LOAD>
    fprintf("\n\tStep 3: Already completed (subj: %s)\n", upper(code_parti))
end

if flag_check == 1
    figure % plot the resulting headmodel for all tissue types of interest
    ft_plot_headmodel(headmodel_dipoli, 'facealpha', 0.5)
    view(90, 0)
end

%% Create sourcemodel for every subject
fprintf("\n\tStep 4: Create sourcemodel (subj: %s)\n", upper(code_parti))
filename_sourcemodel = fullfile(wdir, 'sourcemodel', ...
    sprintf('sourcemodel_%s.mat', code_parti));

if ~exist(filename_sourcemodel, 'file')
    fprintf("\n\t ... Processing ... \t\t\n")
    cfg             = [];
    cfg.headmodel   = headmodel_dipoli; % used to estimate extent of grid
    cfg.resolution  = 0.01; % a source per 0.008 m -> .8cm
    cfg.inwardshift = 0.001; % moving sources 2 mm inwards from the skull
    sourcemodel = ft_prepare_sourcemodel(cfg);
    fprintf("\n\t ... Done! \t\t\n")
    save(filename_sourcemodel, 'sourcemodel', '-v7.3');
else
    load(filename_sourcemodel) %#ok<*LOAD>
    fprintf("\n\tStep 4: already completed (subj: %s)\n", upper(code_parti))
end

if flag_check == 1
    inside = sourcemodel; outside = sourcemodel;
    inside.pos = sourcemodel.pos(sourcemodel.inside, :);
    outside.pos = sourcemodel.pos(~sourcemodel.inside, :);
    
    figure; hold on
    ft_plot_mesh(inside, 'vertexsize', 20, 'vertexcolor', 'red');
    ft_plot_mesh(outside, 'vertexsize', 20)
    ft_plot_headmodel(headmodel_dipoli, 'facealpha', 0.1)
    view(125, 10)
end

%% Re-align electrodes and create a leadfield

% Prepare layout according to the 1005 template from fieldtrip
elec_default = ft_read_sens('standard_1005.elc');
elec_default = ft_convert_units(elec_default, 'm');

if strcmp(headmodel_dipoli.type, 'bemcp')
    scalp_index = 3;
elseif strcmp(headmodel_dipoli.type, 'dipoli') % quite a bit funny here!
    scalp_index = 1;
end

cfg = [];
cfg.method = 'project'; % onto scalp surface
cfg.headshape = headmodel_dipoli.bnd(scalp_index); % scalp surface
elec_realigned = ft_electroderealign(cfg, elec_default);
elec_realigned = ft_convert_units(elec_realigned, 'm');

filename_realignedEEG = fullfile(wdir, 'sourcemodel', ...
    sprintf('elecRealigned_%s.mat', code_parti));
save(filename_realignedEEG , 'elec_realigned', '-v7.3');

if flag_check == 1
    figure; hold on;
    ft_plot_sens(elec_realigned, 'elecsize', 40);
    ft_plot_headmodel(headmodel_dipoli, 'facealpha', 0.5);
    view(90, 0)
end

% Create leadfield
filename_leadfield = fullfile(wdir, 'leadfield', ...
    sprintf('leadfield_%s.mat', code_parti));
if ~exist(fullfile(wdir, "leadfield"), 'dir'); mkdir(fullfile(wdir, "leadfield")); end

cfg = [];
cfg.sourcemodel = sourcemodel;
cfg.headmodel   = headmodel_dipoli;
cfg.elec        = elec_realigned;
leadfield = ft_prepare_leadfield(cfg);

save(filename_leadfield , 'leadfield', '-v7.3');

if flag_check == 1
    sanity_check_leadfield(elec_realigned, sourcemodel_and_leadfield, ...
        headmodel_dipoli)
end

end

function check_for_template(wdir)

%   This function checks for whether the templates were already processed;
%   in case segmentedMRI is missing, processing is started

%   ## Version 1.1; changed the unnecessary tpm option

%   Copyright (C) October 2021, modified December 2021
%   D. Pedrosa
%   University Hospital of Gießen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

files2checkfor = {fullfile(wdir, 'mri_preprocessed', ...
    'segmentedMRI_template.mat')};
allFiles = dir(fullfile(wdir, 'mri_preprocessed'));
logical_files = ismember(files2checkfor,{allFiles.name});
if all(logical_files) % continue as all files already present
    return
else
    filename_mri_template = fullfile(wdir, 'templateMRI', ...
        'MRItrem_template0.nii.gz');
    mri = ft_read_mri(filename_mri_template);
    mri.coordsys = 'lps';
    
    %  Segment (template) MRI to different compartments
    filename_segmented = fullfile(wdir, 'mri_preprocessed', ...
        'segmented_MRItemplate.mat');
    
    cfg.output      = {'gray','white','csf','skull','scalp'};
    mri_segmented   = ft_volumesegment(cfg, mri);
    save(filename_segmented, 'mri_segmented', '-v7.3');
     
    % Create mesh
    filename_mesh = fullfile(wdir, 'sourcemodel', ...
        sprintf('mesh_template.mat'));
    if ~exist(fullfile(wdir, "sourcemodel"), 'dir'); mkdir(fullfile(wdir, "sourcemodel")); end
    
    compartments = {'brain', 'skull', 'scalp'};
    vertices     = [3000, 2000, 1000];
    mesh_brain = cell(3,1);
    cfg             = [];
    for k = 1:numel(compartments)
        cfg.method      = 'projectmesh';
        cfg.tissue      = compartments{k};
        cfg.numvertices = vertices(k);
        mesh_brain{k} = ft_prepare_mesh(cfg, mri_segmented);
        mesh_brain{k} = ft_convert_units(mesh_brain{k}, 'm'); % Use SI Units % Use SI Units
    end
    mesh_eeg = [mesh_brain{1,1} mesh_brain{2} mesh_brain{3}]; % concatenate results
    save(filename_mesh, 'mesh_eeg', '-v7.3');
    
    % Create headmodel
    filename_headmodel = fullfile(wdir, 'sourcemodel', ...
        'headmodel_template.mat');
    
    if ispc
        fprintf("\nHeadmodel cannot be constructed with dipoli option on a window machine. Please reconsider other method or machine!")
        return
    end
    
    cfg              = [];
    cfg.method       = 'dipoli'; % will not work with Windows
    cfg.conductivity = [1 1/80 1] * (1/3); % S/m
    headmodel_dipoli = ft_prepare_headmodel(cfg, mesh_eeg);
    headmodel_dipoli = ft_convert_units(headmodel_dipoli, 'm'); % Use SI Units
    save(filename_headmodel, 'headmodel_dipoli', '-v7.3');
    
    % Create sourcemodel
    filename_sourcemodel = fullfile(wdir, 'sourcemodel', ...
        'sourcemodel_template.mat');
    if ~exist(fullfile(wdir, "sourcemodel"), 'dir'); mkdir(fullfile(wdir, "sourcemodel")); end
    
    cfg             = [];
    cfg.headmodel   = headmodel_dipoli; % used to estimate extent of grid
    cfg.resolution  = 0.01; % a source per 0.008 m -> .8cm
    cfg.inwardshift = 0.001; % moving sources 2 mm inwards from the skull
    sourcemodel = ft_prepare_sourcemodel(cfg);
    fprintf("\n\t ... Done! \t\t\n")
    save(filename_sourcemodel, 'sourcemodel', '-v7.3');
    
    % Realign electrodes to headmodel
    filename_realignedEEG = fullfile(wdir, 'sourcemodel', ...
        'elecRealigned_template.mat');
    
    elec_default = ft_read_sens('standard_1005.elc');
    elec_default = ft_convert_units(elec_default, 'm');
    
    if strcmp(headmodel_dipoli.type, 'bemcp')
        scalp_index = 3;
    elseif strcmp(headmodel_dipoli.type, 'dipoli') % quite a bit funny here!
        scalp_index = 1;
    end
    
    cfg = [];
    cfg.method = 'project'; % onto scalp surface
    cfg.headshape = headmodel_dipoli.bnd(scalp_index); % scalp surface
    elec_realigned = ft_electroderealign(cfg, elec_default);
    elec_realigned = ft_convert_units(elec_realigned, 'm');
    
    save(filename_realignedEEG , 'elec_realigned', '-v7.3');
    
    % Create leadfield
    filename_leadfield = fullfile(wdir, 'leadfield', ...
        'leadfield_template.mat');
    if ~exist(fullfile(wdir, "leadfield"), 'dir'); mkdir(fullfile(wdir, "leadfield")); end
    
    cfg = [];
    cfg.sourcemodel = sourcemodel;
    cfg.headmodel   = headmodel_dipoli;
    cfg.elec        = elec_realigned;
    
    leadfield = ft_prepare_leadfield(cfg);
    save(filename_leadfield , 'leadfield', '-v7.3'); 
    
end
end

function sanity_check_leadfield(elec_realigned, ...
    sourcemodel_and_leadfield, headmodel)

%   This function displays leadfields along with the headmodel and
%   the electrodes to figure out problems with the estimation; adapted from:
%   https://www.fieldtriptoolbox.org/workshop/oslo2019/forward_modeling/#head-models-component-1

%   ## Version 1.0

%   Copyright (C) November 2021
%   D. Pedrosa
%   University Hospital of Gießen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.


figure('units', 'normalized', 'outerposition', [0 0 0.5 0.5])
source_index = 251; %% random source
sensory_dipole_current = 100e-9; % Am (realistic)

n_sensors = length(elec_realigned.label);

inside_sources = find(sourcemodel_and_leadfield.inside);
inside_index = inside_sources(source_index);
lead = sourcemodel_and_leadfield.leadfield{inside_index};
xs = zeros(1, n_sensors);
ys = zeros(1, n_sensors);
zs = zeros(1, n_sensors);
voltages = zeros(1, n_sensors);
titles = {'Lead field (x)' 'Lead field (y)' 'Lead field (z)'};

% get the xyz and norm

for sensor_index = 1:n_sensors
    this_x = lead(sensor_index, 1);
    this_y = lead(sensor_index, 2);
    this_z = lead(sensor_index, 3);
    this_norm = norm(lead(sensor_index, :));
    xs(sensor_index) = this_x * sensory_dipole_current;
    ys(sensor_index) = this_y * sensory_dipole_current;
    zs(sensor_index) = this_z * sensory_dipole_current;
    voltages(sensor_index) = this_norm * sensory_dipole_current;
end

% plot xyz
axes = {xs ys zs};

for axis_index = 1:3
    this_axis = axes{axis_index};
    subplot(1, 3, axis_index)
    hold on
    ft_plot_topo3d(elec_realigned.chanpos, this_axis, 'facealpha', 0.8)
    if strcmp(headmodel.type, 'dipoli')
        caxis([-10e-6, 10e-6])
    end
    c = colorbar('location', 'southoutside');
    c.Label.String = 'Lead field (V)';
    axis tight
    ft_plot_mesh(mesh_eeg, 'facealpha', 0.10);
    ft_plot_sens(elec_realigned, 'elecsize', 20);
    title(titles{axis_index})
    plot3(sourcemodel_and_leadfield.pos(inside_index, 1), ...
        sourcemodel_and_leadfield.pos(inside_index, 2), ...
        sourcemodel_and_leadfield.pos(inside_index, 3), 'bo', ...
        'markersize', 20, 'markerfacecolor', 'r')
end

% plot norm
figure('units', 'normalized', 'outerposition', [0 0 0.5 0.85])
hold on
ft_plot_topo3d(elec_realigned.chanpos, voltages, 'facealpha', 0.8)
if strcmp(headmodel.type, 'dipoli')
    caxis([0, 10e-6])
end
c = colorbar('location', 'eastoutside');
c.Label.String = 'Lead field (V)';
axis tight
ft_plot_mesh(mesh_eeg, 'facealpha', 0.10);
ft_plot_sens(elec_realigned, 'elecsize', 20);
title('Leadfield magnitude')
plot3(sourcemodel_and_leadfield.pos(inside_index, 1), ...
    sourcemodel_and_leadfield.pos(inside_index, 2), ...
    sourcemodel_and_leadfield.pos(inside_index, 3), 'bo', ...
    'markersize', 20, 'markerfacecolor', 'r')
view(-90, 0)
keyboard
end