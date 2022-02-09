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
            mri.coordsys = 'mni';
            
            % Convert to MNI space, as this makes later steps easier
            cfg = [];
            cfg.nonlinear   = 'no';
            mri             = ft_volumenormalise(cfg, mri);
            
            cfg = [];
            cfg.resolution = 1;
            mri = ft_volumereslice(cfg, mri);
            
        catch
            fprintf("\nMRI for subj: %s not existing, will use template instead\n", upper(code_parti))
            check_for_template(wdir)
            return % in case templates are necessary, this should be performed once, so that a return command is set
        end
end


%%  Segment MRI to different compartments
fprintf("\n\tStep 1: Segmenting MRI with different compartments for subj: %s\n", upper(code_parti))
filename_segmented = fullfile(wdir, 'mri_preprocessed', ...
    sprintf('segmentedMRI_%s.mat', code_parti));

filename_segmented2 = fullfile(wdir, 'mri_preprocessed', ...
    sprintf('segmentedMRI_detailed_%s.mat', code_parti));

if ~exist(filename_segmented, 'file')
    fprintf("\n\t ... Processing ... \t\t\n")
    
    cfg = [];
    cfg.spmmethod           = 'new';
    mri_segmented_detailed  = ft_volumesegment(cfg, mri);
    
    cfg.spmmethod           ='old';
    cfg.output              = {'skull', 'scalp', 'brain'};
    mri_segmented           = ft_volumesegment(cfg, mri);
    
    manual_changes = 0;
    if manual_changes == 1
        mri_segmented.skull = imerode(mri_segmented.skull, strel_bol(1));
        mri_segmented.scalp = imdilate(mri_segmented.scalp, strel_bol(2));

        mri_segmented.scalp = volumethreshold(mri_segmented.scalp, ...
            .5, 'scalpmask');

        a1 = volumefillholes(mri_segmented.scalp, 1);
        a2 = volumefillholes(mri_segmented.scalp, 2);
        a3 = volumefillholes(mri_segmented.scalp, 3);
        mri_segmented.scalp = a1 | a2 | a3;

        % threshold again to remove little parts outside of head
        mri_segmented.scalp = volumethreshold(mri_segmented.scalp, ...
            .5, 'scalpmask');

        scalpMask = mri_segmented.scalp;
        for i = 1:size(scalpMask,3)
            scalpMask(:,:,i) = imfill(squeeze(scalpMask(:,:,i)),8,'holes');
        end

        for i = 1:size(scalpMask,2)
            scalpMask(:,i,:) = imfill(squeeze(scalpMask(:,i,:)),8,'holes');
        end

        for i = 1:size(scalpMask,1)
            scalpMask(i,:,:) = imfill(squeeze(scalpMask(i,:,:)),8,'holes');
        end
        mri_segmented.scalp = scalpMask;
        mri_segmented.scalp = volumethreshold(mri_segmented.scalp, ...
            .5, 'scalpmask');
        
       
    end    
        
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
    save(filename_segmented2, 'mri_segmented_detailed', '-v7.3');
else
    load(filename_segmented); load(filename_segmented2)
    fprintf("\n\t\tStep 1: Already completed for subj: %s\n", upper(code_parti))
end

%% Create mesh for important compartments
fprintf("\n\tStep 2: Create mesh for compartments (subj: %s)\n", upper(code_parti))
filename_mesh = fullfile(wdir, 'sourcemodel', ...
    sprintf('mesh_%s.mat', code_parti));
if ~exist(fullfile(wdir, "sourcemodel"), 'dir'); mkdir(fullfile(wdir, "sourcemodel")); end

if ~exist(filename_mesh, 'file')
    fprintf("\n\t ... Processing ... \t\t\n")
    addpath(genpath('/opt/fieldtrip/external/iso2mesh'))
    
    repair_mesh = 1;
    compartments = {'brain', 'skull', 'scalp'};
    if repair_mesh == 1
        vertices = [1 1 1].*10000;
    else
        vertices = [3000, 2000, 1000];
    end
    
    mesh_brain = cell(length(compartments),1);
    cfg             = [];
    for k = 1:numel(compartments)
        cfg.method      = 'projectmesh';
        cfg.tissue      = compartments{k};
        cfg.numvertices = vertices(k);
        mesh_brain{k} = ft_prepare_mesh(cfg, mri_segmented);
    end
    
    mesh_eeg = [mesh_brain{1} mesh_brain{2} mesh_brain{3}]; % concatenate results
    
    if repair_mesh == 1 % according to https://github.com/meeg-cfin/nemolab/blob/master/basics/nemo_mriproc.m
        targetsize = [3000, 2000, 1000]; % downsample to the desired number of vertices
        bnd = mesh_eeg;
        for layer = 1:numel(bnd)
            for ii = 1:length(bnd)
                [bnd(ii).pos, bnd(ii).tri] = meshresample(bnd(ii).pos, ...
                    bnd(ii).tri, targetsize(ii)/size(bnd(ii).pos,1));
                [bnd(ii).pos, bnd(ii).tri] = ...
                    meshcheckrepair(bnd(ii).pos, bnd(ii).tri, 'dup');
                [bnd(ii).pos, bnd(ii).tri] = ...
                    meshcheckrepair(bnd(ii).pos, bnd(ii).tri, 'isolated');
                [bnd(ii).pos, bnd(ii).tri] = ...
                    meshcheckrepair(bnd(ii).pos, bnd(ii).tri, 'deep');
                [bnd(ii).pos, bnd(ii).tri] = ...
                    meshcheckrepair(bnd(ii).pos, bnd(ii).tri, 'meshfix');
            end
        end
        mesh_eeg = bnd;
    end
    
    fprintf("\n\t ... Done! \t\t\n")
    save(filename_mesh, 'mesh_eeg', '-v7.3');
    
    clear mesh_brain compartments vertices
else
    load(filename_mesh)
    fprintf("\n\t\tStep 2: Already completed (subj: %s)\n", upper(code_parti))
end

if flag_check == 1
    figure % plot the created mesh for all tissue types of interest
    ft_plot_mesh(mesh_eeg(1), 'edgecolor', 'black', 'facecolor', 'r') % brain
    ft_plot_mesh(mesh_eeg(2), 'edgecolor', 'none', 'facecolor', 'g') % skull
    ft_plot_mesh(mesh_eeg(3), 'edgecolor', 'none', 'facecolor', 'b') % scalp
    alpha 0.3; view(132, 14)
end

%% Create headmodel (returns error on a Windows machine)
fprintf("\n\tStep 3: Create headmodel for subj: %s.\n", upper(code_parti))
filename_headmodel = fullfile(wdir, 'sourcemodel', ...
    sprintf('headmodel_%s.mat', code_parti));

if ~exist(filename_headmodel, 'file')
    fprintf("\n\t ... Processing ... \t\t\n")
    if isunix
        setenv('PATH', ['/opt/openMEEG/bin:' getenv('PATH')]);
        setenv('LD_LIBRARY_PATH', ['/opt/openMEEG/lib:' getenv('LD_LIBRARY_PATH')]);
    end
    
    cfg = [];
    cfg.method = 'openmeeg';
    try
        headmodel = ft_prepare_headmodel(cfg,mesh_eeg);
    catch
        fprintf("\nManual refinement necessary for the mesh/headmodel! PLease double-check")
        keyboard
    end
    headmodel = ft_convert_units(headmodel, 'mm'); % Use SI Units
    fprintf("\n\t ... Done! \t\t\n")
    
    save(filename_headmodel, 'headmodel', '-v7.3');
else
    load(filename_headmodel) %#ok<*LOAD>
    fprintf("\n\t\tStep 3: Already completed (subj: %s)\n", upper(code_parti))
end

if flag_check == 1
    figure % plot the resulting headmodel for all tissue types of interest
    ft_plot_headmodel(headmodel, 'facealpha', 0.5)
    view(90, 0)
end

%% Load reference (template) sourcemodel (template_grid, cf.
% http://old.fieldtriptoolbox.org/tutorial/sourcemodel) and
% {create_template_grid.m}

filename_grid = fullfile(wdir, 'templateMRI', ...
    'template_grid.mat'); load(filename_grid);

if flag_check == 1
    figure; hold on
    headmodel = ft_convert_units(headmodel, 'mm');
    ft_plot_headmodel(headmodel, 'facecolor', 'cortex', ...
        'edgecolor', 'none'); alpha 0.5; camlight;
    ft_plot_mesh(template_grid.pos(template_grid.inside,:)); %#ok<NODEF>
end


%% Create sourcemodel for every subject
fprintf("\n\tStep 4: Create sourcemodel (subj: %s)\n", upper(code_parti))
filename_sourcemodel = fullfile(wdir, 'sourcemodel', ...
    sprintf('sourcemodel_%s.mat', code_parti));
template_grid = ft_convert_units(template_grid, 'mm');

if ~exist(filename_sourcemodel, 'file')
    fprintf("\n\t ... Processing ... \t\t\n")
    cfg             = [];
    cfg.template  = template_grid;
    cfg.nonlinear = 'yes';
    cfg.mri       = mri;
    cfg.unit      ='mm';
    
    cfg.headmodel   = headmodel; % used to estimate extent of grid
    cfg.inwardshift = 3; % moving sources 3 mm inwards from the skull
    cfg.method      = 'basedonmni';
    %cfg.template    = template_grid;
    
    sourcemodel = ft_prepare_sourcemodel(cfg);
    fprintf("\n\t ... Done! \t\t\n")
    save(filename_sourcemodel, 'sourcemodel', '-v7.3');
else
    load(filename_sourcemodel) %#ok<*LOAD>
    fprintf("\n\t\tStep 4: already completed (subj: %s)\n", upper(code_parti))
end

if flag_check == 1
    inside = sourcemodel; outside = sourcemodel;
    inside.pos = sourcemodel.pos(sourcemodel.inside, :);
    outside.pos = sourcemodel.pos(~sourcemodel.inside, :);
    
    figure; hold on;
    headmodel  = ft_convert_units(headmodel, 'mm');
    ft_plot_mesh(inside, 'vertexsize', 20, 'vertexcolor', 'red');
    ft_plot_mesh(outside, 'vertexsize', 20)
    ft_plot_mesh(template_grid.pos(template_grid.inside,:), 'vertexsize', 20, 'vertexcolor', 'blue')
    ft_plot_headmodel(headmodel, 'facealpha', 0.1)
    view(125, 10)
end

%% Re-align electrodes and create a leadfield

% Prepare layout according to the 1005 template from fieldtrip
elec_default = ft_read_sens('standard_1005.elc');
elec_default = ft_convert_units(elec_default, 'mm');

if strcmp(headmodel.type, 'bemcp')
    scalp_index = 3;
elseif strcmp(headmodel.type, 'dipoli') % quite a bit funny here!
    scalp_index = 1;
else
    scalp_index = 1;
end

cfg = [];
cfg.method = 'project'; % onto scalp surface
cfg.headshape = headmodel.bnd(scalp_index); % scalp surface
elec_realigned = ft_electroderealign(cfg, elec_default);
elec_realigned = ft_convert_units(elec_realigned, 'mm');

filename_realignedEEG = fullfile(wdir, 'sourcemodel', ...
    sprintf('elecRealigned_%s.mat', code_parti));
save(filename_realignedEEG , 'elec_realigned', '-v7.3');

if flag_check == 1
    figure; hold on;
    ft_plot_sens(elec_realigned, 'elecsize', 40);
    ft_plot_headmodel(headmodel, 'facealpha', 0.5);
    view(90, 0)
end

% Create leadfield
filename_leadfield = fullfile(wdir, 'leadfield', ...
    sprintf('leadfield_%s.mat', code_parti));
if ~exist(fullfile(wdir, "leadfield"), 'dir'); mkdir(fullfile(wdir, "leadfield")); end

if ~exist(filename_leadfield, 'file')
    cfg = [];
    cfg.sourcemodel.pos = sourcemodel.pos;
    cfg.headmodel   = headmodel;
    cfg.elec        = elec_realigned;
    leadfield = ft_prepare_leadfield(cfg);
    
    save(filename_leadfield , 'leadfield', '-v7.3');
else
    fprintf("\n\t\tStep 5: already completed (subj: %s)\n", upper(code_parti))
end

if flag_check == 1
    sanity_check_leadfield(elec_realigned, leadfield, ...
        headmodel, mesh_eeg)
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

files2checkfor = {fullfile(wdir, 'templateMRI', ...
    'segmentedMRI_template.mat')};
allFiles = dir(fullfile(wdir, 'templateMRI'));
logical_files = ismember(files2checkfor,{allFiles.name});
if all(logical_files) % continue as all files already present
    return
else
    filename_mri_template = fullfile(wdir, 'templateMRI', ...
        'MRItrem_template0.nii.gz');
    mri = ft_read_mri(filename_mri_template);
    mri.coordsys = 'lps';
    
    % Convert to MNI space, as this makes later steps easier
    cfg = [];
    cfg.nonlinear   = 'yes';
    mri             = ft_volumenormalise(cfg, mri);
    
    %  Segment (template) MRI to different compartments
    filename_segmented = fullfile(wdir, 'templateMRI', ...
        'segmented_MRItemplate.mat');
    
    cfg = [];
    cfg.spmmethod='new';
    mri_segmented = ft_volumesegment(cfg, mri);
    save(filename_segmented, 'mri_segmented', '-v7.3');
    
    % Create mesh
    filename_mesh = fullfile(wdir, 'templateMRI', ...
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
        mesh_brain{k} = ft_convert_units(mesh_brain{k}, 'mm'); % Use SI Units
    end
    mesh_eeg = [mesh_brain{1,1} mesh_brain{2} mesh_brain{3}]; % concatenate results
    save(filename_mesh, 'mesh_eeg', '-v7.3');
    
    
    % Create sourcemodel
    filename_sourcemodel = fullfile(wdir, 'templateMRI', ...
        'sourcemodel_template.mat');
    if ~exist(fullfile(wdir, "sourcemodel"), 'dir'); mkdir(fullfile(wdir, "sourcemodel")); end
    
    % Construct dipole grid in template brain coordinates
    cfg             = [];
    cfg.grid.tight  = 'yes';
    cfg.headmodel   = headmodel; % used to estimate extent of grid
    cfg.resolution  = 8; % a source per 0.008 m -> .8cm
    cfg.inwardshift = -1.5; % moving sources 1.5 mm inwards from the skull
    sourcemodel = ft_prepare_sourcemodel(cfg);
    fprintf("\n\t ... Done! \t\t\n")
    save(filename_sourcemodel, 'sourcemodel', '-v7.3');
    
    % Realign electrodes to headmodel
    filename_realignedEEG = fullfile(wdir, 'templateMRI', ...
        'elecRealigned_template.mat');
    
    elec_default = ft_read_sens('standard_1005.elc');
    elec_default = ft_convert_units(elec_default, 'mm');
    
    if strcmp(headmodel.type, 'bemcp')
        scalp_index = 3;
    elseif strcmp(headmodel.type, 'dipoli') % quite a bit funny here!
        scalp_index = 1;
    else
        scalp_index = 3;
    end
    
    cfg = [];
    cfg.method = 'project'; % onto scalp surface
    cfg.headshape = headmodel.bnd(scalp_index); % scalp surface
    elec_realigned = ft_electroderealign(cfg, elec_default);
    elec_realigned = ft_convert_units(elec_realigned, 'mm');
    
    save(filename_realignedEEG , 'elec_realigned', '-v7.3');
    
    % Create leadfield
    filename_leadfield = fullfile(wdir, 'templateMRI', ...
        'leadfield_template.mat');
    if ~exist(fullfile(wdir, "leadfield"), 'dir'); mkdir(fullfile(wdir, "leadfield")); end
    
    cfg = [];
    cfg.sourcemodel = sourcemodel;
    cfg.headmodel   = headmodel;
    cfg.elec        = elec_realigned;
    
    leadfield = ft_prepare_leadfield(cfg);
    save(filename_leadfield , 'leadfield', '-v7.3');
    
end
end

function sanity_check_leadfield(elec_realigned, ...
    leadfield, headmodel, mesh_eeg)

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
source_index = 1270; %% random source
sensory_dipole_current = 100e-9; % Am (realistic)

n_sensors = length(elec_realigned.label);

inside_sources = find(leadfield.inside);
inside_index = inside_sources(source_index);
lead = leadfield.leadfield{inside_index};
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
    elseif strcmp(headmodel.type, 'openmeeg')
        caxis([0, max(voltages)])
    end
    c = colorbar('location', 'southoutside');
    c.Label.String = 'Lead field (V)';
    axis tight
    ft_plot_mesh(mesh_eeg, 'facealpha', 0.10);
    ft_plot_sens(elec_realigned, 'elecsize', 20);
    title(titles{axis_index})
    plot3(leadfield.pos(inside_index, 1), ...
        leadfield.pos(inside_index, 2), ...
        leadfield.pos(inside_index, 3), 'bo', ...
        'markersize', 20, 'markerfacecolor', 'r')
end

% plot norm
figure('units', 'normalized', 'outerposition', [0 0 0.5 0.85])
hold on
ft_plot_topo3d(elec_realigned.chanpos, voltages, 'facealpha', 0.8)
if strcmp(headmodel.type, 'dipoli')
    caxis([0, 10e-6])
elseif strcmp(headmodel.type, 'openmeeg')
    caxis([0, max(voltages)])
end
c = colorbar('location', 'eastoutside');
c.Label.String = 'Lead field (V)';
axis tight
ft_plot_mesh(mesh_eeg, 'facealpha', 0.10);
ft_plot_sens(elec_realigned, 'elecsize', 20);
title('Leadfield magnitude')
plot3(leadfield.pos(inside_index, 1), ...
    leadfield.pos(inside_index, 2), ...
    leadfield.pos(inside_index, 3), 'bo', ...
    'markersize', 20, 'markerfacecolor', 'r')
view(-90, 0)
keyboard
end