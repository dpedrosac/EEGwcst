function create_template_grid(t1_filename, wdir, flag_check)

%   This function creates a template based grid which may be used for
%   the creation of a grid to use later

%   ## Version 1.1, changed to SST instead of MNI template

%   Copyright (C) September 2021
%   D. Pedrosa, University Hospital of Gießen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

if nargin < 3
    flag_check = 0;
end

filename_segmented = fullfile(wdir, 'templateMRI', ...
    'segmentedMRI_template.mat');
filename_segmented2 = fullfile(wdir, 'templateMRI', ...
    'segmentedMRI_template_detailed.mat');


%% Create templare according to http://old.fieldtriptoolbox.org/tutorial/sourcemodel
template = ft_read_mri(t1_filename);
template.coordsys = 'mni';

% Convert to MNI space, as this makes later steps easier
cfg = [];
cfg.nonlinear   = 'no';
template        = ft_volumenormalise(cfg, template);

cfg = [];
cfg.resolution = 1;
template = ft_volumereslice(cfg, template);

save(fullfile(wdir, 'templateMRI', 'template_affineSPM.mat'), ...
    'template', '-v7.3')

% Segment template brain and construct volume conduction model (i.e. head model):
cfg          = [];
cfg.spmmethod='new';
mri_segmented_detailed = ft_volumesegment(cfg, template);

cfg.spmmethod='old';
cfg.output      = {'skull', 'scalp', 'brain'};
mri_segmented = ft_volumesegment(cfg, template);


if flag_check == 1
    seg_i = ft_datatype_segmentation(mri_segmented, ...
        'segmentationstyle', 'indexed');
    cfg              = [];
    cfg.funparameter = 'tissue';
    cfg.funcolormap  = ...
        lines(numel(seg_i.tissuelabel)); % colors according to tissuelabel
    cfg.location     = 'center';
    ft_sourceplot(cfg, seg_i);
    
    cfg = [];
    template_plot = template;
    template_plot.test = mri_segmented.scalp;
    cfg.funparameter  = 'test';
    ft_sourceplot(cfg, template_plot);
    
end

save(filename_segmented, 'mri_segmented', '-v7.3');
save(filename_segmented2, 'mri_segmented_detailed', '-v7.3');

% Create mesh for template brain:
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

if flag_check == 1
    figure % plot the created mesh for all tissue types of interest
    ft_plot_mesh(mesh_eeg(1), 'edgecolor', 'black', 'facecolor', 'r') % brain
    ft_plot_mesh(mesh_eeg(2), 'edgecolor', 'none', 'facecolor', 'g') % skull
    ft_plot_mesh(mesh_eeg(3), 'edgecolor', 'none', 'facecolor', 'b') % scalp
    alpha 0.3; view(132, 14)
end

save(fullfile(wdir, 'templateMRI', 'mesh_template.mat'), ...
    'mesh_eeg', '-v7.3')

% Create the headmodel usinng the openMEEG toolbox

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

save(fullfile(wdir, 'templateMRI', 'headmodel_template.mat'), ...
    'headmodel', '-v7.3')

% Old method, deprecated
% cfg          = [];
% cfg.method   = 'dipoli'; %'singleshell';
% template_headmodel = ft_prepare_headmodel(cfg, mri_segmented);
% template_headmodel = ft_convert_units(template_headmodel, 'mm'); % Convert the vol to mm, because the CTF convenction is to express everything in cm.


% Construct dipole grid in template brain coordinates
cfg = [];
cfg.grid.resolution = 8;
cfg.grid.tight      = 'yes';
cfg.inwardshift     = -1.5;                                                 % negative inwardshift means an outward shift of the brain surface for inside/outside detection
cfg.headmodel       = headmodel;
template_grid       = ft_prepare_sourcemodel(cfg);
template_grid       = ft_convert_units(template_grid, 'mm'); % Convert the vol to cm, because the CTF convenction is to express everything in cm.


%% Create figure with template head model and dipole grid
if flag_check == 1
    figure; hold on
    ft_plot_headmodel(headmodel, 'facecolor', 'cortex', ...
        'edgecolor', 'none');alpha 0.5; camlight;
    ft_plot_mesh(template_grid.pos(template_grid.inside,:));
    keyboard
end

save(fullfile(wdir, 'templateMRI', 'template_grid.mat'), ...
    'template_grid', '-v7.3')


end
