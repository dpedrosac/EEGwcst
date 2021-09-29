function create_template_grid(t1_filename, wdir, flag_check)

%   This function creates a template based grid which may be used for
%   the creation of a grid to use later

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

%% Create templare according to http://old.fieldtriptoolbox.org/tutorial/sourcemodel
template = ft_read_mri(t1_filename);
template.coordsys = 'acpc'; % so that FieldTrip knows how to interpret the coordinate system

% Segment template brain and construct volume conduction model (i.e. head model):
cfg          = [];
template_seg = ft_volumesegment(cfg, template);

cfg          = [];
cfg.method   = 'dipoli'; %'singleshell';
template_headmodel = ft_prepare_headmodel(cfg, template_seg);
template_headmodel = ft_convert_units(template_headmodel, 'cm'); % Convert the vol to cm, because the CTF convenction is to express everything in cm.

% Construct dipole grid in template brain coordinates
cfg = [];
cfg.grid.resolution = 1;
cfg.grid.tight      = 'yes';
cfg.inwardshift     = -1.5;                                                 % negative inwardshift means an outward shift of the brain surface for inside/outside detection
cfg.headmodel       = template_headmodel;
template_grid       = ft_prepare_sourcemodel(cfg);

%% Create figure with the template head model and dipole grid
if flag_check == 1
    figure
    hold on
    ft_plot_headmodel(template_headmodel, 'facecolor', 'cortex', ...
        'edgecolor', 'none');alpha 0.5; camlight;
    ft_plot_mesh(template_grid.pos(template_grid.inside,:));
    keyboard
end

save(fullfile(wdir, 'templateMRI', 'template_grid.mat'), ...
    'template_grid', '-v7.3')

end
