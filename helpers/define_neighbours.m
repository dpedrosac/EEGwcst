function neighbours = ...
    define_neighbours(data_noica, loaddir, dist, flag_check)

%   This function runs the routine to obtain the neighbours later needed to 
%   interpolate the 'bad channels'

%   Copyright (C) January 2018
%   D. Pedrosa, University Hospital of Gieﬂen and Marburg

%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.


data_noica.elec = ...                           % this block assigns the electrode positions according to BRAINVISION to the data
    ft_read_sens([loaddir 'brainamp2.sfp']);
data_noica.elec = ...                           % for consistency, data is converted to 'mm'
    ft_convert_units(data_noica.elec, 'mm');

cfg                 = [];                       % the next lines provide a neighbour structure
cfg.method          = 'triangulation';               % for every channel in the dataset according
cfg.template        = 'elec1005_neighb.mat';    % to the standard differences provided in the
cfg.channel         = 'EEG';                    % template file; same definition will be
cfg.trials          = 'all';
cfg.feedback        = 'no';                     % used in both conditions (hold/rest)
cfg.neighbourdist   = dist;                       % in mm in this case
neighbours          = ...
    ft_prepare_neighbours(cfg, data_noica);
save(strcat(loaddir, 'neighbours.mat'), ...     % saves data, if not present already
    'neighbours', '-v7.3');

if flag_check == 1                              % possibility to visualise the neighbouring channels in thee fieldtrip toolbox
    cfg = [];
    cfg.neighbours  = neighbours;
    ft_neighbourplot(cfg, data_noica)
end