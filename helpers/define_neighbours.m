function neighbours = ...
    define_neighbours(data, loaddir, dist, flag_check)

%   This function runs the routine to obtain the neighbours later needed to
%   interpolate the 'bad channels'

%   Copyright (C) January 2018
%   D. Pedrosa, University Hospital of Gieﬂen and Marburg

%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

data.elec = ft_read_sens([loaddir '\brainamp_montage.sfp']);                % this block assigns the electrode positions according to BRAINVISION to the data
data.elec = ft_convert_units(data.elec, 'mm');                        % for consistency, data is converted to 'mm'

cfg                 = [];                                                   % the next lines provide a neighbour structure
cfg.method          = 'template';                          % for every channel in the dataset according
cfg.template        = 'elec1005_neighb.mat';                                % to the standard differences provided in the
cfg.channel         = 'EEG';                                                % template file; same definition will be
%cfg.trials          = 'all';
cfg.feedback        = 'no';                                                 % used in both conditions (hold/rest)
% cfg.neighbourdist   = 4;                                                 % in mm in this case
neighbours          = ft_prepare_neighbours(cfg, data);
save(strcat(loaddir, '\neighbours.mat'), 'neighbours', '-v7.3');            % saves data to HDD
    

if flag_check == 1                                                          % visualise neighbouring channels
    cfg             = [];
    cfg.neighbours  = neighbours;
    ft_neighbourplot(cfg, data)
end