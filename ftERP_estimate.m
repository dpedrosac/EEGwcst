function ERP_avg = ftERP_estimate(data_final, trlsnum)
%   This function estimates ERP using the fieldtrip toolbox with the
%   settings defined below and saves the data in a file

%   Copyright (C) March 2018
%   D. Pedrosa, University Hospital of Gieﬂen and Marburg

%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.


cfg = [];
cfg.trials      = find(ismember(data_final.trialinfo, trlsnum));            % use only the
cfg.feedback    = 'no';
data_select     = ft_selectdata(cfg,data_final);

cfg = [];
cfg.keeptrials          = 'yes';
cfg.covariance          = 'yes';
cfg.covariancewindow    = [-inf -.05];
cfg.feedback            = 'no';
ERP_avg                 = ft_timelockanalysis(cfg, data_select);
end