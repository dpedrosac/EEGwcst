function [wdir, ROOTDIR] = EEGwcst_defaults(opt)

%   This is a set of default settings intended to facilitate the scripts
%   applied in this project

%   ## Version 1.0

%   Copyright (C) June 2021
%   D. Pedrosa University Hospital of Gießen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.

if nargin == 0
    opt = 0;
end

if opt == 0
    restoredefaultpath
    close all; clear; clc;
end

global ft_default
ft_default.showlogo = 'no';

if isunix
    ROOTDIR = '/media/storage/skripte/lambda/';                             % adds the folder with all scripts to the wdir
    wdir = fullfile(ROOTDIR, 'data');                                        % defines the working directory
    addpath('/opt/fieldtrip/'); ft_defaults               % set fieldtrip defaults
elseif ispc
    if strcmp(getenv('username'), 'dpedr')
        ROOTDIR = 'D:\skripte\lambda\';
        wdir = fullfile(ROOTDIR, 'data');                                        % defines the working directory
        addpath('D:\skripte\fieldtrip'); ft_defaults;
    else
        warning("Please specific folders to 'EEGwcst_default.m' fitting to your settings")
    end
end
addpath(genpath(ROOTDIR));
