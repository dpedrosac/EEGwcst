function preprocess_sourceanalysis(subj, ROOTDIR, wdir, type)

%   Wrapper function in order to perform all necessary preprocessing steps
%   for the creation of necesary source analyses.

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
cd(wdir)
loaddir     = fullfile(ROOTDIR, 'mri');                         %#ok<NASGU>
load(fullfile(wdir, 'patdat.mat'));                             %#ok<LOAD>  % this file loads the meta data
if strcmp(type, 'p'); tolom = subj{2}; else; tolom = subj{1}; end           % selects whether (p) pateints or controls (c) are analysed (see (type))
debug       = 0;                                                %#ok<NASGU>

% Define and create necessary output directories (outdir)
outdir_mriprocessed = fullfile(wdir, 'mri_preprocessed');                          % directory in which data will be saved after processing 'raw MRI'
if ~exist(outdir_mriprocessed, 'dir'), mkdir(outdir_mriprocessed); end                    % create directory, if not present already

% Start with preparations for source analyses
for np = 1:numel(tolom)
    temp = control; seq = 'subject';
    if strcmp(type, 'p'); temp = patient; seq = 'patient'; end
    fprintf('the %s actually being computed is %s\n', seq, temp(np).name)
    code_parti = upper(temp(np).code);                                      % relevant information for later in the next few lines
    code_mri = upper(temp(np).mrt);
    if isempty(code_mri); continue; end %TODO: generic MRI (template?!) should be used
    
    
    %% Create template grid for source analyses later
    filename_template_grid = fullfile(wdir, 'templateMRI', ...
        'template_grid.mat');
    if ~exist(filename_template_grid, 'file')
        create_template_grid(fullfile(wdir, 'templateMRI', ...
            'mni_icbm152_t1_tal_nlin_sym_09c.nii'), wdir);                  % creates a template grid according
    end
    
    %% Prepare MRI for further processing
    filename_mriprocessed = fullfile(outdir_mriprocessed, ...
        sprintf('MRIprocessed_%s.nii.gz', code_parti));                     % this is needed to check for existance to avoid redundancy
    if ~exist(filename_mriprocessed, 'file')
        filename_erpfinal = fullfile(wdir, 'data_final', ...
            sprintf('data_final_erp_%s.mat', code_parti));                  % this is needed to merge MRI and sensor data
        
        dir_mrifilenames = dir(fullfile(wdir, 'raw_MRI'));
        idx_mri = find(~cellfun(@isempty, regexp({dir_mrifilenames.name}, ...
            sprintf('^(MRItrem_template0tANAT_)+(%s)+[0-9A-Za-z_-]*.nii.gz', ...
            num2str(code_mri)), 'match')));
        mri_file = fullfile(wdir, 'raw_MRI', 'MRItrem_template0.nii.gz'); %#ok<NASGU>
        if ~isempty(idx_mri)
            mri_file = dir_mrifilenames(idx_mri).name;
        else
            fprintf('\n\t no MRI available for subj: %s, using template\n', code_parti);
            continue
        end % use subject specific MRI if avialble else skip
        
        prepare_MRI(code_parti, filename_erpfinal, mri_file, ...
            filename_mriprocessed)
    else
        fprintf('\n\t MRI already prepared for subj: %s !\n', code_parti);
    end
    
    %% Remove channels with artefacts from data after visual inspection
    %     filename_clean = fullfile(outdir_clean, ...
    %         strcat('data_clean_', file_prefix, '_', conds{c}, ...
    %         sprintf('_%s', experiment), '.mat'));                             % removes the artefacts detected in the ICA
    %     if ~exist(filename_clean, 'file')
    %         filename_noica = strcat('data_noica_', file_prefix, '_', ...
    %             conds{c}, sprintf('_%s', experiment), '.mat');                  % filename under which data will be saved in the (outdir) directory
    %         filename_noica = fullfile(wdir, 'data_noica', filename_noica);
    %         remove_badchannels(subj{s}, conds{c}, filename_noica, ...
    %             filename_clean);
    %     else
    %         fprintf('\n\t interpolation of bad channels already finished for subj: %s in the %s (Hz) condition!\n', subj{s}, conds{c});
    %     end
    
    
    
end