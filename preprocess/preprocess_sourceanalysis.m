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

%% Check for SST in the template folder
if ~exist(fullfile(wdir, 'templateMRI', 'MRItrem_template0.nii.gz'), 'file')  % Looks for study specific template
    fprintf('\n\t no study specific template available in folder ~/data/templateMRI/. Please run /helpers/CreateStudySpecificTemplate.sh!\n');
    return
end

%% General settings
cd(wdir)
loaddir         = fullfile(ROOTDIR, 'mri');                     %#ok<NASGU>
load(fullfile(wdir, 'patdat.mat'));                             %#ok<LOAD>  % this file loads the meta data
if strcmp(type, 'p'); tolom = subj{2}; else; tolom = subj{1}; end           % selects whether (p) pateints or controls (c) are analysed (see (type))
debug           = 0;                                            %#ok<NASGU>

% Define and create necessary inout/output directories (outdir)
dir_mrifiles    = dir(fullfile(wdir, 'raw_MRI'));                           % input directory
outdir_mriprocessed = fullfile(wdir, 'mri_preprocessed');                   % directory in which data will be saved after processing 'raw MRI'
if ~exist(outdir_mriprocessed, 'dir'), mkdir(outdir_mriprocessed); end      % create directory, if not present already

% Start preparations for source analyses
for np = 1:numel(tolom)
    temp = control; seq = 'subject';                                        % get metadata correctly along with next line
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
        'MRItrem_template0.nii.gz'), wdir);                                 % creates template grid
    %   'mni_icbm152_t1_tal_nlin_sym_09c.nii'), wdir);                      % creates template grid according
    end
    
    %% Prepare MRI for further processing
    filename_mriprocessed = fullfile(outdir_mriprocessed, ...
        sprintf('segmentedMRI_%s.mat', code_parti));    % this is needed to check for existance to avoid redundancy
    if ~exist(filename_mriprocessed, 'file')
        filename_erpfinal = fullfile(wdir, 'data_final', ...
            sprintf('data_final_erp_%s.mat', code_parti));                  % this is needed to merge MRI and sensor data
        
        idx_mri = find(~cellfun(@isempty, regexp({dir_mrifiles.name}, ...
            sprintf('^(MRItrem_template0tANAT_)+(%s)+[0-9A-Za-z_-]*.nii.gz', ...
            num2str(code_mri)), 'match')));
        mri_subj = fullfile(wdir, 'raw_MRI', 'MRItrem_template0.nii.gz'); %#ok<NASGU>
        
        if ~isempty(idx_mri)
            mri_subj = dir_mrifiles(idx_mri).name;
        else
            fprintf('\n\t no MRI for subj: %s, use template\n', code_parti);
            continue
        end
        
        %% Run necessary steps to get MRI data processed
        prepare_MRI(code_parti, filename_erpfinal, mri_subj)
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