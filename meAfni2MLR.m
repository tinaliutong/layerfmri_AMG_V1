% meAfni2MLR.m
%
%        $Id$
%      usage: meAfni2MLR
%         by: eli merriam
%       date: 04/01/2019
%    purpose: copy the output of afni's multi echo pipeline to MLR
%             run this from the directory that you ran mePreproc from
function [] = meAfni2MLR()

if ~exist('afnipreproc')
    disp('UHOH: afnipreproc does not exist');
    return;
end

% move into the preproc directory
cd('afnipreproc');

% make MLR directories
if ~exist('./Raw'); mkdir('Raw/'); end
if ~exist('./Raw/TSeries'); mkdir('Raw/TSeries'); end
if ~exist('./Etc'); mkdir('Etc'); end
if ~exist('./Anatomy'); mkdir('Anatomy'); end

% copy the anatomy file to Anatomy
if ~isempty(dir('*anat*'))
    system('cp *anat*.nii Anatomy');
end

if ~isempty(dir('*pure*'))
    system('cp *pure*.nii Anatomy');
end

if ~isempty(dir('*volregMean*'))
    system('cp *volregMean*.nii Etc');
end

% get the directory name
preprocdir = sprintf('%s.results', getLastDir(pwd));

% now copy the 'combine' files to Raw
combFiles = dir(sprintf('%s/*combine+orig.HEAD', preprocdir));

for iCombFile = 1:length(combFiles)

    % get run number
    sep = strfind(combFiles(iCombFile).name, 'combine');
    runNum = str2num(combFiles(iCombFile).name(sep-3:sep-2));

    [status, result] = system(sprintf('3dinfo -history %s/%s | grep mr_', combFiles(iCombFile).folder, combFiles(iCombFile).name));
    nameStart = strfind(result, 'mr_0');
    nameEnd = strfind(result, 'chan')-2;
    cmdStr = sprintf('3dcopy %s/%s Raw/TSeries/%s_aproc_r%02i.nii', combFiles(iCombFile).folder, combFiles(iCombFile).name, result(nameStart:nameEnd), runNum);
    disp(cmdStr);
    system(cmdStr);

    % get a listing of the *.1D files
    onedfiles = dir(sprintf('%s/*r%02i*.1D', preprocdir, runNum));
    for iDfile = 1:length(onedfiles)
        cmdStr = sprintf('cp %s/%s Etc/%s_%s', onedfiles(iDfile).folder, onedfiles(iDfile).name, result(nameStart:nameEnd), onedfiles(iDfile).name);
        disp(cmdStr);
        system(cmdStr);
    end       
end

% copy session-wide params
system(sprintf('cp %s/dfile_rall.1D Etc/', preprocdir));
system(sprintf('cp %s/outcount_rall.1D Etc/', preprocdir));
system(sprintf('cp %s/motion_*_enorm.1D Etc/', preprocdir));


% clean up
system(sprintf('mv %s.results afnipreproc', getLastDir(pwd)));
system(sprintf('mv proc.%s afnipreproc', getLastDir(pwd)));
system(sprintf('mv output.proc.%s afnipreproc', getLastDir(pwd)));

system('mv GERT_Reco_dicom* afnipreproc');
system('mv dimon.files.run* afnipreproc');
system('mv *nii afnipreproc');
