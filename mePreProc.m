function [] = mePreProc(filename)

minFileSize = 1000000; %in bytes. If less than 4 TRs will get error. Zvi

% identify the ME file
meFiles = dir('*chan_001.nii');

% identify the blip file
blipFiles = dir('*blip*');
blipFiles = blipFiles(end);

% identify the anatomy file
anatFiles = dir('*anat*');
if isempty(anatFiles)
    anatFiles = dir('*pure*');
end

% create the command string for the list of datasets
dsetString = [];
counter = 1;
for iFile = 1:length(meFiles)
    tempStr = meFiles(iFile).name;
    sep = strfind(tempStr, '.');
    tempStr = tempStr(1:sep-2);

    if meFiles(iFile).bytes > minFileSize
        dsetString = cat(2, dsetString, sprintf(' -dsets_me_run %s?.nii', tempStr));
        % add this filename to the list
        listOfOriginalFilenames{counter} = meFiles(iFile).name;
        counter = counter+1;
    end
    
end
save('origFilenames.mat','listOfOriginalFilenames');

% create the command string for the blip correction
sep = strfind(blipFiles.name, '_');
blipNum = str2num(blipFiles.name(sep(1)+1:sep(2)-1));
ind = [];
for iFile = 1:length(meFiles);    
    sep = strfind(meFiles(iFile).name, '_');
    meNum = str2num(meFiles(iFile).name(sep(1)+1:sep(2)-1));
    ind(iFile) = abs(meNum - blipNum);    
end

[foo,ind] = min(ind);
blipForwardDset = meFiles(ind).name;
blipString = sprintf('-blip_forward_dset %s -blip_reverse_dset %s', blipForwardDset, blipFiles.name);
%blipString = [];

% create the command string for volreg
for iFile=1:length(anatFiles)
    sep = strfind(anatFiles(iFile).name, '_');
    anatNum(iFile) = str2num(anatFiles(iFile).name(sep(1)+1:sep(2)-1));
end

ind = [];
for iFile = 1:length(meFiles);    
    sep = strfind(meFiles(iFile).name, '_');
    meNum = str2num(meFiles(iFile).name(sep(1)+1:sep(2)-1));
    ind(iFile) = abs(meNum - min(anatNum));    
end
[foo,ind] = min(ind);


% run 3dvolreg
volRegBase = meFiles(ind).name;
% compute the last frame number
h = cbiReadNiftiHeader(volRegBase);
lastFrame = h.dim(5)-1;
volRegStr = sprintf('3dvolreg -input %s -prefix volreg.nii -base %i', volRegBase, lastFrame);
system(volRegStr);
% get the mean over time
system('3dTstat -prefix volregMean.nii volreg.nii');
volregString = sprintf('-volreg_base_dset volregMean.nii');

subjID = getLastDir(pwd);

doit = sprintf('afni_proc.py -subj_id %s %s -reg_echo 1 -echo_times 14.2 30.1 46 %s -blocks tshift volreg combine -tcat_remove_first_trs 4 %s -combine_method OC -execute', subjID, dsetString, blipString, volregString);

%doit = sprintf('afni_proc.py -subj_id %s %s -reg_echo 1 -echo_times 14.2 30.1 46 %s -blocks tshift volreg mask combine -tcat_remove_first_trs 4 %s -combine_method tedana -execute', subjID, dsetString, blipString, volregString);
%keyboard
system(doit); %!!!



