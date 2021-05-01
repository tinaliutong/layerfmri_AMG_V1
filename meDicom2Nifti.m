function [] = meDicom2Nifti()

% first step is to unzip the 
tarfile = dir('*tgz');
if ~isempty(tarfile)
    commandString = sprintf('tar xfz %s', tarfile(1).name);
else
    tarfile = dir('*tar.gz');
    if ~isempty(tarfile)
        commandString = sprintf('tar xfz %s', tarfile(1).name);
    else
        tarfile = dir('*tar');
        if ~isempty(tarfile)
            commandString = sprintf('tar xf %s', tarfile(1).name);
        end
    end
end

if isempty(tarfile)
    disp(sprintf('(meDicom2Nifti) Uhoh: There is no tar file here'));
    return;
end

disp(sprintf('(meDicom2Nifti) Unpacking tar file...'));
system(commandString)

% second step is to convert to nifti file format

% parse the elements of the file
tarfile = tarfile.name;
sep = findstr(tarfile, '-');
subName = tarfile(1:sep(1)-1);
mrnNum = tarfile(sep(1)+1:sep(2)-1);
studyDate = tarfile(sep(2)+1:sep(3)-1);
randNum = tarfile(sep(3)+1:sep(4)-1);
fileSuf = tarfile(sep(4)+1:end);

scanDirs = dir(sprintf('%s-%s/%s-%s', subName, mrnNum, studyDate, randNum));


disp(sprintf('(meDicom2Nifti) Looping over runs, converting to nifti format...'));
for iScan = 1:length(scanDirs)
    if strfind(scanDirs(iScan).name, 'mr_')
        % get the list of file names
        scanFileNames = dir(sprintf('%s/%s/*dcm', scanDirs(iScan).folder, scanDirs(iScan).name));
        % check what sort of file we have

        sep = findstr(scanFileNames(1).name, '_');
        dicomfolder = getLastDir(scanFileNames(1).folder);
        dicomname = scanFileNames(1).name(1:sep(1)-1);
        niftiprefix = sprintf('%s_%s', dicomfolder, dicomname);  
        
        if strfind(scanFileNames(1).name, 'anat')  
            % we have a phase encode reversed file
            commandString = sprintf('Dimon -infile_prefix %s/%s -no_wait -sort_method geme_index -gert_create_dataset -gert_write_as_nifti  -gert_to3d_prefix %s',  scanFileNames(1).folder, dicomname, niftiprefix);
            system(commandString);            

        elseif strfind(scanFileNames(1).name, 'rage')
            % we have a phase encode reversed file
            commandString = sprintf('Dimon -infile_prefix %s/%s -no_wait -sort_method geme_index -gert_create_dataset -gert_write_as_nifti  -gert_to3d_prefix %s',  scanFileNames(1).folder, dicomname, niftiprefix);
            system(commandString);            

        elseif strfind(scanFileNames(1).name, 'me3')
            % we have a real multiecho dataset
            commandString = sprintf('Dimon -infile_prefix %s/%s -no_wait -num_chan 3 -sort_method geme_index -gert_create_dataset -gert_write_as_nifti  -gert_to3d_prefix %s',  scanFileNames(1).folder, dicomname, niftiprefix);
            system(commandString);
            
        elseif strfind(scanFileNames(1).name, 'blip')
            % we have a phase encode reversed file
            commandString = sprintf('Dimon -infile_prefix %s/%s -no_wait -sort_method geme_index -gert_create_dataset -gert_write_as_nifti  -gert_to3d_prefix %s',  scanFileNames(1).folder, dicomname, niftiprefix);
            system(commandString);
        else
            % file type that we do not want to convert
            disp(sprintf('UHOH: You have a %s, not converting it', scanFileNames(1).name));
        end
        
    end
end





