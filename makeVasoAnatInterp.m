% makeVasoAnat.m
%
%      usage: makeVasoAnatInterp()
%         by: eli merriam
%       date: 10/25/18
%    purpose: 
%
function retval = makeVasoAnatInterp()

% check arguments
if ~any(nargin == [0])
  help makeVasoAnat
  return
end

% create a new view and get some info
v = newView;
nScans = viewGet(v, 'nScans', 'BoldMC');
hdr = viewGet(v, 'niftihdr', 1);
scaleFactor = 3;

for iScan = 1 %:nScans
  % load the bold scans
  v = viewSet(v, 'curGroup', 'BoldMC'); 
  M = viewGet(v, 'transforms');
  v = viewSet(v, 'curGroup', 'Bold');   
  tSeriesBold = loadTSeries(v, iScan);
  dims = size(tSeriesBold);
  interpBold = zeros(dims(1)*scaleFactor, dims(2)*scaleFactor, dims(3)*scaleFactor, length(M));

  disppercent(-inf, 'Interpolating BOLD scans');
  for iFrame=1:length(M)
    thisM = M{iFrame};
    thisM = thisM*diag([1./[scaleFactor scaleFactor scaleFactor] 1]);
    %interpBold(:,:,:,iFrame) = warpAffine3(tSeriesBold(:,:,:,iFrame), thisM, NaN, 0, 'cubic', dims(1:3)*scaleFactor);
    interpBold(:,:,:,iFrame) = interpVolume(tSeriesBold(:,:,:,iFrame), thisM, dims(1:3)*scaleFactor);
    disppercent(iFrame/length(M));
  end
  disppercent(inf);
  
  % load the nulled scans
  v = viewSet(v, 'curGroup', 'NulledMC');    
  M = viewGet(v, 'transforms');
  v = viewSet(v, 'curGroup', 'Nulled');    
  tSeriesNulled = loadTSeries(v, iScan);
  dims = size(tSeriesNulled);
  interpNulled = zeros(dims(1)*scaleFactor, dims(2)*scaleFactor, dims(3)*scaleFactor, length(M));
  disppercent(-inf, 'Interpolating Nulled scans');
  for iFrame=1:length(M)
    thisM = M{iFrame};
    thisM = thisM*diag([1./[scaleFactor scaleFactor scaleFactor] 1]);
    %interpNulled(:,:,:,iFrame) = warpAffine3(tSeriesNulled(:,:,:,iFrame), thisM, NaN, 0, 'cubic', dims(1:3)*scaleFactor);
    interpNulled(:,:,:,iFrame) = interpVolume(tSeriesNulled(:,:,:,iFrame), thisM, dims(1:3)*scaleFactor);
   disppercent(iFrame/length(M));
  end
  disppercent(inf);
  
  % concat bold and nulled tSeries together
  tSeries = cat(4, interpBold, interpNulled);
  % compute 1 / (std/abs(mean))
  cvar(:,:,:,iScan) = 1 / (nanstd(tSeries,[],4) ./ abs(nanmean(tSeries,4)));
end

% average across runs
vasoAnat = nanmedian(cvar,4);

% trim bad values
vasoAnat(isinf(vasoAnat)) = 0;
vasoAnat(isnan(vasoAnat)) = 0;
vasoAnat(vasoAnat>20) = 0;

% set the qform/sform
hdr = cbiSetNiftiQform(hdr,hdr.qform44*diag([1./[scaleFactor scaleFactor scaleFactor] 1]));
hdr = cbiSetNiftiSform(hdr,hdr.sform44*diag([1./[scaleFactor scaleFactor scaleFactor] 1]));
hdr.dim(2:4) = [size(vasoAnat)];


% % upsample some more
% scaleFactor = 3;
% xform = diag([1./[scaleFactor scaleFactor scaleFactor 1]]);
% vasoAnat = warpAffine3(vasoAnat, xform, NaN, 0, 'linear', size(vasoAnat)*scaleFactor);
% hdr = cbiSetNiftiQform(hdr,hdr.qform44*diag([1./[scaleFactor scaleFactor scaleFactor] 1]));
% hdr = cbiSetNiftiSform(hdr,hdr.sform44*diag([1./[scaleFactor scaleFactor scaleFactor] 1]));
% hdr.dim(2:4) = [size(vasoAnat)];
% 
% % blur a bit to deal with interpolation error
% for iSlice = 1:size(vasoAnat,3);
%   vasoAnat(:,:,iSlice) = blur(vasoAnat(:,:,iSlice));
% end

% set the file extension
niftiFileExtension = '.nii';

% Path
pathStr = fullfile(viewGet(view,'anatomydir'),['vasoAnatInterp',niftiFileExtension]);

% Write, though check for over-writing
saveFlag = 'Yes';
if exist([pathStr,'.mat'],'file')
    if confirm
        saveFlag = questdlg([pathStr,' already exists. Overwrite?'],...
            'Save Overlay?','Yes','No','No');
	end
end
if strcmp(saveFlag,'Yes')
    fprintf('Saving %s...\n',pathStr);
    [byteswritten,hdr] = cbiWriteNifti(pathStr, vasoAnat, hdr);
    fprintf('done\n');
else
    fprintf('Anatomy not saved...');
end

% Delete temporary views
deleteView(v);

return;


