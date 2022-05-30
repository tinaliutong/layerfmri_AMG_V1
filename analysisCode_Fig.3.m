%% prep: download the data, unzip the folder 
% 0: caveat: don't do motion correction before splitting bold and vaso
% 1: dovasoatNIH('splitVaso');
% (cbiWriteNifti) Scaling data from single to float32
% (viewSet) framePeriod set to: 4.825000
% set scanParams.framePeriod = 2.4125
% 2: run full motion comp on Nulled (110) and Bold (111) groups
% make a new group: creating nulledMC and boldMC
% 3: dovasoatNIH('upSampleBold')
% 4: dovasoatNIH('upSampleNulled')
% 5: junkFrame 
    % under upsampleBold and upsampleNulled
    % go to Group/Edit Group, num Frames
    % for 221 volumes,  enter 2 in the Junk frames, then enter 217 in the
    % Num frames. This methods throw out the first and last 2 farmes. 
% 6: dovasoatNIH('computeVaso2')
% 7: concatenate across vaso runs 
    % edit VASO group and change junkframe to 0, total frame = 217 for a total of 221 volume)
% 8: concatenate across upSampleBold runs 
% 9: do regular GLM analysis %  
% 10: make VasoAnatomy
% 11: mrQuit then mrAlign (for benson atlas)
% 12: re makevasoAnat and makeVasoAnatInterpZ (upsampled VasoAnatomy)

%%
function[] =  dovasoatNIH(fignum, varargin)

getArgs(varargin, [], 'verbose=0');

switch (fignum)
    case {'splitVaso'}
        splitVaso();
    case {'upSampleBold'}
        upSampleBold();
    case {'upSampleNulled'}
        upSampleNulled();
    case {'computeVaso2'}
        computeVaso2();
    otherwise
        disp(sprintf('You need to specifiy an analysis to run'))
end
end

function splitVaso()

v = newView;
nScans = viewGet(v, 'nScans');

v = viewSet(v, 'newGroup', 'Bold');
v = viewSet(v, 'newGroup', 'Nulled');

for iScan = 1:nScans
    % load the data from raw
    % iScan = 1
    v = viewSet(v, 'curGroup', 'Raw');
    tSeries = loadTSeries(v, iScan);
    
    % subset the correct frames
    bold  =  tSeries(:,:,:,1:2:end);
    nulled = tSeries(:,:,:,2:2:end);
    
    % load the scan params and nifti header
    scanParams = viewGet(v, 'scanparams', iScan);
    %scanParams.framePeriod = scanParams.framePeriod*2;
    scanParams.framePeriod = 2.4125;

    % load the nifti header
    hdr = viewGet(v, 'niftihdr', iScan);
    % adjust the TR
    % hdr.pixdim(5) = hdr.pixdim(5)*2; = 2.7366x2 
    hdr.pixdim(5) = 2.4125*2;

    v = viewSet(v, 'curGroup', 'Nulled');
    scanParams.description = [scanParams.description ' nulled'];
    scanParams.nFrames = size(nulled,4);
    v = saveNewTSeries(v, nulled, scanParams, hdr);
    v = viewSet(v, 'auxparam', 'volTrigRatio', 1/2, iScan);
    
    v = viewSet(v, 'curGroup', 'Bold');
    scanParams.description = [scanParams.description ' bold'];
    scanParams.nFrames = size(bold,4);
    v = saveNewTSeries(v, bold, scanParams, hdr);
    v = viewSet(v, 'auxparam', 'volTrigRatio', 1/2, iScan);
    
end
deleteView(v);

end

function  upSampleBold()
v = newView;
v = viewSet(v, 'newGroup', 'upsampledBold');
v = viewSet(v, 'curGroup', 'boldMC');
%v = viewSet(v, 'curGroup', 'MotionComp');

nScans = viewGet(v, 'nscans');
for iScan = 1:nScans
    v = viewSet(v, 'curGroup', 'boldMC');
   %v = viewSet(v, 'curGroup', 'MotionComp');

    tSeries = loadTSeries(v, iScan);
    %frameperiod = 5;
    frameperiod = viewGet(v, 'frameperiod');

    % original size
    [Nx Ny Nz Nt] = size(tSeries);
    
    % upsample in time
    % for odd number of volumes
    [xgrid,ygrid,zgrid,tgrid] = ndgrid(1:Nx,1:Ny,1:Nz,1:0.5:Nt);
    % for even number of volumes
    % [xgrid,ygrid,zgrid,tgrid] = ndgrid(1:Nx,1:Ny,1:Nz,1:0.5:Nt+0.5);
    
    tSeries = interpn(tSeries, xgrid, ygrid, zgrid, tgrid, 'cubic');
       
    % get the nifti header and scan params
    hdr = viewGet(v, 'niftihdr', 1);
    hdr.dim(5) = size(tSeries,4);
    hdr.pixdim(5) = hdr.pixdim(5)/2;
    
    scanParams = viewGet(v, 'scanparams', 1);
    scanParams.description = 'upsampledBold';
    scanParams.fileName = [];
    scanParams.totalFrames = size(tSeries,4);
    % switch the the upsampledBold group
    v = viewSet(v, 'curGroup', 'upsampledBold');
    
    % write out the new tSeries
    v = saveNewTSeries(v, tSeries, scanParams, hdr);
    v = viewSet(v, 'auxparam', 'volTrigRatio', 1, iScan);
    v = viewSet(v, 'frameperiod', frameperiod/2, iScan);
    saveSession();
    
end
deleteView(v);

end

function  upSampleNulled()
v = newView;
v = viewSet(v, 'newGroup', 'upsampledNulled');
v = viewSet(v, 'curGroup', 'nulledMC');
%v = viewSet(v, 'curGroup', 'MotionComp');

nScans = viewGet(v, 'nscans');
for iScan = 1:nScans
    v = viewSet(v, 'curGroup', 'nulledMC');
   % v = viewSet(v, 'curGroup', 'MotionComp');

    tSeries = loadTSeries(v, iScan);
    frameperiod = viewGet(v, 'frameperiod');
    
    % original size
    [Nx Ny Nz Nt] = size(tSeries);
    
    % upsample in time
    % for odd number of volumes
     [xgrid,ygrid,zgrid,tgrid] = ndgrid(1:Nx,1:Ny,1:Nz,0.5:0.5:Nt+0.5);
    % for even number of volumes
    % [xgrid,ygrid,zgrid,tgrid] = ndgrid(1:Nx,1:Ny,0.5:Nz,0.5:0.5:Nt);
    tSeries = interpn(tSeries, xgrid, ygrid, zgrid, tgrid, 'cubic');
    %tSeries = permute(tSeries, [2 1 3 4]);
    
    % get the nifti header and scan params
    hdr = viewGet(v, 'niftihdr', 1);
    hdr.dim(5) = size(tSeries,4);
    hdr.pixdim(5) = hdr.pixdim(5)/2;
    
    scanParams = viewGet(v, 'scanparams', 1);
    scanParams.description = 'upsamplednulled';
    scanParams.fileName = [];
    scanParams.totalFrames = size(tSeries,4);
    % switch the the Vaso group
    v = viewSet(v, 'curGroup', 'upsampledNulled');
    
    % write out the new tSeries
    v = saveNewTSeries(v, tSeries, scanParams, hdr);
    v = viewSet(v, 'auxparam', 'volTrigRatio', 1, iScan);
    v = viewSet(v, 'frameperiod', frameperiod/2, iScan);
    saveSession();
    
end
deleteView(v);

end

function  computeVaso2()
v = newView;
v = viewSet(v, 'newGroup', 'Vaso');
nScans = viewGet(v, 'nScans');

%for iScan = [1:7]
for iScan = 1:nScans   % iScan = 1
    v = viewSet(v, 'curGroup', 'upSampleBold');
    dBold = loadTSeries(v, iScan);
    dims = size(dBold);
    dBold = reshape(dBold, dims(1)*dims(2)*dims(3), dims(4));
    junkFrames = viewGet(v, 'junkFrames', iScan);
    nFrames = viewGet(v, 'nFrames', iScan);
    dBold = dBold(:,junkFrames+1:nFrames+junkFrames);
    
    v = viewSet(v, 'curGroup', 'upsampleNulled');
    dNulled = loadTSeries(v, iScan);
    dims = size(dNulled);
    dNulled = reshape(dNulled, dims(1)*dims(2)*dims(3), dims(4));
    junkFrames = viewGet(v, 'junkFrames', iScan);
    nFrames = viewGet(v, 'nFrames', iScan);
    dNulled = dNulled(:,junkFrames+1:nFrames+junkFrames);
    
    if size(dBold,2) ~= size(dNulled,2);
        disp(sprintf('UHOH: size of bold and nulled do not match'));
        keyboard
    end
    
    % loop over frames
    vaso = zeros(size(dNulled));
    for iFrame = 1:size(dNulled,2);
        vaso(:,iFrame) = dNulled(:,iFrame) ./ dBold(:,iFrame);
    end    
    vaso = reshape(vaso, dims(1), dims(2), dims(3), size(dNulled,2));
    
    % get the nifti header and scan params
    hdr = viewGet(v, 'niftihdr', 1);
    scanParams = viewGet(v, 'scanparams', 1);
    scanParams.description = 'vaso';
    scanParams.fileName = [];
    scanParams.totalJunkedFrames = junkFrames;
    
    % switch the the Vaso group
    v = viewSet(v, 'curGroup', 'Vaso');
    
    % write out the new tSeries
    v = saveNewTSeries(v, vaso, scanParams, hdr);
    v = viewSet(v, 'auxparam', 'volTrigRatio', 1, iScan);
    
end

saveSession;

deleteView(v);
end
