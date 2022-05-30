% makeTemplateOverlays
%
% Takes a list of subjects and creates Benson atlas, converts the files
% into nifty, and creates group called templates for visualizing/filterning
% by eccentricity and phase.
%
% mr 7/2018



% step 0: Get Benson atlas for all subjects
subList = {'s0052'};
% 's0005b','s0016b','s0008', 's0040','s0005', 's0022', 's0032', 's0039', 's0041', 's0045', 's0042', 's0046', 's0049'

for iSub = 1:length(subList)
	
	subNum = subList{iSub};
	% Invert the right hemisphere
	system(sprintf('xhemireg --s %s',subNum));
	% Register the left hemisphere to fsaverage_sym
	system(sprintf('surfreg --s %s --t fsaverage_sym --lh',subNum));
	% Register the inverted right hemisphere to fsaverage_sym
	system(sprintf('surfreg --s %s --t fsaverage_sym --lh --xhemi',subNum));
	
	system(sprintf('docker run -ti --rm -v $SUBJECTS_DIR/%s:/input nben/occipital_atlas:latest', subNum));
	
	surfPath = sprintf('/misc/data58/merriamep/data/freesurfer/%s/mri/',subNum);
	cd(surfPath);
	
	% step 1: convert mgz files to nifti
	
	system('mri_convert native.template_eccen.mgz ../surfRelax/template_eccen.nii');
	system('mri_convert native.template_angle.mgz ../surfRelax/template_angle.nii');
	system('mri_convert native.template_areas.mgz ../surfRelax/template_areas.nii');
	
	cd('/misc/data58/merriamep/data/freesurfer/');
end


for subjNum = length(subjDirs) %1:length(subjDirs)
	mrQuit
	fullSubSess = subjDirs{subjNum};
	cd(fullSubSess)
	
	% step 2: create a new mrLoadRet group for the templates
	v = newView;
	v = viewSet(v, 'newGroup', 'templates');
	v = viewSet(v, 'curGroup', 'templates');
	
	nScans = viewGet(v, 'nScans');
	
	% Step 3: create a 'tSeries', this is just for header information...
	subject = sprintf('s%s',fullSubSess(2:5)) 
    % for 7T directory, subject = sprintf('s%s',fullSubSess(2:5),'_7T');
	pathtotemplate = fullfile('/misc/data58/merriamep/data/freesurfer/', subject, 'surfRelax');
	
	if nScans == 0
		[v,params] = importTSeries(v,[],'justGetParams=0','defaultParams=1',sprintf('pathname=%s/template_eccen.nii', pathtotemplate));
	end
	
	iScan = 1;	
	s = viewGet(v, 'stimfile', iScan);
	groupNum = viewGet(v, 'curGroup');
	mrSetPref('overwritePolicy','Merge');
	
    % subject = 's0070'
	d = cbiReadNifti(fullfile('/misc/data58/merriamep/data/freesurfer/', subject, 'surfRelax', 'template_eccen.nii'));
	[v templateRet] = mrDispOverlay(d,iScan,groupNum,v,'overlayName=eccen','saveName=templateRet', 'cmap', hsv);
	
	d = cbiReadNifti(fullfile('/misc/data58/merriamep/data/freesurfer/', subject, 'surfRelax', 'template_angle.nii'));
	[v templateRet] = mrDispOverlay(d,iScan,groupNum,v,'overlayName=angle','saveName=templateRet', 'cmap', hsv);
	
	d = cbiReadNifti(fullfile('/misc/data58/merriamep/data/freesurfer/', subject, 'surfRelax', 'template_areas.nii'));
	[v templateRet] = mrDispOverlay(d,iScan,groupNum,v,'overlayName=areas','saveName=templateRet', 'cmap', flag);
	
	deleteView(v);

	cd(amypath);
end


