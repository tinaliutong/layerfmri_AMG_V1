%% open mrloadRet

%% under 'GLM_all' in the faceLoc concatenation 
rois = loadROITSeries(getMLRView, []); % load all 50 lh+rh.Wang2015atlas ROIs

for i=1:length(rois)
scanCoords{i} = rois{1,i}.scanCoords;
end

%size=[160,160,42]; % scan dim for 7T BOLD
size=[64,64,30]; % scan dim for 3T BOLD

d1 = viewGet(getMLRView, 'd');
R2_data = viewGet(getMLRView, 'overlay')
R2 = R2_data.data{1}; % dim = [160,160,42] for 7T, or [64,64,30] for 3T

max(R2(:))
nanmean(R2(:))

for roi=1:numel(scanCoords)
if ~isempty(scanCoords{roi})
roi2{roi}=scanCoords{roi};
linearInd{roi} = sub2ind(size,roi2{roi}(1,:),roi2{roi}(2,:),roi2{roi}(3,:));
r2_roi2{roi}=R2(linearInd{roi});
roi2_threshold{roi}=roi2{roi}(:,r2_roi2{roi}>0.1);
end
end

for j=1:25
    roi2_threshold_hemi{j}=[roi2_threshold{j},roi2_threshold{j+25}];
end

% FEF, hV4, IPS0
for k=1:3
roi2_threshold_combinedROI2{k}=roi2_threshold_hemi{k};
end

% IPS1-5
roi2_threshold_combinedROI2{4}=[roi2_threshold_hemi{4},roi2_threshold_hemi{5},roi2_threshold_hemi{6},roi2_threshold_hemi{7},roi2_threshold_hemi{8}];
% LO = LO1+LO2
roi2_threshold_combinedROI2{5}=[roi2_threshold_hemi{9},roi2_threshold_hemi{10}];
% PHC = PHC1+PHC2
roi2_threshold_combinedROI2{6}=[roi2_threshold_hemi{11},roi2_threshold_hemi{12}];
% SPL
roi2_threshold_combinedROI2{7}=[roi2_threshold_hemi{13}];
% TO = TO1+TO2
roi2_threshold_combinedROI2{8}=[roi2_threshold_hemi{14},roi2_threshold_hemi{15}];
% V1 = V1d+V1v
roi2_threshold_combinedROI2{9}=[roi2_threshold_hemi{16},roi2_threshold_hemi{17}];
% V2 = V2d+V2v
roi2_threshold_combinedROI2{10}=[roi2_threshold_hemi{18},roi2_threshold_hemi{19}];
% V3AB = V3A+V3v
roi2_threshold_combinedROI2{11}=[roi2_threshold_hemi{20},roi2_threshold_hemi{21}];
% V3 = V3d+V3v
roi2_threshold_combinedROI2{12}=[roi2_threshold_hemi{22},roi2_threshold_hemi{23}];
% VO = VO1+VO2
roi2_threshold_combinedROI2{13}=[roi2_threshold_hemi{24},roi2_threshold_hemi{25}];


%% under 'GLM_F-N' in the amyV1 task
d2 = viewGet(getMLRView, 'd');
rois_AmyV1 = loadROITSeries(getMLRView, []); % load all 50 lh and rh.Wang2015atlas ROIs (same order!)
%rois_AmyV1 = rois

%% for hemi combined
for j=1:25
   rois_AmyV1_hemi{j}.tSeries=[rois_AmyV1{j}.tSeries; rois_AmyV1{j+25}.tSeries];
   rois_AmyV1_hemi{j}.scanCoords=[rois_AmyV1{j}.scanCoords, rois_AmyV1{j+25}.scanCoords];
end

for k=1:3
rois_AmyV1_hemi_combinedROI2{k}.tSeries = rois_AmyV1_hemi{k}.tSeries;
rois_AmyV1_hemi_combinedROI2{k}.scanCoords=rois_AmyV1_hemi{k}.scanCoords;
end

%% 13 ROIs
rois_AmyV1_hemi_combinedROI2{4}.scanCoords =[rois_AmyV1_hemi{4}.scanCoords,rois_AmyV1_hemi{5}.scanCoords,rois_AmyV1_hemi{6}.scanCoords,rois_AmyV1_hemi{7}.scanCoords,rois_AmyV1_hemi{8}.scanCoords];
rois_AmyV1_hemi_combinedROI2{5}.scanCoords =[rois_AmyV1_hemi{9}.scanCoords,rois_AmyV1_hemi{10}.scanCoords];
rois_AmyV1_hemi_combinedROI2{6}.scanCoords =[rois_AmyV1_hemi{11}.scanCoords,rois_AmyV1_hemi{12}.scanCoords];
rois_AmyV1_hemi_combinedROI2{7}.scanCoords =[rois_AmyV1_hemi{13}.scanCoords];
rois_AmyV1_hemi_combinedROI2{8}.scanCoords =[rois_AmyV1_hemi{14}.scanCoords,rois_AmyV1_hemi{15}.scanCoords];
rois_AmyV1_hemi_combinedROI2{9}.scanCoords =[rois_AmyV1_hemi{16}.scanCoords,rois_AmyV1_hemi{17}.scanCoords];
rois_AmyV1_hemi_combinedROI2{10}.scanCoords=[rois_AmyV1_hemi{18}.scanCoords,rois_AmyV1_hemi{19}.scanCoords];
rois_AmyV1_hemi_combinedROI2{11}.scanCoords=[rois_AmyV1_hemi{20}.scanCoords,rois_AmyV1_hemi{21}.scanCoords];
rois_AmyV1_hemi_combinedROI2{12}.scanCoords=[rois_AmyV1_hemi{22}.scanCoords,rois_AmyV1_hemi{23}.scanCoords];
rois_AmyV1_hemi_combinedROI2{13}.scanCoords=[rois_AmyV1_hemi{24}.scanCoords,rois_AmyV1_hemi{25}.scanCoords];

rois_AmyV1_hemi_combinedROI2{4}.tSeries =[rois_AmyV1_hemi{4}.tSeries;rois_AmyV1_hemi{5}.tSeries;rois_AmyV1_hemi{6}.tSeries;rois_AmyV1_hemi{7}.tSeries;rois_AmyV1_hemi{8}.tSeries];
rois_AmyV1_hemi_combinedROI2{5}.tSeries =[rois_AmyV1_hemi{9}.tSeries;rois_AmyV1_hemi{10}.tSeries];
rois_AmyV1_hemi_combinedROI2{6}.tSeries =[rois_AmyV1_hemi{11}.tSeries;rois_AmyV1_hemi{12}.tSeries];
rois_AmyV1_hemi_combinedROI2{7}.tSeries =[rois_AmyV1_hemi{13}.tSeries];
rois_AmyV1_hemi_combinedROI2{8}.tSeries =[rois_AmyV1_hemi{14}.tSeries;rois_AmyV1_hemi{15}.tSeries];
rois_AmyV1_hemi_combinedROI2{9}.tSeries =[rois_AmyV1_hemi{16}.tSeries;rois_AmyV1_hemi{17}.tSeries];
rois_AmyV1_hemi_combinedROI2{10}.tSeries=[rois_AmyV1_hemi{18}.tSeries;rois_AmyV1_hemi{19}.tSeries];
rois_AmyV1_hemi_combinedROI2{11}.tSeries=[rois_AmyV1_hemi{20}.tSeries;rois_AmyV1_hemi{21}.tSeries];
rois_AmyV1_hemi_combinedROI2{12}.tSeries=[rois_AmyV1_hemi{22}.tSeries;rois_AmyV1_hemi{23}.tSeries];
rois_AmyV1_hemi_combinedROI2{13}.tSeries=[rois_AmyV1_hemi{24}.tSeries;rois_AmyV1_hemi{25}.tSeries];



%% 13 ROIs
for i=1:length(rois_AmyV1_hemi_combinedROI2)
    if ~isempty(rois_AmyV1_hemi_combinedROI2{i}.tSeries())& ~isempty(roi2_threshold_combinedROI2{i}')
rois_AmyV1_hemi_combinedROI2_tS{i} = rois_AmyV1_hemi_combinedROI2{i}.tSeries();
rois_AmyV1_hemi_combinedROI2_tS_tr{i}=rois_AmyV1_hemi_combinedROI2_tS{i}(ismember(rois_AmyV1_hemi_combinedROI2{i}.scanCoords',roi2_threshold_combinedROI2{i}','rows'),:);
resps_combinedROI2{i} = getr2timecourse(mean(rois_AmyV1_hemi_combinedROI2_tS_tr{i}), d2.nhdr, 1, d2.scm);
    end
end

for i=1:length(rois_AmyV1_hemi_combinedROI2)
%model{i} = resps{i}.scm .* resps{i}.ehdr';
 if ~isempty(resps_combinedROI2{i})
data_13ROI(i,:)= resps_combinedROI2{i}.ehdr';
errlow_13ROI(i,:) = resps_combinedROI2{i}.ehdrste';
errhigh_13ROi(i,:) = resps_combinedROI2{i}.ehdrste';
end
end
