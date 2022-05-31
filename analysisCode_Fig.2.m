clear all
%% open mrloadRet

%% under 'GLM_all' in the faceLoc concatenation 
rois = loadROITSeries(getMLRView, []); % load all 50 lh+rh.Wang2015atlas ROIs
rois_Amy_faceLoc = loadROITSeries(getMLRView, 'amy');
rois_FFA_faceLoc = loadROITSeries(getMLRView, 'ffa');
rois_V1_parafoveal = loadROITSeries(getMLRView, 'V1_parafoveal');
rois_V1_peripheral = loadROITSeries(getMLRView, 'V1_peripheral');

% 50 from wang atlas + amy + FFA + V1 parafoveal + V1 peripheral= 54 rois
rois{51} = rois_Amy_faceLoc; rois{52} = rois_FFA_faceLoc; rois{53}= rois_V1_parafoveal; rois{54} = rois_V1_peripheral;

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

% Same roi label combined across hemispheres, wang atlas 50 roi labels
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
% Amy 
roi2_threshold_combinedROI2{14}=roi2_threshold{51};
% FFA 
roi2_threshold_combinedROI2{15}=roi2_threshold{52};
% V1_parafoveal
roi2_threshold_combinedROI2{16}=roi2_threshold{53};
% V1_peripheral 
roi2_threshold_combinedROI2{17}=roi2_threshold{54};
% roi2_threshold_combinedROI2 (1x17 rois)

%% under 'event-related analysis' in the amyV1 task
d2 = viewGet(getMLRView, 'd');

rois_AmyV1_FC = loadROITSeries(getMLRView, []); % load all 50 lh and rh.Wang2015atlas ROIs (same order!)
rois_Amy = loadROITSeries(getMLRView, 'amy');
rois_FFA = loadROITSeries(getMLRView, 'ffa');% 50 from wang atlas + amy + FFA = 52 rois
rois_V1_para = loadROITSeries(getMLRView, 'V1_parafoveal');
rois_V1_peri = loadROITSeries(getMLRView, 'V1_peripheral');

%% for hemi combined
for j=1:25
   rois_AmyV1_FC_hemi{j}.tSeries=[rois_AmyV1_FC{j}.tSeries; rois_AmyV1_FC{j+25}.tSeries];
   rois_AmyV1_FC_hemi{j}.scanCoords=[rois_AmyV1_FC{j}.scanCoords, rois_AmyV1_FC{j+25}.scanCoords];
end

%% 15 ROIs
for k=1:3
rois_AmyV1_FC_hemi_combinedROI2{k}.tSeries = rois_AmyV1_FC_hemi{k}.tSeries;
rois_AmyV1_FC_hemi_combinedROI2{k}.scanCoords=rois_AmyV1_FC_hemi{k}.scanCoords;
end

rois_AmyV1_FC_hemi_combinedROI2{4}.scanCoords =[rois_AmyV1_FC_hemi{4}.scanCoords,rois_AmyV1_FC_hemi{5}.scanCoords,rois_AmyV1_FC_hemi{6}.scanCoords,rois_AmyV1_FC_hemi{7}.scanCoords,rois_AmyV1_FC_hemi{8}.scanCoords];
rois_AmyV1_FC_hemi_combinedROI2{5}.scanCoords =[rois_AmyV1_FC_hemi{9}.scanCoords,rois_AmyV1_FC_hemi{10}.scanCoords];
rois_AmyV1_FC_hemi_combinedROI2{6}.scanCoords =[rois_AmyV1_FC_hemi{11}.scanCoords,rois_AmyV1_FC_hemi{12}.scanCoords];
rois_AmyV1_FC_hemi_combinedROI2{7}.scanCoords =[rois_AmyV1_FC_hemi{13}.scanCoords];
rois_AmyV1_FC_hemi_combinedROI2{8}.scanCoords =[rois_AmyV1_FC_hemi{14}.scanCoords,rois_AmyV1_FC_hemi{15}.scanCoords];
rois_AmyV1_FC_hemi_combinedROI2{9}.scanCoords =[rois_AmyV1_FC_hemi{16}.scanCoords,rois_AmyV1_FC_hemi{17}.scanCoords];
rois_AmyV1_FC_hemi_combinedROI2{10}.scanCoords=[rois_AmyV1_FC_hemi{18}.scanCoords,rois_AmyV1_FC_hemi{19}.scanCoords];
rois_AmyV1_FC_hemi_combinedROI2{11}.scanCoords=[rois_AmyV1_FC_hemi{20}.scanCoords,rois_AmyV1_FC_hemi{21}.scanCoords];
rois_AmyV1_FC_hemi_combinedROI2{12}.scanCoords=[rois_AmyV1_FC_hemi{22}.scanCoords,rois_AmyV1_FC_hemi{23}.scanCoords];
rois_AmyV1_FC_hemi_combinedROI2{13}.scanCoords=[rois_AmyV1_FC_hemi{24}.scanCoords,rois_AmyV1_FC_hemi{25}.scanCoords];
rois_AmyV1_FC_hemi_combinedROI2{14}.scanCoords=[rois_Amy.scanCoords];
rois_AmyV1_FC_hemi_combinedROI2{15}.scanCoords=[rois_FFA.scanCoords];
rois_AmyV1_FC_hemi_combinedROI2{16}.scanCoords=[rois_V1_para.scanCoords];
rois_AmyV1_FC_hemi_combinedROI2{17}.scanCoords=[rois_V1_peri.scanCoords];

rois_AmyV1_FC_hemi_combinedROI2{4}.tSeries =[rois_AmyV1_FC_hemi{4}.tSeries;rois_AmyV1_FC_hemi{5}.tSeries;rois_AmyV1_FC_hemi{6}.tSeries;rois_AmyV1_FC_hemi{7}.tSeries;rois_AmyV1_FC_hemi{8}.tSeries];
rois_AmyV1_FC_hemi_combinedROI2{5}.tSeries =[rois_AmyV1_FC_hemi{9}.tSeries;rois_AmyV1_FC_hemi{10}.tSeries];
rois_AmyV1_FC_hemi_combinedROI2{6}.tSeries =[rois_AmyV1_FC_hemi{11}.tSeries;rois_AmyV1_FC_hemi{12}.tSeries];
rois_AmyV1_FC_hemi_combinedROI2{7}.tSeries =[rois_AmyV1_FC_hemi{13}.tSeries];
rois_AmyV1_FC_hemi_combinedROI2{8}.tSeries =[rois_AmyV1_FC_hemi{14}.tSeries;rois_AmyV1_FC_hemi{15}.tSeries];
rois_AmyV1_FC_hemi_combinedROI2{9}.tSeries =[rois_AmyV1_FC_hemi{16}.tSeries;rois_AmyV1_FC_hemi{17}.tSeries];
rois_AmyV1_FC_hemi_combinedROI2{10}.tSeries=[rois_AmyV1_FC_hemi{18}.tSeries;rois_AmyV1_FC_hemi{19}.tSeries];
rois_AmyV1_FC_hemi_combinedROI2{11}.tSeries=[rois_AmyV1_FC_hemi{20}.tSeries;rois_AmyV1_FC_hemi{21}.tSeries];
rois_AmyV1_FC_hemi_combinedROI2{12}.tSeries=[rois_AmyV1_FC_hemi{22}.tSeries;rois_AmyV1_FC_hemi{23}.tSeries];
rois_AmyV1_FC_hemi_combinedROI2{13}.tSeries=[rois_AmyV1_FC_hemi{24}.tSeries;rois_AmyV1_FC_hemi{25}.tSeries];
rois_AmyV1_FC_hemi_combinedROI2{14}.tSeries=[rois_Amy.tSeries];
rois_AmyV1_FC_hemi_combinedROI2{15}.tSeries=[rois_FFA.tSeries];
rois_AmyV1_FC_hemi_combinedROI2{16}.tSeries=[rois_V1_para.tSeries];
rois_AmyV1_FC_hemi_combinedROI2{17}.tSeries=[rois_V1_peri.tSeries];

%% run the deconvolution on the mean time series
nhdr = 3; % 3 predictors
hdrlen = 15; % 3T BOLD: 15 TR x 2s = 30s
%hdrlen = 18; % 7T BOLD: 18 TR x 1.5s = 27s

%% take the average tSeries across all voxels in an ROI
for i=1:length(rois_AmyV1_FC_hemi_combinedROI2) 
    if ~isempty(rois_AmyV1_FC_hemi_combinedROI2{i}.tSeries()) && ~isempty(roi2_threshold_combinedROI2{i}')
        rois_AmyV1_FC_hemi_combinedROI2_tSeries{i} = rois_AmyV1_FC_hemi_combinedROI2{i}.tSeries();
        rois_AmyV1_FC_hemi_combinedROI2_tSeries_threshold{i}=rois_AmyV1_FC_hemi_combinedROI2_tSeries{i}(ismember(rois_AmyV1_FC_hemi_combinedROI2{i}.scanCoords',roi2_threshold_combinedROI2{i}','rows'),:);
        tSeries{i} = mean(rois_AmyV1_FC_hemi_combinedROI2_tSeries_threshold{i});
        resps_combinedROI2{i} = getr2timecourse(tSeries{i}, nhdr, hdrlen, d2.scm);
    end
end

%% recover the model from the resulting d structure
% for 3T
clear model
for i = 1:length(rois_AmyV1_FC_hemi_combinedROI2)
    if ~isempty(resps_combinedROI2{i})
        model{i}(1,:) = conv(d2.scm(:,1), resps_combinedROI2{i}.ehdr(1,:));
        model{i}(2,:) = conv(d2.scm(:,16), resps_combinedROI2{i}.ehdr(2,:));
        model{i}(3,:) = conv(d2.scm(:,31), resps_combinedROI2{i}.ehdr(3,:));
        model{i} = model{i}(:,1:length(tSeries{15}));
    end
end
% % for 7T
% clear model
% for i = 1:length(rois_AmyV1_FC_hemi_combinedROI2)
%     if ~isempty(resps_combinedROI2{i})
%         model{i}(1,:) = conv(d2.scm(:,1), resps_combinedROI2{i}.ehdr(1,:));
%         model{i}(2,:) = conv(d2.scm(:,19), resps_combinedROI2{i}.ehdr(2,:));
%         model{i}(3,:) = conv(d2.scm(:,37), resps_combinedROI2{i}.ehdr(3,:));
%         model{i} = model{i}(:,1:length(tSeries{15}));
%     end
% end
% 
%% plot the data and best fitting model and confirm that they look good.
FigHandle = figure('Position', [100, 100, 1800, 1500]);
for i = 1:length(rois_AmyV1_FC_hemi_combinedROI2)
    pctTseries{i} = 100*(tSeries{i}-1);
    subplot(4,5,i)
    plot(pctTseries{i}, 'k');
    hold on
    plot(sum(model{i}), 'r');
end

% now we are going to subtract them to get the residuals
FigHandle = figure('Position', [100, 100, 1800, 1500]);
for i = 1:length(rois_AmyV1_FC_hemi_combinedROI2)
    residual{i} = pctTseries{i}-sum(model{i});
    subplot(4,5,i)
    plot(residual{i}, 'b');
end


%% check V1 and AMG: measured timeseries and model fits 
FigHandle = figure('Position', [100, 100, 1600, 800]);
pctTseries{9} = 100*(tSeries{9}-1);
subplot(2,1,1) % V1
plot(pctTseries{9}, 'color',[0 .8 0],'LineWidth',2); % green
hold on
plot(sum(model{9}),'color', [1 174/255 66/255],'LineWidth',2); %yellow
hold on 
plot(residual{9}, 'color',[168/255 4/255 225/255],'LineWidth',2); %violet
xlim([0,168])
ylim([-3 3]) % for para
title('V1','FontSize',20,'Color','k')
legend('measured time series','mean stimulus response','residual time series')

subplot(2,1,2) % Amygdala
plot(pctTseries{14}, 'color',[0 .8 0],'LineWidth',2); % green
hold on
plot(sum(model{14}),'color', [1 174/255 66/255],'LineWidth',2); %yellow
hold on 
plot(residual{14}, 'color',[168/255 4/255 225/255],'LineWidth',2); %violet
xlim([0,168])
ylim([-1.5 1.5]) % for para
title('Amygdala','FontSize',20,'Color','k')
legend('measured time series','mean stimulus response','residual time series')

%% run correlation
clear r p 
for i =1:length(rois_AmyV1_FC_hemi_combinedROI2)
    for j=1:length(rois_AmyV1_FC_hemi_combinedROI2)
        if  isempty(resps_combinedROI2{i})|| isempty(resps_combinedROI2{j})
            r(i,j)=0; p(i,j)=1;
          else     [r(i,j),p(i,j)]=corr(residual{i}',residual{j}');
        end
    end
end

%% now we are going to seperate the residuals by condition 
stimTiming=d2.stimvol;
  
for i = 1:length(stimTiming{1})
    for k = 1:nhdr
     stimDuration(i,:,k)=[stimTiming{k}(i):stimTiming{k}(i)+hdrlen-1]; % 1-15 TR=15
    end
end

stimDuration_fear = reshape(stimDuration(:,:,1)',1,[]);
stimDuration_happy = reshape(stimDuration(:,:,2)',1,[]);
stimDuration_neutral = reshape(stimDuration(:,:,3)',1,[]);

for i = 1:length(rois_AmyV1_FC_hemi_combinedROI2)
    if ~isempty(tSeries{i})
        residual_fear{i} = residual{i}(stimDuration_fear);
        residual_happy{i} = residual{i}(stimDuration_happy);
        residual_neutral{i} = residual{i}(stimDuration_neutral);
        pctTseries_fear{i} = pctTseries{i}(stimDuration_fear);
        pctTseries_happy{i} = pctTseries{i}(stimDuration_happy);
        pctTseries_neutral{i} = pctTseries{i}(stimDuration_neutral);
        model_sum = sum(model{i});
        model_fear{i} = model_sum(stimDuration_fear);
        model_happy{i} = model_sum(stimDuration_happy);
        model_neutral{i} = model_sum(stimDuration_neutral);
    end
end
 
%% check
length(residual_fear{2}) % 315 for 3T           % 432 for 7T

%% run correlation on negative vs. positive valence
for i =1:length(rois_AmyV1_FC_hemi_combinedROI2)
    for j=1:length(rois_AmyV1_FC_hemi_combinedROI2)
        if  isempty(resps_combinedROI2{i})|| isempty(resps_combinedROI2{j})
            r_fear(i,j)=0; r_happy(i,j)=0;r_neutral(i,j)=0;
            p_fear(i,j)=1; p_happy(i,j)=1;p_neutral(i,j)=1;
        else
            [r_fear(i,j),p_fear(i,j)]=corr(residual_fear{i}',residual_fear{j}');
            [r_happy(i,j),p_happy(i,j)]=corr(residual_happy{i}',residual_happy{j}');
            [r_neutral(i,j),p_neutral(i,j)]=corr(residual_neutral{i}',residual_neutral{j}');
        end
    end
end

%% get r for each scan.participant

%% group stats 
% file_path='/Users/liut7/OneDrive - National Institutes of Health/AmyV1_analysis/FC_individual'
file_path='https://github.com/tinaliutong/layerfmri_AMG_V1/tree/main/SourceData/Fig.2_FC_individual'  

% download this folder and go to this folder
files = dir('*.txt');
N = length(files);
% loop for each file 
for i = 1:N
    thisfilename = files(i).name ;
    thisdata = load(thisfilename); %load just this file
    fulldata(:,:,i) = thisdata; % dim: 14*42*15
    %fprintf( 'File #%d, "%s", max and min value was: %g\n', i, thisfilename, max(thisdata(:)), min(thisdata(:)));
    negaVal(:,:,i) = fulldata(:,1:14,i) - fulldata(:,29:42,i);
    posiVal(:,:,i) = fulldata(:,15:28,i) - fulldata(:,29:42,i);
end
meanFearful = nanmean(fulldata(:,1:14,:),3);
meanHappy = nanmean(fulldata(:,15:28,:),3);
meanNeutral = nanmean(fulldata(:,29:42,:),3);
meanNegaVal = nanmean(negaVal(:,:,:),3);
meanPosiVal = nanmean(posiVal(:,:,:),3);

%% amy
for i=1:15 % subject
    for j = 1:14 % region
        AMGConn_fearful(i,j) = fulldata(j,1,i); 
        AMGConn_happy(i,j) = fulldata(j,15,i); 
        AMGConn_neutral(i,j) = fulldata(j,29,i); 
        % column: residual correlation with each of the 15 regions 
        % row: each subject
    end
end
% V1c
for i=1:15 % subject
    for j = 1:14 % region
        V1cConn_fearful(i,j) = fulldata(j,2,i); 
        V1cConn_happy(i,j) = fulldata(j,16,i); 
        V1cConn_neutral(i,j) = fulldata(j,30,i); 
    end
end
% V1p
for i=1:15 % subject
    for j = 1:14 % region
        V1pConn_fearful(i,j) = fulldata(j,3,i); 
        V1pConn_happy(i,j) = fulldata(j,17,i); 
        V1pConn_neutral(i,j) = fulldata(j,31,i); 
    end
end
% V2
for i=1:15 % subject
    for j = 1:14 % region
        V2Conn_fearful(i,j) = fulldata(j,4,i); 
        V2Conn_happy(i,j) = fulldata(j,18,i); 
        V2Conn_neutral(i,j) = fulldata(j,32,i); 
    end
end
% V3 column, 5th row in figure
for i=1:15 % subject
    for j = 1:14 % region
        V3Conn_fearful(i,j) = fulldata(j,5,i); 
        V3Conn_happy(i,j) = fulldata(j,19,i); 
        V3Conn_neutral(i,j) = fulldata(j,33,i); 
    end
end
% hV4 column, 6th row in figure
for i=1:15 % subject
    for j = 1:14 % region
        hV4Conn_fearful(i,j) = fulldata(j,6,i); 
        hV4Conn_happy(i,j) = fulldata(j,20,i); 
        hV4Conn_neutral(i,j) = fulldata(j,34,i); 
    end
end
%% FFA column, 7th row in figure
for i=1:15 % subject
    for j = 1:14 % region
        FFAConn_fearful(i,j) = fulldata(j,7,i); 
        FFAConn_happy(i,j) = fulldata(j,21,i); 
        FFAConn_neutral(i,j) = fulldata(j,35,i); 
    end
end
%% VO column, 8th row in figure
for i=1:15 % subject
    for j = 1:14 % region
        VOConn_fearful(i,j) = fulldata(j,8,i); 
        VOConn_happy(i,j) = fulldata(j,22,i); 
        VOConn_neutral(i,j) = fulldata(j,36,i); 
    end
end
%% PHC column, 10th row in figure
for i=1:15 % subject
    for j = 1:14 % region
        PHCConn_fearful(i,j) = fulldata(j,9,i); 
        PHCConn_happy(i,j) = fulldata(j,23,i); 
        PHCConn_neutral(i,j) = fulldata(j,37,i); 
    end
end
%% LO column, 10th row in figure
for i=1:15 % subject
    for j = 1:14 % region
        LOConn_fearful(i,j) = fulldata(j,10,i); 
        LOConn_happy(i,j) = fulldata(j,24,i); 
        LOConn_neutral(i,j) = fulldata(j,38,i); 
    end
end
%% TO column, 11th row in figure
for i=1:15 % subject
    for j = 1:14 % region
        TOConn_fearful(i,j) = fulldata(j,11,i); 
        TOConn_happy(i,j) = fulldata(j,25,i); 
        TOConn_neutral(i,j) = fulldata(j,39,i); 
    end
end
%% V3ab
for i=1:15 % subject
    for j = 1:14 % region
        V3ABConn_fearful(i,j) = fulldata(j,12,i); 
        V3ABConn_happy(i,j) = fulldata(j,26,i); 
        V3ABConn_neutral(i,j) = fulldata(j,40,i); 
        % column: residual correlation with each of the 15 regions 
        % row: each subject
    end
end
%% IPS0 column, 13st row in figure
for i=1:15 % subject
    for j = 1:14 % region
        IPS0Conn_fearful(i,j) = fulldata(j,13,i); 
        IPS0Conn_happy(i,j) = fulldata(j,27,i); 
        IPS0Conn_neutral(i,j) = fulldata(j,41,i); 
        % column: residual correlation with each of the 15 regions 
        % row: each subject
    end
end
%% IPS1-5 column, 13st row in figure
for i=1:15 % subject
    for j = 1:14 % region
        IPS15Conn_fearful(i,j) = fulldata(j,14,i); 
        IPS15Conn_happy(i,j) = fulldata(j,28,i); 
        IPS15Conn_neutral(i,j) = fulldata(j,42,i); 
        % column: residual correlation with each of the 15 regions 
        % row: each subject
    end
end
%% amygdala-central V1 column, 1st row in figure
k=2
[p_AmgV1c_Fearful,h_AmgV1c_Fearful,stats_AmgV1c_Fearful] = signrank(AMGConn_fearful(:,k))
[h_AmgV1c_F,p_AmgV1c_F,ci_AmgV1c_F,stats_AmgV1c_F] = ttest(AMGConn_fearful(:,k))
[p_AmgV1c_Happy,h_AmgV1c_Happy,stats_AmgV1c_happy] = signrank(AMGConn_happy(:,k))
[h_AmgV1c_H,p_AmgV1c_H,ci_AmgV1c_H,stats_AmgV1c_H] = ttest(AMGConn_happy(:,k))
[p_AmgV1c_Neutral,h_AmgV1c_Neutral,stats_AmgV1c_Neutral] = signrank(AMGConn_neutral(:,k))
[h_AmgV1c_N,p_AmgV1c_N,ci_AmgV1c_N,stats_AmgV1c_N] = ttest(AMGConn_neutral(:,k))

% paired-sampes t-test
k=2
[h_AmgV1c_negaVal,p_AmgV1c_negaVal,ci_AmgV1c_negaVal,stats_AmgV1c_negaVal] = ttest(AMGConn_fearful(:,k),AMGConn_neutral(:,k))
[h_AmgV1c_posiVal,p_AmgV1c_posiVal,ci_AmgV1c_posiVal,stats_AmgV1c_posiVal] = ttest(AMGConn_happy(:,k),AMGConn_neutral(:,k))
p_AmgV1c_negaVal
p_AmgV1c_negaVal*14

k=3 %AMy-V1_peri
[p_AmgV1p_Fearful,h_AmgV1p_Fearful,stats_AmgV1p_Fearful] = signrank(AMGConn_fearful(:,k))
[h_AmgV1p_F,p_AmgV1p_F,ci_AmgV1p_F,stats_AmgV1p_F] = ttest(AMGConn_fearful(:,k))
[p_AmgV1p_Happy,h_AmgV1p_Happy,stats_AmgV1p_happy] = signrank(AMGConn_happy(:,k))
[h_AmgV1p_H,p_AmgV1p_H,ci_AmgV1p_H,stats_AmgV1p_H] = ttest(AMGConn_happy(:,k))
[p_AmgV1p_Neutral,h_AmgV1p_Neutral,stats_AmgV1p_Neutral] = signrank(AMGConn_neutral(:,k))
[h_AmgV1p_N,p_AmgV1p_N,ci_AmgV1p_N,stats_AmgV1p_N] = ttest(AMGConn_neutral(:,k))

% paired-sampes t-test
k=3 %AMy-V1_peri
[h_AmgV1p_negaVal,p_AmgV1p_negaVal,ci_AmgV1p_negaVal,stats_AmgV1p_negaVal] = ttest(AMGConn_fearful(:,k),AMGConn_neutral(:,k))
[h_AmgV1p_posiVal,p_AmgV1p_posiVal,ci_AmgV1p_posiVal,stats_AmgV1p_posiVal] = ttest(AMGConn_happy(:,k),AMGConn_neutral(:,k))
p_AmgV1p_negaVal*14
 


U=tril(meanFearful)
 
% making 0s as NaNs
for i=1:14
    for j = 1:14
        if U(j,i)== 0
            U(j,i)= nan;
        end      
    end
end

FigHandle = figure('Position', [100, 100, 1600, 600]);
subplot(1,2,1)
imagesc(U)
 [nr,nc] = size(U)
pcolor([U nan(nr,1); nan(1,nc+1)]);
 shading flat;
set(gca, 'ydir', 'reverse');
 
max(U(:))
min(U(:))
set(gca,'TickLength',[0 0])
set(gca,'xaxisLocation','bottom')

caxis([0 1])
% caxis([0 1])
xticks([1.5:1.02:14.8])
yticks([1.28:1:14.28])

yticklabels({'AMG','central V1','peripheral V1','V2','V3','hV4','FFA','VO','PHC','LO','TO','V3AB','IPS0','IPS1-5'});
xticklabels({'AMG','central V1','peripheral V1','V2','V3','hV4','FFA','VO','PHC','LO','TO','V3AB','IPS0','IPS1-5'});
colormap(pink)

xtickangle(45)
ytickangle(45)

ax = gca;
ax.FontSize = 16;
axis square
box off
colorbar

U=tril(meanHappy)
 
% making 0s as NaNs
for i=1:14
    for j = 1:14
        if U(j,i)== 0
            U(j,i)= nan;
        end      
    end
end

subplot(1,2,1)
imagesc(U)
 [nr,nc] = size(U)
pcolor([U nan(nr,1); nan(1,nc+1)]);
 shading flat;
set(gca, 'ydir', 'reverse');
 
max(U(:))
min(U(:))
set(gca,'TickLength',[0 0])
set(gca,'xaxisLocation','bottom')

caxis([0 1])
% caxis([0 1])
xticks([1.5:1.02:14.8])
yticks([1.28:1:14.28])

yticklabels({'AMG','central V1','peripheral V1','V2','V3','hV4','FFA','VO','PHC','LO','TO','V3AB','IPS0','IPS1-5'});
xticklabels({'AMG','central V1','peripheral V1','V2','V3','hV4','FFA','VO','PHC','LO','TO','V3AB','IPS0','IPS1-5'});
colormap(pink)

xtickangle(45)
ytickangle(45)

ax = gca;
ax.FontSize = 16;
axis square
box off
colorbar

 U=tril(meanNeutral)

% making 0s as NaNs
for i=1:14
    for j = 1:14
        if U(j,i)== 0
            U(j,i)= nan;
        end      
    end
end

subplot(1,2,2)
imagesc(U)
 [nr,nc] = size(U)
pcolor([U nan(nr,1); nan(1,nc+1)]);
 shading flat;
set(gca, 'ydir', 'reverse');
 
max(U(:))
min(U(:))
set(gca,'TickLength',[0 0])
set(gca,'xaxisLocation','bottom')

caxis([0 1])
% caxis([0 1])
xticks([1.5:1.02:14.8])
yticks([1.28:1:14.28])

yticklabels({'AMG','central V1','peripheral V1','V2','V3','hV4','FFA','VO','PHC','LO','TO','V3AB','IPS0','IPS1-5'});
xticklabels({'AMG','central V1','peripheral V1','V2','V3','hV4','FFA','VO','PHC','LO','TO','V3AB','IPS0','IPS1-5'});
colormap(pink)

xtickangle(45)
ytickangle(45)

ax = gca;
ax.FontSize = 16;
axis square
box off
colorbar


file_path='/Users/liut7/OneDrive - National Institutes of Health/AmyV1_analysis/FC_individual'
 % go to this folder
files = dir('*.txt'); % you are in the folder of files 
N = length(files);
% loop for each file 
for i = 1:N
    thisfilename = files(i).name ;
    thisdata = load(thisfilename); %load just this file
    fulldata(:,:,i) = thisdata; % dim: 14*42*15
    %fprintf( 'File #%d, "%s", max and min value was: %g\n', i, thisfilename, max(thisdata(:)), min(thisdata(:)));
    negaVal(:,:,i) = fulldata(:,1:14,i) - fulldata(:,29:42,i);
    posiVal(:,:,i) = fulldata(:,15:28,i) - fulldata(:,29:42,i);
end
meanFearful = nanmean(fulldata(:,1:14,:),3)
meanHappy = nanmean(fulldata(:,15:28,:),3)
meanNeutral = nanmean(fulldata(:,29:42,:),3)
meanNegaVal = nanmean(negaVal(:,:,:),3) 
meanPosiVal = nanmean(posiVal(:,:,:),3) 

%% plot valence effect: negative and positive
U=tril(meanNegaVal)
 
% making 0s as NaNs
for i=1:14
    for j = 1:14
        if U(j,i)== 0
            U(j,i)= nan;
        end      
    end
end

FigHandle = figure('Position', [100, 100, 1600, 600]);
subplot(1,2,1)
imagesc(U)
 [nr,nc] = size(U)
pcolor([U nan(nr,1); nan(1,nc+1)]);
 shading flat;
set(gca, 'ydir', 'reverse');
 
max(U(:))
min(U(:))
set(gca,'TickLength',[0 0])
set(gca,'xaxisLocation','bottom')

caxis([-0.01 .07])
% caxis([0 1])
xticks([1.5:1.02:14.8])
yticks([1.28:1:14.28])

yticklabels({'AMG','central V1','peripheral V1','V2','V3','hV4','FFA','VO','PHC','LO','TO','V3AB','IPS0','IPS1-5'});
xticklabels({'AMG','central V1','peripheral V1','V2','V3','hV4','FFA','VO','PHC','LO','TO','V3AB','IPS0','IPS1-5'});
colormap(pink)

xtickangle(45)
ytickangle(45)

ax = gca;
ax.FontSize = 16;
axis square
box off
colormap('hot');
colorbar

U=tril(meanPosiVal)
 
% making 0s as NaNs
for i=1:14
    for j = 1:14
        if U(j,i)== 0
            U(j,i)= nan;
        end      
    end
end

subplot(1,2,2)
imagesc(U)
 [nr,nc] = size(U)
pcolor([U nan(nr,1); nan(1,nc+1)]);
 shading flat;
set(gca, 'ydir', 'reverse');
 
max(U(:))
min(U(:))
set(gca,'TickLength',[0 0])
set(gca,'xaxisLocation','bottom')

caxis([-0.01 .07])
% caxis([0 1])
xticks([1.5:1.02:14.8])
yticks([1.28:1:14.28])

yticklabels({'AMG','central V1','peripheral V1','V2','V3','hV4','FFA','VO','PHC','LO','TO','V3AB','IPS0','IPS1-5'});
xticklabels({'AMG','central V1','peripheral V1','V2','V3','hV4','FFA','VO','PHC','LO','TO','V3AB','IPS0','IPS1-5'});
colormap(pink)

xtickangle(45)
ytickangle(45)

ax = gca;
ax.FontSize = 16;
axis square
box off
colormap('hot');
colorbar

