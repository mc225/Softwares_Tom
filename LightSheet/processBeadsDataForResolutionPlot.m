% processBeadsDataForResolutionPlot(xRange,yRange,zRange,localizationData,testData)
% 
%
function processBeadsDataForResolutionPlot(inputFolder)
    close all;
    if (nargin<1)
%         inputFolder='Z:\RESULTS\20120501_alpha7_int1000_g1_littlemorepower_betterfocus_beads_leftofprobe\backwardScans\';
%         inputFolder='Z:\RESULTS\20120501_alpha7_int1000_g1_littlemorepower_betterfocus_beads_leftofprobe\forwardScans\';
%         inputFolder='Z:\RESULTS\20120501_alpha7_int1000_g1_littlemorepower_betterfocus_beads_leftofprobe\forwardScans2\';
        inputFolder='Z:\RESULTS\20120501_alpha7_int1000_g1_morepower_beads_leftofprobe\';
%         load('Z:\RESULTS\20120501_alpha7_int1000_g1_littlemorepower_betterfocus_beads_leftofprobe\forwardScans\Gaussian_0F.mat');
%         load('Z:\RESULTS\20120501_alpha7_int1000_g1_littlemorepower_betterfocus_beads_leftofprobe\forwardScans\Airy_0F.mat');
%         load('Z:\RESULTS\600nmBeadsApOneThird\B\at0\recording_lambda532nm_alpha7_beta100.mat');
%         load('E:\Stored Files\RESULTS\compressed\offsetBeads\at_25um\recording_lambda532nm_alpha0_beta100.mat');
%            load('E:\Stored Files\RESULTS\2013_01_15_SmallBeadSample\2013-01-15 16_14_11.433_green_0.05_(2)\recording_lambda532nm_alpha0_beta100.mat');
%         load('Z:\RESULTS\2013-04-10_Beads200nm\2013-04-10 16_10_20.097_ApOneHalf\recording_lambda488nm_alpha10_beta100.mat');
%         load('Z:\RESULTS\2013-04-10_Beads200nm\2013-04-10 16_10_20.097_ApOneHalf\recording_lambda488nm_alpha0_beta100.mat');
%         load('E:\Stored Files\RESULTS\compressed\2013-01-29 11_11_55.143_mixedSmallBeads_inAgarose_on_PDMS_green_blue\recording_lambda488nm_alpha5_beta100.mat');
    end
    
    scanType='Bessel5';
    scanDir='B';
    deconvolvedTestData=strcmpi(scanType,'Airy');
    
    [xRange,yRange,zRange,localizationData,testData]=loadFrom(inputFolder,['Airy_0',scanDir,'.mat'],[scanType,'_0',scanDir,'.mat'],deconvolvedTestData);
        
    outputPath=fullfile(inputFolder,['images',scanType,'.mat']);
    
    voxelSize=[diff(xRange(1:2)) diff(yRange(1:2)) diff(zRange(1:2))];
    
    if (deconvolvedTestData)
        shiftInXYPerZ=[0 0];
    else
        shiftInXYPerZ=[0.026 0.03];
    end
    
    threshold=localizationData(1:16:end,1:16:end,1:16:end);
    backgroundLevel=0; %median(threshold(:));
    threshold=quantile(threshold(:),.9999);
    
	hotSpots=localizationData>threshold;
    regions  = regionprops(hotSpots, 'Centroid','BoundingBox');
    boundingBoxes = cat(1, regions.BoundingBox);
    boundingBoxSizes=boundingBoxes(:,4:6).*repmat(voxelSize,[numel(regions) 1]);
    
    localizationBoxSize=[1.25 1.25 1.75]*1e-6;
    maxDevBoundingBoxSize=[1 1 1.5]*1e-6;
    projBoundingBoxSize=[5 5 30]*1e-6;
    displayBoundingBoxSize=[5 5 30]*1e-6;
    
%     % Co-register the testData
%     localizationSliceXZ=squeeze(mean(localizationData(:,abs(yRange)<5e-6,zRange<80e-6),2));
%     testSliceXZ=squeeze(mean(testData(:,abs(yRange)<5e-6,zRange<80e-6),2));
%     localizationSliceXY=squeeze(mean(localizationData(:,abs(yRange)<50e-6,zRange>40e-6 & zRange<80e-6),3));
%     testSliceXY=squeeze(mean(testData(:,abs(yRange)<50e-6,zRange>40e-6 & zRange<80e-6),3));
% %     [shiftXZ Greg] = dftregistration(fft2(localizationSliceXZ),fft2(testSliceXZ),100);
% %     dataShiftXZ=shiftXZ(3:4).*[diff(xRange(1:2)) diff(yRange(1:2))];
%     [shiftXY Greg] = dftregistration(fft2(localizationSliceXY),fft2(testSliceXY),100);
%     dataShiftXY=shiftXY(3:4).*[diff(xRange(1:2)) diff(yRange(1:2))];
%     figure();axs(1)=subplot(1,2,1);showImage(localizationSliceXZ,-1,zRange(zRange<80e-6)*1e6,xRange*1e6); axis equal; axs(2)=subplot(1,2,2);showImage(ifft2(Greg,'symmetric'),-1,zRange(zRange<80e-6)*1e6,xRange*1e6); axis equal; linkaxes(axs)
    
    % Filter out clusters and noise
    regions=regions(all(abs(boundingBoxSizes-repmat(localizationBoxSize,[numel(regions) 1])) < repmat(maxDevBoundingBoxSize,[numel(regions) 1]),2));
    centroids = cat(1, regions.Centroid);
    centroids=centroids(:,[2 1 3]);
    centroids(:,1)=xRange(round(centroids(:,1)));
    centroids(:,2)=yRange(round(centroids(:,2)));
    centroids(:,3)=zRange(round(centroids(:,3)));
    
    centroids=sortrows(centroids,2);
    
    fwhmsZ=zeros(1,numel(regions));
    fwhmsX=fwhmsZ;
    fwhmsY=fwhmsX;
%     shiftXYs=zeros(numel(regions),2);
    samples=[];
    for beadIdx=1:numel(regions)
        centroid=centroids(beadIdx,:);
        
        xSel=abs(xRange-centroid(1))<=projBoundingBoxSize(1)/2;
        ySel=abs(yRange-centroid(2))<=projBoundingBoxSize(2)/2;
        zSel=abs(zRange-centroid(3))<=projBoundingBoxSize(3)/2;
        localizationDataSubCube=localizationData(xSel,ySel,zSel);
        localizationDataSubCubeProjZ=max(localizationDataSubCube,[],3);
        % Translate the test data cube back to origin
        shiftInXY=round(shiftInXYPerZ*centroid(3)./voxelSize(1:2));
        testDataSubCube=testData(circshift(xSel,[0 shiftInXY(1)]),circshift(ySel,[0 shiftInXY(2)]),zSel);
        testDataSubCubeProjZ=max(testDataSubCube,[],3);
        if all(size(localizationDataSubCubeProjZ)==size(testDataSubCubeProjZ))
%             % Test alignment
%             [shiftXY Greg] = dftregistration(fft2(localizationDataSubCubeProjZ),fft2(testDataSubCubeProjZ),10);
%             shiftXYs(beadIdx,:)=shiftXY(3:4);
%             testDataSubCube=circshift(testDataSubCube,1-round(shiftXYs(beadIdx,:))+floor(1+localizationBoxSize(1:2).*voxelSize(1:2)/2));
            testDataProjectionZ=squeeze(mean(mean(testDataSubCube))).';
            fwhmsZ(beadIdx)=calcFullWidthAtHalfMaximum(zRange(zSel),testDataProjectionZ-backgroundLevel,'Linear');
        else
            fwhmsZ(beadIdx)=Inf;
        end
        
        % Find FWHM in lateral dimensions
        testDataProjectionX=squeeze(mean(mean(testDataSubCube,2),3)).';
        fwhmsX(beadIdx)=calcFullWidthAtHalfMaximum(xRange(xSel),testDataProjectionX-backgroundLevel,'Linear');
        testDataProjectionY=squeeze(mean(mean(testDataSubCube,1),3)).';
        fwhmsY(beadIdx)=calcFullWidthAtHalfMaximum(xRange(xSel),testDataProjectionY-backgroundLevel,'Linear');
        
        %
        % Store results
        %
        sample=struct();
        sample.centroid=centroid;
        sample.fwhm=[fwhmsX(beadIdx) fwhmsY(beadIdx) fwhmsZ(beadIdx)];
%         logMessage(sprintf('Propagation position: %0.1f um',sample.centroid(2)*1e6));
        xSelDisplay=abs(xRange-sample.centroid(1))<displayBoundingBoxSize(1)/2;
        ySelDisplay=abs(yRange-sample.centroid(2))<displayBoundingBoxSize(2)/2;
        zSelDisplay=abs(zRange-sample.centroid(3))<displayBoundingBoxSize(3)/2;
        sample.xRange=xRange(xSelDisplay);
        sample.yRange=yRange(ySelDisplay);
        sample.zRange=zRange(zSelDisplay);
        sample.localization=struct();
        sample.localization.proj1=squeeze(max(localizationData(xSelDisplay,ySelDisplay,zSelDisplay),[],1)).';
        sample.localization.proj2=squeeze(max(localizationData(xSelDisplay,ySelDisplay,zSelDisplay),[],2)).';
        sample.localization.proj3=squeeze(max(localizationData(xSelDisplay,ySelDisplay,zSelDisplay),[],3));
        sample.test=struct();
        % Translate the test data cube back to origin
        testDataSubCube=testData(circshift(xSelDisplay,[0 shiftInXY(1)]),circshift(ySelDisplay,[0 shiftInXY(2)]),zSelDisplay);
        sample.test.proj1=squeeze(max(testDataSubCube,[],1)).';
        sample.test.proj2=squeeze(max(testDataSubCube,[],2)).';
        sample.test.proj3=squeeze(max(testDataSubCube,[],3));
        
        if (~isempty(samples))
            samples(beadIdx)=sample;
        else
            samples=sample;
        end
    end
    save(outputPath,'samples');
    
    % Plot results
    scatterFig=figure();
    scatter(abs(centroids(:,2))*1e6,fwhmsZ*1e6,'+');
    xlim([0 150]);ylim([0 50]);
    title([scanType,' ',scanDir]);
    
end


function [xRange yRange,zRange,localizationData,testData]=loadFrom(folderName,localizationDataFileName,fittingDataFileName,loadProcessedData)
    if (nargin<3 || isempty(fittingDataFileName))
        fittingDataFileName=localizationDataFileName;
    end
    if (nargin<4)
        loadProcessedData=true;
    end
    load(fullfile(folderName,localizationDataFileName),'xRange');
    load(fullfile(folderName,localizationDataFileName),'yRange');
    load(fullfile(folderName,localizationDataFileName),'zRange');
    localizationData=load(fullfile(folderName,localizationDataFileName),'restoredDataCube');
    localizationData=localizationData.restoredDataCube;
    if loadProcessedData
        testData=load(fullfile(folderName,fittingDataFileName),'restoredDataCube');
        testData=testData.restoredDataCube;
    else
        testData=load(fullfile(folderName,fittingDataFileName),'recordedImageStack');
        testData=testData.recordedImageStack;
    end
end