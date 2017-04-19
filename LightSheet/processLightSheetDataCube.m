% [restoredDataCubeI xRange yRange zRange tRangeExtended recordedImageStack tRange ZOtf lightSheetPsf lightSheetOtf lightSheetDeconvFilter]
%   =processLightSheetDataCube(inputFileName,stageTranslationStepSize,excitation,detection,sample,detector,focusCoordinates,openFractionOfRadius,alpha,extraIntegrationTime,illuminationCroppingFactors,dataShear,deNoiseFrames,deconvolveWithDetectionPsf)
%
% Processes the datacube obtained after a light sheet capturing experiment.
% 
% Input:
%   inputFileName: The avi file containing the captured frames.
%   stageTranslationStepSize: z-step size of the sample in meters.
%   excitation: struct array with the fields:
%       wavelength
%       numericalApertureInAir
%       magnification
%       tubeLength
%       power
%   detection: struct array with the fields:
%       wavelength
%       numericalApertureInAir
%       magnification
%       tubeLength
%   sample: struct array with the fields:
%       fluorophore: struct array with the fields:
%           maxDensity
%           extinctionCoefficient: m^2
%           quantumYield
%       refractiveIndex: assuming the objective NA is in air
%       signalLevel: used for deconvolution
%       backgroundLevel
%   detector: struct array with the fields:
%       wavelength
%       numericalApertureInAir
%       magnification
%       tubeLength
%   focusCoordinates: The coordinates of the beam focus used to calculate the PSF for the
%          deconvolution. The coordinates are given in pixels and frames
%          (starting from 0)
%   openFractionOfRadius: The open fraction of the radius of the annular
%          aperture to create a Bessel beam, default 1 (full open aperture, top-hat beam).
%   alpha: The cubic phase mask value, default 0 (no cubic phase mask). The
%          optical path length is altered by alpha*u^3 wavelengths, where U is
%          the normalized pupil coordinate.
%   extraIntegrationTime: for Bessel beam
%   illuminationCroppingFactors: 
%   dataShear: The (fractional) pixel shift in X and Y of the projection when moving one step in Z.
%   deNoiseFrames: before deconvolution in Z, default false
%   deconvolveWithDetectionPsf: after deconvolution in Z, default false
%   maxFilesToLoad: The limit to the number of replicated scans that will be loaded if duplicates exist. (default: Inf)
%
% Output:
%   restoredDataCubeI:
%   xRange yRange zRange tRangeExtended:
%   recordedImageStack
%   tRange
%   ZOtf
%   lightSheetPsf
%   lightSheetOtf
%   lightSheetDeconvFilter
%
function [restoredDataCubeI xRange yRange zRange tRangeExtended recordedImageStack tRange ZOtf lightSheetPsf lightSheetOtf lightSheetDeconvFilter]=processLightSheetDataCube(inputFileName,stageTranslationStepSize,excitation,detection,sample,detector,focusCoordinateDecenter,openFractionOfRadius,alpha,extraIntegrationTime,illuminationCroppingFactors,dataShear,deNoiseFrames,deconvolveWithDetectionPsf,maxFilesToLoad)
    if (nargin<1)
        inputFileName='Z:\RESULTS\20120501_alpha7_int1000_g1_littlemorepower_betterfocus_beads_leftofprobe\Airy_0F.avi';
    end
    if (nargin<2)
    %     stageTranslationStepSize=0.00025e-3; % [m/frame]
        stageTranslationStepSize=100e-9; % [m/frame]
    end
    if (nargin<3 || isempty(excitation))
        % Objectives and wavelengths
        excitation={};
        excitation.wavelength=532e-9;
        excitation.objective={};
        excitation.objective.numericalAperture=.42;
        excitation.fractionOfNumericalApertureUsed=1;
        excitation.objective.magnification=20;
        excitation.objective.tubeLength=200e-3; %for Mitutoyo the rear focal length is 200mm
        excitation.power=0.01e-3;
        excitation.tubeLength=200e-3;
    end
    if (nargin<4 || isempty(detection))
        detection={};
        detection.wavelength=612e-9; % Firefli* Fluorescent Red (542/612nm) % 395e-9; for GFP
        detection.objective={};
        detection.objective.numericalAperture=.40; % .28;
        detection.objective.magnification=20; %Actual total magnification 22
        detection.objective.tubeLength=160e-3; % for Newport, %200e-3 for Mitutoyo / Nikon
        detection.tubeLength=1.1*160e-3; 
    end
    if (nargin<5 || isempty(sample))
        % Sample medium specification
        NAvogadro=6.0221417930e23;
        sample={};
        sample.fluorophore={};
        sample.fluorophore.maxDensity=0.01*10^3*NAvogadro; % m^-3, 0.01*mol/L density 
        sample.fluorophore.extinctionCoefficient=30000*(100/10^3)/NAvogadro; % m^2
        sample.fluorophore.quantumYield=0.79;
        sample.refractiveIndex=1.4; %Assuming the objective NA is in air
        sample.signalLevel=0.5; % for deconvolution
        sample.backgroundLevel=0.05; % fraction of dynamic range
    end
    if (nargin<6 || isempty(detector))
        % Detector specification
        detector={};
        detector.integrationTime=2e-3; % s
        detector.quantumEfficiency=0.55;
        detector.darkCurrent=200; % e-/s
        detector.readOutNoise=16; % e-
        detector.wellDepth=20e3; % e-
        detector.numberOfGrayLevels=2^12;
        detector.pixelSize=7.4*[1 1]*1e-6;
        detector.framesPerSecond=1; % [s]
        detector.center=[0 0]*1e-6;
    end
    if (nargin<7)
        focusCoordinateDecenter=[];
    end
    % Choice of light-sheet
    if (nargin<8)
        openFractionOfRadius=1.0; %Fraction of radius
    end
    if (nargin<9)
        alpha=0.0;
    end
    if (nargin<10)
        extraIntegrationTime=1.0; % Assume the bessel beam is illuminated equally
    end
    if (nargin<11 || isempty(illuminationCroppingFactors))
        illuminationCroppingFactors=[1 1; 1 1]*0.06;
    end
    if (nargin<12)
        dataShear=[];
    end
    if (nargin<13 || isempty(deNoiseFrames))
        deNoiseFrames=false;
    end
    if (nargin<14 || isempty(deconvolveWithDetectionPsf))
        deconvolveWithDetectionPsf=false;
    end
    if (nargin<15 || isempty(maxFilesToLoad))
        maxFilesToLoad=Inf;
    end
    
    % General constants
    hPlanck=6.6260695729e-34;
    cLight=299792458;
    
    %% Calculation parameters
    deflectBeamInsteadOfSampleMovement=false;
        
    if (alpha~=0)
        lightSheetName='Airy';
    else
        if (openFractionOfRadius==1.0)
            lightSheetName='Classic';
        else
            if (openFractionOfRadius==0.0)
                lightSheetName='Theoretical Bessel';
            else
                lightSheetName='Bessel';
            end
        end
    end
    
    %Assume that the power is adjusted to compensate for a reduction in open fraction.
    beamPowerAdjustment=1; %1-(1-openFractionOfRadius)^2;
    excitation.power=excitation.power*beamPowerAdjustment; % W
    
    %Load and pre-process data
    if (~strcmpi(inputFileName(end-4:end),'_.avi'))
    	recordedImageStack=readDataCubeFromAviFile(inputFileName);
    else
        pathName = fileparts(inputFileName);
        videoFiles=dir(strcat(inputFileName(1:end-4),'*.avi'));
        positionFiles=dir(strcat(inputFileName(1:end-4),'*.txt'));
        nbFilesToLoad=min(maxFilesToLoad,length(videoFiles));
        logMessage('Loading %u video files of %u ...',[nbFilesToLoad length(videoFiles)]);
        recordedImageStack=[];
        for (fileIdx=1:nbFilesToLoad)
            fileName=videoFiles(fileIdx).name;
            singlePositions=dlmread(strcat(pathName,'/',positionFiles(fileIdx).name),'\t');
            forward=mean(mean(diff(singlePositions)))>0;
            logMessage('Loading video file %s...',fileName);
            singleImageStack=readDataCubeFromAviFile(strcat(pathName,'/',fileName));
%             %Debug 
%             singleImageStack=getTestImage('usaf'); singleImageStack(1,1,10)=0;
%             if (~forward)
%                 singleImageStack=singleImageStack(:,:,[end end 1:end-2]);
%             end
            
            if (~forward)
                singleImageStack=singleImageStack(:,:,end:-1:1);
                singlePositions=singlePositions(end:-1:1,:);
            end
            save(strcat(pathName,'/',fileName(1:end-4),'.mat'),'singleImageStack');
            singleImageStackVectorFFT=fft(squeeze(max(max(singleImageStack,[],1),[],2)));
            if (~isempty(recordedImageStack))
                logMessage('Merging...');
                [output Greg] = dftregistration(recordedImageStackVectorFFT,singleImageStackVectorFFT,100);
                imageShift=output(3);
                logMessage('Detected a shift of %0.2f pixels',imageShift);
                singleImageStack=circshift(singleImageStack,[0 0 round(imageShift)]);
                recordedImageStack=recordedImageStack+singleImageStack;
            else
                recordedImageStack = singleImageStack;
                recordedImageStackVectorFFT=singleImageStackVectorFFT;
                positions=singlePositions(:,1);
            end
        end
%         recordedImageStack=ifftn(recordedImageStack,'symmetric');
    end
    clear singleImageStack;
    
    %Mask saturated beads
%     mask=recordedImageStack>=1;
%     logMessage('Masking %f saturated pixels.',sum(mask(:)));
%     mask=circshift(mask,-[2^4 2^4 2^7]);
%     for (xIdx=2.^[0:5])
%         mask=mask|circshift(mask,[xIdx 0 0]);
%     end
%     for (yIdx=2.^[0:5])
%         mask=mask|circshift(mask,[0 yIdx 0]);
%     end
%     for (zIdx=2.^[0:8])
%         mask=mask|circshift(mask,[0 0 zIdx]);
%     end

%     mask=0*recordedImageStack;
%     mask(:,:,1:710)=1;
%     recordedImageStack=recordedImageStack.*(1-mask);
%     logMessage('Masked %0.6f%% of the data.',100*mean(mask(:)));
%     clear mask;
    
    %Correct shear
    recordedImageStack=recordedImageStack/extraIntegrationTime;
    recordedImageStack=permute(recordedImageStack,[2 1 3]); % Dimension order [x y z]
    if (~isempty(dataShear))
        logMessage('Correcting a shear of (%0.3f,%0.3f)...',dataShear);
        recordedImageStack=correctShear(single(recordedImageStack),dataShear);
    end
    if (isempty(focusCoordinateDecenter))
        focusCoordinateDecenter=[0 0 0];
    end

    % Sample grid specification
    realDetectionMagnification=detection.objective.magnification*detection.objective.tubeLength/detection.tubeLength;
    xRange=detector.pixelSize(1)*([1:size(recordedImageStack,1)]-floor(size(recordedImageStack,1)/2)-1)/realDetectionMagnification; % left/right
    yRange=detector.pixelSize(2)*([1:size(recordedImageStack,2)]-floor(size(recordedImageStack,2)/2)-1)/realDetectionMagnification; % up/down
    tRange=(stageTranslationStepSize*sample.refractiveIndex)*([1:size(recordedImageStack,3)]-floor(size(recordedImageStack,3)/2+1))/detector.framesPerSecond; %Translation range (along z-axis)
    zRange=tRange; % detection axis, zero centered
    xRange=single(xRange); yRange=single(yRange); zRange=single(zRange); tRange=single(tRange);

%     %Crop the first two dimensions by half
%     recordedImageStack=recordedImageStack([1:floor(end/2)]+floor(end/4)+focusCoordinateDecenter(1),[1:floor(end/2)]+floor(end/4)+focusCoordinateDecenter(2),:);
    %Shift in z (the scan direction) after deconvolution
    %Crop the first dimension by 1/8th and the second by a half
    %recordedImageStack=recordedImageStack([1:floor(7*end/8)]+floor(end/16)+focusCoordinateDecenter(1),[1:floor(end/2)]+floor(end/4)+focusCoordinateDecenter(2),:);
    recordedImageStack=recordedImageStack(:,[1:floor(end/2)]+floor(end/4)+focusCoordinateDecenter(2),:);
    if (focusCoordinateDecenter(1)~=0)
        logMessage('Warning: light sheet not assumed to be centered!!');
        %recordedImageStack=circshift(recordedImageStack,[focusCoordinateDecenter(1) 0 0]);
        xRange=xRange-focusCoordinateDecenter(1)*detector.pixelSize(1)/realDetectionMagnification; % left/right;
    end
    
    %% Preparative Calculations
    
    % Define some noise related variables
    detectionObjectiveSolidAngle=2*pi*(1-cos(asin(detection.objective.numericalAperture/sample.refractiveIndex)));
    objectiveEfficiency=detectionObjectiveSolidAngle/(4*pi);
    overallDetectionEfficiency=sample.fluorophore.quantumYield*objectiveEfficiency*detector.quantumEfficiency;
         
    % Pre-calc. detection PSF
    if (deconvolveWithDetectionPsf)
        logMessage('Calculating detection PSF...');
        detectionPsf=calcDetectionPsf(xRange,yRange,zRange,0,0,detection,sample.refractiveIndex); %Assume shift invariance of the detection PSF
        detectionPsf=single(detectionPsf);
        detectionPsf=overallDetectionEfficiency*detectionPsf;
    else
        detectionPsf=[];
    end
    
    %Pre-calc the light-sheet
    logMessage('Calculating theoretical light sheet...');
    photonEnergy=hPlanck*cLight/excitation.wavelength;
    lightSheetPsf=calcLightSheetPsf(xRange,yRange,zRange,0,excitation,alpha,openFractionOfRadius,sample.refractiveIndex,illuminationCroppingFactors);
    lightSheetPsf=lightSheetPsf*excitation.power*detector.integrationTime/photonEnergy; %Convert to photons per voxel
        
    
    %% Image reconstruction
    logMessage('Reconstructing convolved data set...');
    config={}; config.excitation=excitation; config.detection=detection; config.detector=detector; config.sample=sample;
    config.excitation.objective.refractiveIndex=1;
    config.stagePositions={}; config.stagePositions.target=tRange;
    config.modulation={}; config.modulation.alpha=7; config.modulation.beta=1;
    [restoredDataCube lightSheetDeconvFilter lightSheetOtf ZOtf xRange,yRange,zRange tRange lightSheetPsf]=deconvolveRecordedImageStack(recordedImageStack,config);
%     [restoredDataCubeI lightSheetDeconvFilter lightSheetOtf ZOtf tRangeExtended]=reconstructLightSheetDataCube(xRange,yRange,zRange,tRange,recordedImageStack,excitation,detection,lightSheetPsf,detectionPsf,sample.signalLevel,sample.backgroundLevel,deflectBeamInsteadOfSampleMovement,deNoiseFrames);
%     restoredDataCubeI=circshift(restoredDataCubeI,[0 0 round(focusCoordinateDecenter(3))]);
    
    if (nargout==0)
        save('restoredDataCube.mat','restoredDataCubeI','xRange','yRange','zRange','tRange','tRangeExtended');

        %% Display results
        resultFig=figure();
        resultAx(1)=subplot(2,2,1);
        resultAx(2)=subplot(2,2,2);
        resultAx(3)=subplot(2,2,3);
        resultAx(4)=subplot(2,2,4);
        %The recorded image
        showImage(squeeze(sum(recordedImageStack,2)).',.20,xRange*1e6,tRange*1e6,resultAx(1));
        set(gca,'FontWeight','bold','FontName','Times','FontSize',18,'LineWidth',3);
        xlabel('x [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Times','FontSize',18); ylabel('z [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Times','FontSize',18);
        title(sprintf('%s',lightSheetName),'FontWeight','bold','FontName','Times','FontSize',24);
        centerImageView();
        %The restored image
        showImage(squeeze(sum(restoredDataCubeI,2)).',.05,xRange*1e6,tRangeExtended*1e6,resultAx(2));
        set(gca,'FontWeight','bold','FontName','Times','FontSize',18,'LineWidth',3);
        xlabel('x [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Times','FontSize',18); ylabel('z [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Times','FontSize',18);
        title(sprintf('%s restored',lightSheetName),'FontWeight','bold','FontName','Times','FontSize',24);
        centerImageView();
        %The light-sheet
        showImage(squeeze(lightSheetPsf(:,floor(end/2)+1,:)).'./max(abs(lightSheetPsf(:))),[],xRange*1e6,zRange*1e6,resultAx(3));
        set(gca,'FontWeight','bold','FontName','Times','FontSize',18,'LineWidth',3);
        xlabel('x [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Times','FontSize',18); ylabel('z [\mum]','Interpreter','Tex','FontWeight','bold','FontName','Times','FontSize',18);
        title(sprintf('%s light-sheet',lightSheetName),'FontWeight','bold','FontName','Times','FontSize',24);
        centerImageView();
        %The MTF
        axes(resultAx(4));
        plotSliceXPosIdx=floor(length(xRange)/2)+1;
    %     plotSliceXPosIdx=find(xRange==0);
        mtf=abs(squeeze(lightSheetOtf(plotSliceXPosIdx,floor(end/2)+1,[floor(end/2)+1:end,1]) ));
        filterAmplification=abs(squeeze(lightSheetDeconvFilter(plotSliceXPosIdx,floor(end/2)+1,[floor(end/2)+1:end,1])));
        mtfRestored=mtf.*filterAmplification;
        noiseLevel=0.01;
        noiseAmplification=noiseLevel*filterAmplification;
        ZOtfRange=-1e-6*ZOtf(floor(end/2)+1:-1:1);
        area(ZOtfRange,noiseAmplification,'LineWidth',3,'FaceColor',[.5 .5 .5],'EdgeColor','none');
        hold on;
        plot(ZOtfRange,mtf,'Color',[.8 0 0],'LineWidth',3,'LineStyle','-');
        plot(ZOtfRange,mtfRestored,'Color',[0 .75 .1],'LineWidth',2,'LineStyle','-');
        xlim([0 ZOtfRange(end)]);
        ylim([0 1.1]);
        hold off;
        set(gca,'FontWeight','bold','FontName','Times','FontSize',18,'LineWidth',3);
        xlabel('\nu_z [cycles/\mum]','Interpreter','Tex','FontWeight','bold','FontName','Times','FontSize',18); ylabel('MTF','Interpreter','Tex','FontWeight','bold','FontName','Times','FontSize',18);
        title(sprintf('%s Deconvolution MTF',lightSheetName),'FontWeight','bold','FontName','Times','FontSize',24);


    %     showSlices(recordedImageStack./max(abs(recordedImageStack(:))),0.1,xRange,yRange,tRange);
    
        clear restoredDataCubeI;
    else
        if (nargout>4)
            lightSheetPsf=squeeze(lightSheetPsf(:,1+floor(end/2),:));
            lightSheetOtf=squeeze(lightSheetOtf(:,1+floor(end/2),:));
            lightSheetDeconvFilter=squeeze(lightSheetDeconvFilter(:,1+floor(end/2),:));
        end
    end
    
end

function centerImageView()
    xlim([-50 50]);
    ylim([-20 20]);
    axis equal;
    set(gca,'XTick',[-50:10:50]);
    set(gca,'YTick',[-20:10:20]);
end
