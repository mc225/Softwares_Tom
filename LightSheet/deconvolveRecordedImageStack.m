% [restoredDataCube lightSheetDeconvFilter lightSheetOtf ZOtf xRange,yRange,zRange tRange lightSheetPsf]=deconvolveRecordedImageStack(recordedImageStack,config)
%
% recordedImageStack: The recorded values [x y z]=[down right back]
% config: a light sheet microscope set-up configuration file, including the wavelengths and the pupil functions. 
%
function [restoredDataCube lightSheetDeconvFilter lightSheetOtf ZOtf xRange,yRange,zRange tRange lightSheetPsf]=deconvolveRecordedImageStack(recordedImageStack,config)
    % Sample grid specification
    stageTranslationStepSize=norm(median(diff(config.stagePositions.target)));
    
    realMagnification=config.detection.objective.magnification*config.detection.tubeLength/config.detection.objective.tubeLength;
    xRange=-config.detector.center(1)+config.detector.pixelSize(1)*([1:size(recordedImageStack,1)]-floor(size(recordedImageStack,1)/2)-1)/realMagnification; % up/down
    yRange=-config.detector.center(2)+config.detector.pixelSize(2)*([1:size(recordedImageStack,2)]-floor(size(recordedImageStack,2)/2)-1)/realMagnification; % left/right
    tRange=stageTranslationStepSize*([1:size(recordedImageStack,3)]-floor(size(recordedImageStack,3)/2+1))/config.detector.framesPerSecond; %Translation range (along z-axis)
    zRange=tRange; % detection axis, zero centered
%     xRange=single(xRange); yRange=single(yRange); zRange=single(zRange); tRange=single(tRange);
        
    tilt=0;
    logMessage('Calculating light sheet...');
    lightSheetPsf=calcLightSheetPsf(single(xRange),single(yRange),single(zRange),tilt,config.excitation,config.modulation.alpha,config.modulation.beta,config.sample.refractiveIndex);
    
    if (isfield(config.detector,'perspectiveScaling'))
        scaling=config.detector.perspectiveScaling;
    else
        scaling=[];
    end
    
    if (~isempty(scaling))
        % Geometrically correcting recorded data cube and light sheet
        lightSheetPsfOrig=lightSheetPsf;
        logMessage('Geometrically deforming recorded date cube and light sheet with magnification change of [%0.3f%%,%0.3f%%] per micrometer...',scaling*1e-6*100);
        for (zIdx=1:size(recordedImageStack,3))
            zPos=zRange(zIdx);
            sampleXRange=xRange*(1-scaling(1)*zPos);
            sampleYRange=yRange*(1-scaling(2)*zPos);
            recordedImageStack(:,:,zIdx)=interp2(yRange.',xRange,recordedImageStack(:,:,zIdx),sampleYRange.',sampleXRange,'*cubic',0);
            lightSheetPsf(:,:,zIdx)=interp1(yRange,lightSheetPsf(:,:,zIdx),sampleYRange,'*cubic',0);
        end
    end
    
    logMessage('Reconstructing convolved data set...');
    [recordedImageStack lightSheetDeconvFilter lightSheetOtf ZOtf tRange]=deconvolveLightSheet(xRange,yRange,zRange,tRange,recordedImageStack,config,lightSheetPsf);
    restoredDataCube=recordedImageStack; clear recordedImageStack; % This operation does not take extra memory in Matlab
    
    if (~isempty(scaling))
        logMessage('Undoing the geometrically deformation of the date cube and light sheet with magnification change of [%0.3f%%,%0.3f%%] per micrometer...',scaling*1e-6*100);
        for (tIdx=1:size(restoredDataCube,3))
            tPos=tRange(tIdx);
            sampleXRange=xRange*(1-scaling(1)*tPos);
            sampleYRange=yRange*(1-scaling(2)*tPos);
            restoredDataCube(:,:,tIdx)=interp2(sampleYRange.',sampleXRange,restoredDataCube(:,:,tIdx),yRange.',xRange,'*cubic',0);
        end
        lightSheetPsf=lightSheetPsfOrig;
    end
end

%
% Deconvolution in Z after extending the last pixels outward. Note that the
% PSF should occupy maximum 1/2 of the total input z-range to avoid wrapping.
%
function [restoredDataCube lightSheetDeconvFilter lightSheetOtf ZOtf tRange]=deconvolveLightSheet(xRange,yRange,zRange,tRange,recordedImageStack,config,lightSheetPsf)
    % Calculate an extended range over which the deconvolved image can be
    % spread out due to the light sheet diffraction.
    if (length(tRange)>1), dt=diff(tRange([1:2])); else dt=1; end
    tRangePadded=tRange(floor(end/2)+1) + dt*([1:length(tRange)*2]-length(tRange)-1);
    
    inputSize=size(recordedImageStack);
    if (length(inputSize)<4)
        if (length(inputSize)<3)
            inputSize(3)=1;
        end
        inputSize(4)=1;
    end
    
    %Just (sub-pixel)-shift the PSF to the center and pad to double the size 
    if (length(zRange)>1)
        lightSheetDeconvPsf=interp1(-zRange,squeeze(lightSheetPsf(floor(end/2)+1,:,:)).',tRangePadded,'*cubic',0).';
        lightSheetDeconvPsf=permute(lightSheetDeconvPsf,[3 1 2]);
    else
        lightSheetDeconvPsf=lightSheetPsf;
        lightSheetDeconvPsf(:,:,2)=0;
    end
    deconvSize=size(lightSheetDeconvPsf); deconvSize(2)=size(recordedImageStack,2);
    
    %Deconvolve the light-sheet in Z
    ZOtf=([1:deconvSize(3)]-floor(deconvSize(3)/2)-1)/(dt*deconvSize(3));
    numericalApertureInAir=(config.excitation.fractionOfNumericalApertureUsed*config.excitation.objective.numericalAperture)/config.excitation.objective.refractiveIndex;
    excitationOpticalCutOffSpFreq=2*numericalApertureInAir/config.excitation.wavelength;
    excitationNoiseToSignalRatio=config.sample.backgroundLevel*(ZOtf/excitationOpticalCutOffSpFreq)/config.sample.signalLevel;

    lightSheetOtf=fftshift(fft(ifftshift(lightSheetDeconvPsf(1,:,:),3),[],3),3);
    clear lightSheetDeconvPsf;
    lightSheetOtf=lightSheetOtf./max(abs(lightSheetOtf(:))); % Maintain the mean brightness
    %Construct the deconvolution filter:
    lightSheetDeconvFilter=conj(lightSheetOtf)./(abs(lightSheetOtf).^2+repmat(permute(excitationNoiseToSignalRatio.^2,[1 3 2]),[size(lightSheetOtf,1) size(lightSheetOtf,2) 1]));
    
    %Extend edges to double the size and convolve (slice by slice to avoid memory issues)
    % restoredDataCube=zeros([inputSize(1:2), inputSize(3), inputSize(4:end)],'single'); % overwrite recordedImageStack instead to save memory
    for (xIdx=1:size(recordedImageStack,1))
        recordedImageStackSliceFft=fft(recordedImageStack(xIdx,:,[1:end, end*ones(1,floor(end/2)), ones(1,floor((end+1)/2))],:),[],3); % Extend edges
        restoredDataCubeSlice=ifft(recordedImageStackSliceFft.*repmat(ifftshift(lightSheetDeconvFilter,3),[1 1 1 inputSize(4:end)]),[],3,'symmetric');

        %Drop the bit that might have overlap from a wrapped kernel
        recordedImageStack(xIdx,:,:,:)=restoredDataCubeSlice(:,:,1:end/2,:); % overwrite recordedImageStack to save memory
    end
    
    restoredDataCube=recordedImageStack; clear recordedImageStack; % This operation does not take extra memory in Matlab
end
