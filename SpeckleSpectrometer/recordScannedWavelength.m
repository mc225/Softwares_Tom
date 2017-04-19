function recordScannedWavelength(wavelengths,slm,cam,source,outputFileName)
    if (nargin<1 || isempty(wavelengths))
        wavelengths=[500:1:700]*1e-9;
    end
    if (nargin<2 || isempty(slm))
        slm=PhaseSLM(1); %Use the whole of display 1
        slm.referenceDeflectionFrequency=[1/10 1/10];
    end
    if (nargin<3 || isempty(cam))
        cam=BaslerGigECam();
        cam.integrationTime=50e-3;
        cam.gain=50;
        cam.numberOfFramesToAverage=4;
%         cam.acquireBackground();
    end
    if (nargin<4 || isempty(source))
        % Initialize the SuperK
        source=SuperK();
        source.targetPower=1.0;
    end
    if (nargin<5 || isempty(outputFileName))
        outputFileName='scan500to700nmGaPGeEdgeScratch';
    end
    
    emailsToNotify={'adf10@st-andrews.ac.uk','tv2@st-andrews.ac.uk'};
    
    slmMask=@(X,Y) sqrt(X.^2+Y.^2)<200;
    
%     IC=load('C:\Users\nkt\Documents\MATLAB\intensityCorrectionOrig.mat');
%     source.setIntensityCorrection(IC.wavelengths,IC.intensityCorrection);
    
    if (nargin==0)
        waitForUserSpecifiedTime();
    end
    
    nbWavelengths=numel(wavelengths);
    
    % Adjust the grating so that the first order diffraction always goes
    % through the same spot in the iris, independently of the wavelength
    baseWavelength=500e-9;
    baseDeflectionFrequency=slm.referenceDeflectionFrequency*baseWavelength;
    baseTwoPiEquivalent=slm.twoPiEquivalent/baseWavelength;
        
    % Open the video output file
    recordingObj = VideoWriter([outputFileName,'.avi'],'Uncompressed AVI'); %'Motion JPEG AVI'); %'Uncompressed AVI');
    recordingObj.FrameRate=25;
    open(recordingObj);
    
    detectorNoise=[];
    differencesWithInitial=zeros(1,nbWavelengths);
    differencesWithPrevious=zeros(1,nbWavelengths);
    measuredIntensities=zeros(1,nbWavelengths);
    images=zeros([cam.regionOfInterest(3:4) nbWavelengths],'single');
    for (wavelengthIdx=1:nbWavelengths)
        wavelength=wavelengths(wavelengthIdx);
        logMessage('Doing wavelength %0.1f.',wavelength*1e9);
        %Adjust the SLM for the new wavelength
        slm.referenceDeflectionFrequency=baseDeflectionFrequency/wavelength; % Keep tilt the same on the sample
        slm.twoPiEquivalent=baseTwoPiEquivalent*wavelength; % Keep efficiency identical
        slm.modulate(slmMask);
        % Set the source to the new wavelength
        source.setWavelengths(wavelength);
        img=cam.acquire();
        writeVideo(recordingObj,max(0,img./max(img(:))));
        images(:,:,wavelengthIdx)=img;
        measuredIntensities(wavelengthIdx)=mean(img(:));
        if (~isempty(detectorNoise))
            differenceImageWithInitial=img./norm(img(:))-initialImage./norm(initialImage(:));
            differencesWithInitial(wavelengthIdx)=norm(differenceImageWithInitial(:));
            differenceImageWithPrevious=img./norm(img(:))-previousImage./norm(previousImage(:));
            differencesWithPrevious(wavelengthIdx)=norm(differenceImageWithPrevious(:));
        else
            initialImage=img;
            nbNoiseImages=100;
            logMessage('Recording %d images to measure detection noise.',nbNoiseImages);
            noiseImages=initialImage;
            noiseImages(:,:,nbNoiseImages)=0;
            for imgIdx=2:nbNoiseImages,
                noiseImages(:,:,imgIdx)=cam.acquire();
            end
            detectorNoise=mean(mean(std(noiseImages,0,3)));
            clear noiseImages;
        end
        previousImage=img;
    end
    
    amplificationLimit=10;
    measuredIntensitiesWithoutCorrection=measuredIntensities./source.intensityCorrection(wavelengths);
    maxIntensity=max(abs(measuredIntensitiesWithoutCorrection(:)));
    intensityCorrection=1./max(1,amplificationLimit*measuredIntensitiesWithoutCorrection./maxIntensity);
    
    save([outputFileName,'.mat'],'images','wavelengths','measuredIntensities','measuredIntensitiesWithoutCorrection','intensityCorrection','differencesWithInitial','differencesWithPrevious','detectorNoise');
    
    close(recordingObj);
    
    logMessage('Wavelength range recording done.',[],emailsToNotify);
    
    % Normalize for power fluctuations
    intensityNorms=squeeze(sqrt(sum(sum(abs(images).^2))));
    images=images./repmat(permute(intensityNorms,[3 2 1]),[size(images,1) size(images,2) 1]);
    meanImage=mean(images,3);
    images=images-repmat(meanImage,[1 1 size(images,3)]);
    unbiasedIntensityNorms=squeeze(sqrt(sum(sum(abs(images).^2))));
    images=images./repmat(permute(unbiasedIntensityNorms,[3 2 1]),[size(images,1) size(images,2) 1]);
    
    images=reshape(images,[],size(images,3));
    angleDeviationFromOrthogonal=pi/2-acos(max(-1,min(1,images'*images)));
    figure;
    showImage(angleDeviationFromOrthogonal/(pi/2)+1e-6i); axis equal; title('Deviation from orthogonal.');
    
    save([outputFileName,'.mat'],'angleDeviationFromOrthogonal','-append');
    
end

function waitForUserSpecifiedTime()
    minutesToWait=input('Enter the number of minutes to wait before starting the experiment: ');
    if (~isempty(minutesToWait))
        logMessage('Waiting for %0.0f minutes to start experiment.',minutesToWait);
        pause(minutesToWait*60);
        logMessage('Starting measurement...');
    else
        logMessage('Starting measurement immediately...');
    end
end