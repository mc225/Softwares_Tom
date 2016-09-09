function recordScannedWavelength(wavelengths,slm,cam,source,outputFileName)
    if (nargin<1 || isempty(wavelengths))
        wavelengths=[450:1:750]*1e-9;
    end
    if (nargin<2 || isempty(slm))
        slm=PhaseSLM(1); %Use the whole of display 1
        slm.referenceDeflectionFrequency=[1/10 1/10];
    end
    if (nargin<3 || isempty(cam))
        cam=BaslerGigECam();
        cam.integrationTime=30e-3;
        cam.gain=10;
    end
    if (nargin<4 || isempty(source))
        % Initialize the SuperK
        source=SuperK();
        source.targetPower=1.00;
    end
    if (nargin<5 || isempty(outputFileName))
        outputFileName='scan500to700nm_inAir.avi';
    end
    
    IC=load('C:\Users\nkt\Documents\MATLAB\intensityCorrectionOrig.mat');
    source.setIntensityCorrection(IC.wavelengths,IC.intensityCorrection);
    
    waitForUserSpecifiedTime();
    
    % Adjust the grating so that the first order diffraction always goes
    % through the same spot in the iris, independently of the wavelength
    baseWavelength=500e-9;
    baseDeflectionFrequency=slm.referenceDeflectionFrequency*baseWavelength;
    baseTwoPiEquivalent=slm.twoPiEquivalent/baseWavelength;
        
    % Open the video output file
    recordingObj = VideoWriter(outputFileName,'Uncompressed AVI'); %'Motion JPEG AVI'); %'Uncompressed AVI');
    recordingObj.FrameRate=25;
    open(recordingObj);
    
    measuredIntensities=[];
    for (wavelength=wavelengths)
        logMessage('Doing wavelength %0.1f.',wavelength*1e9);
        %Adjust the SLM for the new wavelength
        slm.referenceDeflectionFrequency=baseDeflectionFrequency/wavelength; % Keep tilt the same on the sample
        slm.twoPiEquivalent=baseTwoPiEquivalent*wavelength; % Keep efficiency identical
        slm.modulate(1);
        % Set the source to the new wavelength
        source.setWavelengths(wavelength);
        img=cam.acquire();
        writeVideo(recordingObj,max(0,img./max(img(:))));
        measuredIntensities(end+1)=mean(img(:));
    end
    
    amplificationLimit=10;
    measuredIntensitiesWithoutCorrection=measuredIntensities./source.intensityCorrection(wavelengths);
    maxIntensity=max(abs(measuredIntensitiesWithoutCorrection(:)));
    intensityCorrection=1./max(1,amplificationLimit*measuredIntensitiesWithoutCorrection./maxIntensity);
    save('intensityCorrection.mat','wavelengths','measuredIntensities','measuredIntensitiesWithoutCorrection','intensityCorrection');
    
    close(recordingObj);
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