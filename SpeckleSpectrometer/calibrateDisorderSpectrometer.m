% calibrateSpectrometer(wavelengths,slm,cam,outputFileName)
%
% 
function calibrateDisorderSpectrometer(wavelengths,slm,cam,source,outputFileName,emailsToNotify)
    if (nargin<1 || isempty(wavelengths))
        wavelengths=[550 650 600]*1e-9;
        % Put them in a rectangle
        factors=factor(numel(wavelengths));
        gridSize(1)=prod(factors([1:floor(end/4) (floor(3*end/4)+1):end]));
        gridSize(2)=numel(wavelengths)/gridSize(1);
        wavelengths=reshape(wavelengths,gridSize).';
        logMessage('Creating a grid of %dx%d for %d wavelengths.',[gridSize numel(wavelengths)]);
    end
    if (nargin<2 || isempty(slm))
        slm=PhaseSLM(1); %Use the whole of display 1
        slm.referenceDeflectionFrequency=[1/10 1/10];
    end
    if (nargin<3 || isempty(cam))
        cam=BaslerGigECam();
        cam.integrationTime=3e-3;
        cam.gain=1;
    end
    if (nargin<4 || isempty(source))
        % Initialize the SuperK
        source=SuperK();
        source.targetPower=0.50; % Half power
    end
    if (nargin<5 || isempty(outputFileName))
        outputFileName=[pwd(),'/spectrometerCalibration_',datestr(now(),'YYYY-mm-DD_HH-MM'),'.mat'];
    end
    progressEmailAddresses={'tv2@st-andrews.ac.uk','adf10@st-andrews.ac.uk'};
    if (nargin<6 || isempty(emailsToNotify))
        emailsToNotify=progressEmailAddresses;
    end
    
    nbWavelengths=size(wavelengths);
    
    baseWavelength=500e-9;
    
    probeGridSize=[25 25 3];
    
    function cont=progressFunctor(fractionDone)
        cont=true;
        if (floor(fractionDone*100)>floor(prevFractionDone*100))
            logMessage('%0.0f%% done.',100*fractionDone);
            prevFractionDone=fractionDone;
        end
    end

    baseDeflectionFrequency=slm.referenceDeflectionFrequency*baseWavelength;
    baseTwoPiEquivalent=slm.twoPiEquivalent/baseWavelength;
    
    initialCorrection=slm.correctionFunction;
    
    gridSpacing=[1 1]*min(cam.regionOfInterest(3)/nbWavelengths(1),cam.regionOfInterest(4)/nbWavelengths(2));
    [targetPosY targetPosX]=ndgrid(round(cam.regionOfInterest(3)/2+([1:nbWavelengths(1)]-nbWavelengths(1)/2-1/2)*gridSpacing(1)),round(cam.regionOfInterest(4)/2+([1:nbWavelengths(2)]-nbWavelengths(2)/2-1/2)*gridSpacing(2)));
    
    minutesToWait=input('Enter the number of minutes to wait before starting the experiment: ');
    if (~isempty(minutesToWait))
        logMessage('Waiting for %0.0f minutes to start experiment.',minutesToWait,emailsToNotify);
        originalSetPower=source.targetPower;
        pause(minutesToWait*60);
        logMessage('Starting measurement...');
        source.targetPower=originalSetPower;
    else
        logMessage('Starting measurement immediately...');
    end
    
    source.targetPower=0.50; % Half power
    
    transferMatrix=zeros([slm.regionOfInterest(3:4) nbWavelengths]);
    for (wavelengthIdx=1:prod(nbWavelengths))
        wavelength=wavelengths(wavelengthIdx);
        targetPos=[targetPosY(wavelengthIdx) targetPosX(wavelengthIdx)];
        %Adjust the SLM for the new wavelength
        slm.referenceDeflectionFrequency=baseDeflectionFrequency/wavelength; % Keep tilt the same on the sample
        slm.twoPiEquivalent=baseTwoPiEquivalent*wavelength; % Keep efficiency identical
        
        logMessage('Setting the wavelength to %0.3f nm',wavelength*1e9);
        source.setWavelengths(wavelength,0.5); % half the modulation amplitude
    
        prevFractionDone=0;
        [measuredPupilFunction eigenVector probeField4DMatrix sampleX sampleY]=aberrationMeasurement(slm,probeGridSize,@() probeFunctor(cam,targetPos),@progressFunctor);
        transferMatrix(:,:,wavelengthIdx)=measuredPupilFunction./initialCorrection;

        save(outputFileName,'wavelengths','transferMatrix');
        logMessage('Probing for wavelength %0.0fnm done.',wavelength*1e9,progressEmailAddresses);
        logMessage('Partial output written to %s.',outputFileName);
    end
    
    % Calculate the final mask
    amplificationLimit=1;
    spectrometerSuperpositionMask=zeros(slm.regionOfInterest(3:4));
    for (wavelengthIdx=1:prod(nbWavelengths))
        pupilFunctionCorrection=calcCorrectionFromPupilFunction(transferMatrix(:,:,wavelengthIdx).*conj(initialCorrection),amplificationLimit);
        spectrometerSuperpositionMask=spectrometerSuperpositionMask+pupilFunctionCorrection;
    end
    spectrometerSuperpositionMask=spectrometerSuperpositionMask./max(abs(spectrometerSuperpositionMask(:)));
    
    % Save the results
    pupilFunctionCorrection=spectrometerSuperpositionMask;
    initialCorrection=slm.correctionFunction;
    referenceDeflectionFrequency=slm.referenceDeflectionFrequency;
    slmRegionOfInterest=slm.regionOfInterest;
    twoPiEquivalent=slm.twoPiEquivalent;
    save(outputFileName,'wavelengths','transferMatrix','pupilFunctionCorrection','referenceDeflectionFrequency','slmRegionOfInterest','twoPiEquivalent','initialCorrection','gridSpacing','targetPosX','targetPosY');
    logMessage('Experiment done, output written to %s.',outputFileName,emailsToNotify); %,outputFileName);
    
    
%     logMessage('Done, displaying results now... (Close figure window to exit)');
%     testDisorderSpectrometer(transferMatrix,wavelengths,cam,slm,source);
    
    source.delete(); % Shut down
end

function [value auxValues]=probeFunctor(cam,centerPos)
    probeSize=[3 3];
    img=cam.acquire();
    vals=img(centerPos(1)-cam.regionOfInterest(1)+[-floor(probeSize(1)/2):floor((probeSize(1)-1)/2)],centerPos(2)-cam.regionOfInterest(2)+[-floor(probeSize(1)/2):floor((probeSize(1)-1)/2)]);
    value=mean(vals(:));
    
    auxValues=[];
end