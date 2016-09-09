%
%
function testDisorderSpectrometer(transferMatrix,wavelengths,cam,slm,source)
    if (nargin<1)
        logMessage('Please specify the transferMatrix!');
        return;
    end
    if (nargin<2 || isempty(wavelengths))
        nbWavelengths=[size(transferMatrix,3) size(transferMatrix,4)];
        wavelengths=[500+[1:prod(nbWavelengths)]]*1e-9;
        wavelengths=reshape(wavelengths,nbWavelengths);
    end
    if (nargin<3 || isempty(slm))
        slm=PhaseSLM(1); %Use the whole of display 1
        slm.referenceDeflectionFrequency=[1/10 1/10];
    end
    if (nargin<4 || isempty(cam))
        cam=BaslerGigECam();
        cam.integrationTime=10e-3;
        cam.gain=10;
    end
    if (nargin<5 || isempty(source))
        % Initialize the SuperK
        source=SuperK();
        source.targetPower=0.50; % Half power
    end
    
    % Adjust the grating so that the first order diffraction always goes
    % through the same spot in the iris, independently of the wavelength
    baseWavelength=500e-9;
    baseDeflectionFrequency=slm.referenceDeflectionFrequency*baseWavelength;
    baseTwoPiEquivalent=slm.twoPiEquivalent/baseWavelength;
    
    % Calculate the final mask
    nbWavelengths=size(wavelengths);
    amplificationLimit=3;
    spectrometerSuperpositionMask=zeros(slm.regionOfInterest(3:4));
    for (wavelengthIdx=1:prod(nbWavelengths))
        pupilFunctionCorrection=calcCorrectionFromPupilFunction(transferMatrix(:,:,wavelengthIdx),amplificationLimit);
        spectrometerSuperpositionMask=spectrometerSuperpositionMask+pupilFunctionCorrection;
    end
    
    % Load the mask permanently onto the SLM
    slm.correctionFunction=spectrometerSuperpositionMask;
          
    % Display
    fig=figure;
    try
        idx=1;
        while(ishandle(fig))
            wavelength=wavelengths(idx);
            %Adjust the SLM for the new wavelength
            slm.referenceDeflectionFrequency=baseDeflectionFrequency/wavelength; % Keep tilt the same on the sample
            slm.twoPiEquivalent=baseTwoPiEquivalent*wavelength; % Keep efficiency identical
            slm.modulate(1); % Update the SLM with the new grating
            % Set the source wavelength now
            source.setWavelengths(wavelength);
            img=cam.acquire();
            if (ishandle(fig))
                figure(fig);
                imagesc(img);
                title(sprintf('Wavelength: %0.1f nm',wavelength*1e9));
                drawnow();
                pause(.05);
                idx=1+mod(idx,numel(wavelengths));
            end
        end
    catch Exc
        logMessage(Exc);
    end
end
