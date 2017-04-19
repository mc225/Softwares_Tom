function measureDiffusion(wavelength,offsetRange,slm,cam,source,outputFileName)
    if (nargin<1 || isempty(wavelength))
        wavelength=600e-9;
    end
    if (nargin<2 || isempty(offsetRange))
        offsetRange=[-90:1:40]; %[-75:1:75];
    end
    if (nargin<3 || isempty(slm))
        slm=PhaseSLM(1); %Use the whole of display 1
        slm.referenceDeflectionFrequency=[1/10 1/10];
    end
    if (nargin<4 || isempty(cam))
        cam=BaslerGigECam();
        cam.integrationTime=50e-3;
        cam.gain=50;
        cam.numberOfFramesToAverage=4;
%         cam.acquireBackground();
    end
    if (nargin<5 || isempty(source))
        % Initialize the SuperK
        source=SuperK();
        source.targetPower=1.0;
    end
    if (nargin<6 || isempty(outputFileName))
        outputFileName='scanVertical600nmGaPGeEdgeScratchDiffusion6_rectangle20micron';
    end
    lineHeight=20; % pixels are 2 micron on sample
    lineWidth=100; % pixels
    
    nbVerticalOffsets=numel(offsetRange);
    
    % Adjust the grating so that the first order diffraction always goes
    % through the same spot in the iris, independently of the wavelength
    baseWavelength=500e-9;
    baseDeflectionFrequency=slm.referenceDeflectionFrequency*baseWavelength;
    baseTwoPiEquivalent=slm.twoPiEquivalent/baseWavelength;
    % Set the source to the wavelength
    source.setWavelengths(wavelength);
    %Adjust the SLM for the wavelength
    slm.referenceDeflectionFrequency=baseDeflectionFrequency/wavelength; % Keep tilt the same on the sample
    slm.twoPiEquivalent=baseTwoPiEquivalent*wavelength; % Keep efficiency identical
        
    % Open the video output file
    recordingObj = VideoWriter([outputFileName,'.avi'],'Uncompressed AVI'); %'Motion JPEG AVI'); %'Uncompressed AVI');
    recordingObj.FrameRate=25;
    open(recordingObj);
    
    images=zeros([cam.regionOfInterest(3:4) nbVerticalOffsets],'single');
    measuredIntensities=zeros(1,nbVerticalOffsets);
    for (offsetIdx=1:length(offsetRange))
        offset=offsetRange(offsetIdx);
        logMessage('Doing offset %0.0f.',offset);
        %slm.modulate(@(X,Y) sqrt(X.^2+(Y-offset).^2)<50);
        slm.modulate(@(X,Y) abs(X)<round(lineWidth/2)&abs(Y-offset)<round(lineHeight/2));
        img=cam.acquire();
        writeVideo(recordingObj,max(0,img./max(img(:))));
        images(:,:,offsetIdx)=img;
        measuredIntensities(offsetIdx)=mean(img(:));
    end
    
    save([outputFileName,'.mat'],'images','wavelength','measuredIntensities','offsetRange','lineWidth','lineHeight');
    
    close(recordingObj);    
end