function calibrateSpeckleSpectrometer(wavelengths)
    if (nargin<1 || isempty(wavelngths))
        wavelengths=[500:50:600]*1e-9;
    end
    
    cam=BaslerGigECam();
    cam.integrationTime=20e-3;
    cam.gain=1;
    centerPos=floor(1+cam.maxSize/2);
    
    slm=PhaseSLM(1);
    slm.referenceDeflectionFrequency=[1/10 1/10];
    
    cam=selectRegionOfInterestAroundPeakIntensity(cam,slm,[128 128],centerPos);
    
    for wavelengthIdx=1:length(wavelengths)
        wavelength=wavelenths(wavelengthIdx);
        logMessage('Please set the laser to %0.3fnm',wavelength);
        pause();
        
        cam=selectRegionOfInterestAroundPeakIntensity(cam,slm,min(cam.regionOfInterest(3:4),[1 1]*32),centerPos);
    
    end

end