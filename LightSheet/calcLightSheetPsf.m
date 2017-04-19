% [psf psfTwoPhotonSwiped psfTwoPhotonCylindricalLens]=calcLightSheetPsf(xRange,yRange,zRange,tilt,excitation,alpha,openFractionOfRadius,refractiveIndexOfSample)
% Calculates the 2D intensity profile created by a swiped light-sheet.
%
% Inputs:
%     tilt: a scalar indicating the tilt of the pupil, or thus the lateral
%           position of the light sheet.
%     excitation: struct with fields wavelength and objective, the latter
%                 must contain the fields magnification, tubeLength, and optionally fractionOfNumericalApertureUsed
%     alpha: the alpha factor of the cubic phase modulation
%     openFractionOfRadius: the open fraction of the annular aperture
%     refractiveIndexOfSample: the refractive index of the sample medium
%
% Returns:
%     psf = the single photon point-spread function (PSF)
%
%  TODO: Check if these arguments still work:
%     psfTwoPhotonSwiped = the two photon PSF assuming that the sheet is created by scanning the beam
%     psfTwoPhotonCylindricalLens = the two photon PSF assuming that the sheet is created with a cylindrical lens
%
%
% Example:
%     zRange=[-50:.1:50]*1e-6;
%     xRange=zRange;
%     yRange=[0:20:100]*1e-6;
%     excitation=struct();
%     excitation.wavelength=532e-9;
%     excitation.objective=struct();
%     excitation.objective.numericalAperture=0.42;
%     excitation.objective.refractiveIndex=1.0;
%     excitation.objective.magnification=20;
%     excitation.objective.tubeLength=200e-3;
%     excitation.objective.illuminationClippingFactors=[1 1; 1 1]*0.0;
%     excitation.gaussianIlluminationStd=2/3;
%     
%     lightSheet=calcLightSheetPsf(xRange,yRange,zRange,0,excitation,0,1,1.4);
%     lightSheet=squeeze(lightSheet).';
% 
%     lightSheetN=lightSheet./repmat(max(lightSheet),[length(zRange) 1]);
%     imagesc(yRange*1e6,zRange*1e6,lightSheetN); axis equal;
%     xlabel('y (propagation) [\mum]');
%     ylabel('z (scan) [\mum]');
%
%     fwhm=[];
%     for idx=1:length(yRange),
%         fwhm(idx)=calcFullWidthAtHalfMaximum(zRange,lightSheet(:,idx),'BiasedLinear');
%     end;
%     close all;
%     plot(yRange*1e6,fwhm*1e6)
%
function [psf psfTwoPhotonSwiped psfTwoPhotonCylindricalLens]=calcLightSheetPsf(xRange,yRange,zRange,tilt,excitation,alpha,openFractionOfRadius,refractiveIndexOfSample)
    if (nargin<5 || isempty(excitation))
        excitation=struct('wavelength',532e-9,...
            'objective',struct('numericalAperture',0.80,'magnification',40,'tubeLength',200e-3),...
            'fractionOfNumericalApertureUsed',1.0);
    end
    numericalAperture=excitation.objective.numericalAperture;
    if (isfield(excitation,'fractionOfNumericalApertureUsed'))
        numericalAperture=numericalAperture*excitation.fractionOfNumericalApertureUsed;
    end
    if (nargin<8 || isempty(refractiveIndexOfSample))
        refractiveIndexOfSample=1.0;
    end
    if (nargin<1 || isempty(xRange))
        xRange=(7.4*1e-6/excitation.objective.magnification)*[-200:199];
        xRange=single(xRange);
    end
    if (nargin<2 || isempty(yRange))
        yRange=(7.4*1e-6/excitation.objective.magnification)*[-200:199];
        yRange=single(yRange);
    end
    if (nargin<3 || isempty(yRange))
        stageTranslationStepSize=0.1*1e-6;
        zRange=(stageTranslationStepSize*refractiveIndexOfSample)*[-250:249]; %Translation range (along z-axis)
        zRange=single(zRange);
    end
    if (nargin<4 || isempty(tilt))
        tilt=0;
    end
    if (nargin<6 || isempty(alpha))
        alpha=0;
    end
    if (nargin<7 || isempty(openFractionOfRadius))
        openFractionOfRadius=1;
    end
    
    % ! Code to support old data sets !
    % Check if the SLM is smaller than the back aperture, if so, simulate the limited illumination
    if (~isfield(excitation,'illuminationClippingFactors'))
        illuminationClippingFactors=0*[1 1; 1 1];
    else
        illuminationClippingFactors=excitation.illuminationClippingFactors;
    end
    if (~isfield(excitation,'gaussianIlluminationStd'))
        gaussianIlluminationStd=[];
    else
        gaussianIlluminationStd=excitation.gaussianIlluminationStd;
    end
    if (~isfield(excitation,'beamAngle'))
        beamAngle=[];
    else
        beamAngle=excitation.beamAngle;
    end
    
    projectionDimension=1;
    if (length(xRange)<=1)
        %Nyquist sampling is sufficient, we will do a projection later anyway
        projRangeLength=512;
        xRangeForProj=([1:projRangeLength]-floor(projRangeLength/2)-1)*0.5*0.5*excitation.wavelength/numericalAperture;
    else
        xRangeForProj=xRange;
    end

    if (openFractionOfRadius>0)
        crop=@(U,V) 1.0*(U>=-(1-illuminationClippingFactors(1,1))&U<=(1-illuminationClippingFactors(1,2)) & V>=-(1-illuminationClippingFactors(2,1))&V<=(1-illuminationClippingFactors(2,2)));
        if (~isempty(gaussianIlluminationStd))
            apodization=@(U,V) exp(-(U.^2+V.^2)./(2*gaussianIlluminationStd^2));
        else
            apodization=@(U,V) 1;
        end
        
        if (~isempty(beamAngle))
            phaseMask=@(U,V) alpha*((cos(beamAngle)*U-sin(beamAngle)*V).^3 + (sin(beamAngle)*U+cos(beamAngle)*V).^3);
        else
            phaseMask=@(U,V) alpha*(U.^3+V.^3);
        end
        
        %Rotate the coordinate system so that X and Z are interchanged
        %pupilFunctor=@(U,V) crop(U,V).*(sqrt(U.^2+V.^2)>=(1-openFractionOfRadius)).*exp(2i*pi*(alpha*V.^3 + (tilt-alpha*3/5)*V));
        pupilFunctor=@(U,V) crop(U,V).*apodization(U,V).*(sqrt(U.^2+V.^2)>=(1-openFractionOfRadius)).*exp(2i*pi*(phaseMask(U,V) + tilt*V));
        %Assume circular polarization propagating along the x-axis
        cache=Cache(); % Get the default cache
        key={xRangeForProj,zRange,yRange,{excitation.wavelength,1/sqrt(2),1i/sqrt(2),illuminationClippingFactors,gaussianIlluminationStd,beamAngle,alpha,openFractionOfRadius,tilt},numericalAperture,refractiveIndexOfSample,excitation.objective.magnification,excitation.objective.tubeLength,projectionDimension};
        if (isfield(cache,key))
            logMessage('Loading PSF projection from cache...');
            psf=cache(key);
        else
            startTime=clock();
            psf=calcVectorialPsf(xRangeForProj,zRange,yRange,excitation.wavelength,@(U,V) pupilFunctor(U,V)/sqrt(2),@(U,V) 1i*pupilFunctor(U,V)/sqrt(2),numericalAperture,refractiveIndexOfSample,excitation.objective.magnification,excitation.objective.tubeLength,projectionDimension);
            if (etime(clock(),startTime)>10) % Only cache slow calculations
                cache(key)=psf;
            end
        end
        if (nargout>=2)
            keyTwoPhoton={key,'psfTwoPhotonSwiped'};
            if (isfield(cache,keyTwoPhoton))
                logMessage('Loading two photon PSF from cache...');
                psfTwoPhotonSwiped=cache(keyTwoPhoton);
            else
                startTime=clock();
                [psfTwoPhoton,~,psfTwoPhotonSwiped]=calcVectorialPsf(xRangeForProj,zRange,yRange,2*excitation.wavelength,@(U,V) pupilFunctor(U,V)/sqrt(2),@(U,V) 1i*pupilFunctor(U,V)/sqrt(2),numericalAperture,refractiveIndexOfSample,excitation.objective.magnification,excitation.objective.tubeLength,projectionDimension);
                if (nargout>=3)
                    psfTwoPhotonCylindricalLens=psfTwoPhoton.^2;
                end
                clear psfTwoPhoton;
                if (etime(clock(),startTime)>10) % Only cache slow calculations
                    cache(keyTwoPhoton)=psfTwoPhotonSwiped;
                end
            end
        end
    else
        error('calcLightSheetPsf.m Not implemented for openFractionOfRadius==0.');
    end
    %Return to original coordinate system
    psf=permute(psf,[1 3 2]);
    if (nargout>=2)
        psfTwoPhotonSwiped=permute(psfTwoPhotonSwiped,[1 3 2]);
        if (nargout>=3)
            psfTwoPhotonCylindricalLens=permute(psfTwoPhotonCylindricalLens,[1 3 2]);
        end
    end
    
    if (nargout==0 && numel(psf)>100)
        figure;
        tmp=squeeze(psf(:,1,:)).';
        tmp=tmp./repmat(max(tmp),[size(tmp,1) 1]);
        imagesc(xRange*1e6,zRange*1e6,tmp);axis equal
        clear psf;
    end
end