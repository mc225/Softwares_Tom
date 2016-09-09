% [psf psfTwoPhotonSwiped psfTwoPhotonCylindricalLens]=calcLightSheetPsf(xRange,yRange,zRange,tilt,excitation,alpha,openFractionOfRadius,refractiveIndexOfSample,illuminationCroppingFactors)
% Calculates the 2D intensity profile created by a swiped light-sheet.
%
% Inputs:
%     excitation: struct with fields wavelength and objective, the latter
%     must contain the fields magnification and tubeLength
%
% Returns:
%     psf = the single photon point-spread function (PSF)
%
%  The following output arguments are obsolete now!
%
%     psfTwoPhotonSwiped = the two photon PSF assuming that the sheet is created by scanning the beam
%     psfTwoPhotonCylindricalLens = the two photon PSF assuming that the sheet is created with a cylindrical lens
%
function [psf psfTwoPhotonSwiped psfTwoPhotonCylindricalLens]=calcLightSheetPsf(xRange,yRange,zRange,tilt,excitation,alpha,openFractionOfRadius,refractiveIndexOfSample,illuminationCroppingFactors)
    if (nargin<5 || isempty(excitation))
        excitation=struct('wavelength',532e-9,...
            'objective',struct('numericalAperture',0.80,'magnification',40,'tubeLength',200e-3),...
            'fractionOfNumericalApertureUsed',1.0);
    end
    numericalAperture=excitation.objective.numericalAperture*excitation.fractionOfNumericalApertureUsed;
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
    if (nargin<9 || isempty(illuminationCroppingFactors))
        illuminationCroppingFactors=0*[1 1; 1 1];
    end
    
    if (nargout<=1)
        projectionDimension=1;
        %Nyquist sampling is sufficient, we will do a projection later anyway
        xRangeForProj=[-128:127]*0.5*0.5*excitation.wavelength/numericalAperture;
    else
        if (length(xRange)<=1)
            xRangeForProj=[-128:127]*0.5*0.5*excitation.wavelength/numericalAperture;
            projectionDimension=2;
        else
            error('calcLightSheetPsf: multiple output arguments are obsolete for non-singleton xRange!');
        end
    end

    if (openFractionOfRadius>0)
        crop=@(U,V) 1.0*(U>=-(1-illuminationCroppingFactors(1,1))&U<=(1-illuminationCroppingFactors(1,2)) & V>=-(1-illuminationCroppingFactors(2,1))&V<=(1-illuminationCroppingFactors(2,2)));
        %Rotate the coordinate system so that X and Z are interchanged
        %pupilFunctor=@(U,V) crop(U,V).*(sqrt(U.^2+V.^2)>=(1-openFractionOfRadius)).*exp(2i*pi*(alpha*V.^3 + (tilt-alpha*3/5)*V));
        pupilFunctor=@(U,V) crop(U,V).*(sqrt(U.^2+V.^2)>=(1-openFractionOfRadius)).*exp(2i*pi*(alpha*V.^3 + (tilt-0*alpha*3/5)*V));
        %Assume circular polarization propagating along the x-axis
        psf=calcVectorialPsf(xRangeForProj,zRange,yRange,excitation.wavelength,@(U,V) pupilFunctor(U,V)/sqrt(2),@(U,V) 1i*pupilFunctor(U,V)/sqrt(2),numericalAperture,refractiveIndexOfSample,excitation.objective.magnification,excitation.objective.tubeLength,projectionDimension);
        if (nargout>=2)
            [psfTwoPhoton psfField psfTwoPhotonSwiped]=calcVectorialPsf(xRangeForProj,zRange,yRange,2*excitation.wavelength,@(U,V) pupilFunctor(U,V)/sqrt(2),@(U,V) 1i*pupilFunctor(U,V)/sqrt(2),numericalAperture,refractiveIndexOfSample,excitation.objective.magnification,excitation.objective.tubeLength,projectionDimension);
            clear psfField;
            if (nargout>=3)
                psfTwoPhotonCylindricalLens=psfTwoPhoton.^2;
            end
            clear psfTwoPhoton;
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