% [psf, psfField, varargout]=calcVectorialPsf(xRange,yRange,zRange,wavelength,...
%                                pupilFunctorH,pupilFunctorV,...
%                                objectiveNumericalAperture,refractiveIndexOfSample,
%                                objectiveMagnification,objectiveTubeLength,...
%                                projectionDimensions)
%
% Calculates the 3D point spread function at the grid specified by
% x/y/zRange for a pupil function given by pupilFunctorH/V(U,V) where U and V
% are normalized carthesian pupil coordinates, in a medium with
% refractive index refractiveIndexOfSample and an objective. The horizontal
% polarization given by pupilFunctorH is along the first dimension, and the
% vertical given by pupilFunctorV is along the second dimension of the
% output.
%
% This function gives approximately the same results as PSFLab:
% http://onemolecule.chem.uwm.edu/software , though it is significantly
% faster. Please contact Tom Vettenburg in case you suspect discrepancies
% with the theoretical or simply for more information.
%
% Output:
%     psf: the single photon intensity
%     psfField: the electric field components (x,y,z listed in the fourth
%               dimension of the matrix). For non-vectorial calculations
%               (pupilFunctorV==[]), this matrix is of a maximum of three
%               dimensions.
%     varargout: the higher order nonlinear intensities. This is
%     effectivelly the same as psf.^(N-1) unless projectionDimensions is
%     specified.
%
% Inputs:
%     x/y/zRange: The sample position in metric cathesian coordinates.
%                 Singleton dimensions x and y are treated for normalization
%                 as if they were Nyquist sampled.
%     wavelength: The wavelength in the same metricf units.
%     pupilFunctorH: A function returning the complex horizontal pupil field (along X) 
%                    as a function of carthesian coordinates normalized to the pupil radius.
%                    When a scalar is specified, a constant field of this value is assumed for the whole circular pupil.
%                    Default: unity transmission inside the pupil
%     pupilFunctorV: A function returning the complex vertical pupil field (along Y)
%                    as a function of carthesian coordinates normalized to the pupil radius.
%                    A scalar calculation will be done instead when nothing
%                    or the empty list is give. Specify 0 when only
%                    horizontal input fields are required. Don't use the
%                    empty matrix unless scalar calculations are required.
%                    Default: []: scalar calculation.
%     numericalApertureInSample: The objective's numerical aperture (including refractive index, n, of the medium: n*sin(theta)).
%     refractiveIndexOfSample: The refractive index of the medium at the focus.
%                    Cover slip correction is assumed in the calculation, hence
%                    this only scales the sample grid.
%     objectiveMagnification: The objective magnification (default 100x)
%     objectiveTubeLength: The focal length of the tube lens (default: 200mm)
%     projectionDimensions:
%               -when empty ([]), the full data cube is returned
%               -when a vector with integers, the data is integrated along
%                        the dimensions indicated in the vector.
%                   Default: no projection ([])
%
% Example:
%     xRange=[-10:.1:10]*1e-6;yRange=[-10:.1:10]*1e-6;zRange=[-10:.1:10]*1e-6;
%     objectiveNumericalAperture=asin(1./sqrt(2));
%     pupilFunctor=@(U,V) sqrt(U.^2+V.^2)>.9; %Bessel beam with 10% open fraction
%     [psf psfField]=calcVectorialPsf(xRange,yRange,zRange,500e-9,@(U,V) pupilFunctor(U,V)/sqrt(2),@(U,V) 1i*pupilFunctor(U,V)/sqrt(2),objectiveNumericalAperture,1.0,20,200e-3);
%     psfProj=calcVectorialPsf(xRange,yRange,zRange,500e-9,@(U,V) pupilFunctor(U,V)/sqrt(2),@(U,V) 1i*pupilFunctor(U,V)/sqrt(2),objectiveNumericalAperture,1.0,20,200e-3,[2]);
%
%     img=squeeze(psf(:,1+floor(end/2),:)).';
%     img=img./repmat(mean(img,2),[1 size(img,2)]);
%     imagesc(xRange*1e6,zRange*1e6,img);axis equal;xlabel('x [\mu m]');ylabel('z [\mu m]');
%
function [psf, psfField, varargout]=calcVectorialPsf(xRange,yRange,zRange,wavelength,pupilFunctorH,pupilFunctorV,objectiveNumericalAperture,refractiveIndexOfSample,objectiveMagnification,objectiveTubeLength,projectionDimensions)
    if (nargin<1 || isempty(xRange))
        xRange=[-5000:50:5000]*1e-9;
    end
    if (nargin<2 || isempty(yRange))
        yRange=[-5000:50:5000]*1e-9;
    end
    if (nargin<3 || isempty(zRange))
        zRange=[-5000:500:5000]*1e-9;
    end
    if (nargin<4 || isempty(wavelength))
        wavelength=532e-9;
    end
    if (nargin<5 || isempty(pupilFunctorH))
        % Don't change the following, it is a sensible default!
        pupilFunctorH=@(normalU,normalV) 0.0; % Along the first dimension
    end
    if (nargin<6)
        pupilFunctorV=[]; %Along the second dimension
    end
    if (nargin<7 || isempty(objectiveNumericalAperture))
        objectiveNumericalAperture=.80;
    end
    if (nargin<8 || isempty(refractiveIndexOfSample))
        refractiveIndexOfSample=1.33;
    end
    if (nargin<9 || isempty(objectiveMagnification))
        objectiveMagnification=40;
    end
    if (nargin<10 || isempty(objectiveTubeLength))
        objectiveTubeLength=200e-3;
    end
    if (nargin<11)
        projectionDimensions=[];
    end
    %When scalars are specified instead of function, assume constant input fields
    if (~isa(pupilFunctorH,'function_handle') && isscalar(pupilFunctorH))
        pupilFunctorH=@(normalU,normalV) pupilFunctorH;
    end
    if (~isa(pupilFunctorV,'function_handle') && isscalar(pupilFunctorV))
        pupilFunctorV=@(normalU,normalV) pupilFunctorV;
    end
    vectorialCalculation=~isempty(pupilFunctorV);
    if (~vectorialCalculation)
        logMessage('Starting a scalar calculation of the PSF...');
    end
    
    % Check how many multi-photon orders of the intensity have to be calculated
    highestOrderIntensityRequired=max(1,nargout-1);
    
    focalLengthInSample=objectiveTubeLength/objectiveMagnification; %TODO: Check for correctness
    
    objectiveSinMaxHalfAngleInSample=objectiveNumericalAperture/refractiveIndexOfSample;
    
    % Determine the requested step size for each non-singleton dimension
    sampleDelta=zeros(1,3);
    if (length(xRange)>1), sampleDelta(1)=diff(xRange(1:2)); end
    if (length(yRange)>1), sampleDelta(2)=diff(yRange(1:2)); end
    if (length(zRange)>1), sampleDelta(3)=diff(zRange(1:2)); end
    
    minPupilSize=[1 1]*256; % To ensure that rapid pupil changes are properly sampled
    maxPupilSize=[1 1]*1024; % To avoid memory problems
    
    %The minimum pupil size to avoid PSF replication
    requiredPupilSize=2*ceil([length(xRange) length(yRange)].*sampleDelta(1:2)/(wavelength/(objectiveSinMaxHalfAngleInSample*refractiveIndexOfSample)));
    if (requiredPupilSize>maxPupilSize)
        logMessage('Limiting pupil sampling grid size to (%0.0f,%0.0f) while (%0.0f,%0.0f) would be required in principle.\nThis will cause replicas.',[maxPupilSize,requiredPupilSize]);
    end
    
    %Check if pupil sampling not too sparse for the defocus we intend to simulate
    maxSampleNA=objectiveSinMaxHalfAngleInSample*(1-0.25/(max(maxPupilSize)/2)); % Use of 'max' because the sampling rate near the edge for NA=1 diverges
    minPupilSizeToHandleDefocus=minPupilSize+[1 1]*max(abs(zRange/(wavelength/refractiveIndexOfSample)))*4*objectiveSinMaxHalfAngleInSample*maxSampleNA/sqrt(1-maxSampleNA^2);
    if (minPupilSize<minPupilSizeToHandleDefocus)
        logMessage('A minimum pupil size of (%0.0f,%0.0f) is required to handle the specified defocus.',minPupilSizeToHandleDefocus);
        requiredPupilSize=max(requiredPupilSize,minPupilSizeToHandleDefocus);
    end
    if (minPupilSizeToHandleDefocus>maxPupilSize)
        logMessage('Limiting pupil sampling grid size to (%0.0f,%0.0f) while (%0.0f,%0.0f) would be required in principle. This can cause artefacts.',[maxPupilSize,minPupilSizeToHandleDefocus]);
    end
    
    pupilSize=ceil(min(maxPupilSize,max(minPupilSize,requiredPupilSize)));
    logMessage('Setting the pupil sampling grid size to (%0.0f,%0.0f)',pupilSize);
    
    %Choose the pupil grid
    wavelengthInSample=wavelength/refractiveIndexOfSample;
    uRange=objectiveSinMaxHalfAngleInSample*2*[-floor(pupilSize(1)/2):floor((pupilSize(1)-1)/2)]/pupilSize(1);
	vRange=objectiveSinMaxHalfAngleInSample*2*[-floor(pupilSize(2)/2):floor((pupilSize(2)-1)/2)]/pupilSize(2);
    
    [U,V]=ndgrid(uRange,vRange);
    sinApAngle2=U.^2+V.^2;
    apertureFieldTransmission=double(sinApAngle2<objectiveSinMaxHalfAngleInSample^2);
%     apertureArea=sum(sum(double(sinApAngle2<objectiveSinMaxHalfAngleInSample^2)));
    apertureArea=numel(U)*pi/4;
    sinApAngle=sqrt(apertureFieldTransmission.*sinApAngle2);
    cosApAngle=apertureFieldTransmission.*sqrt(1-apertureFieldTransmission.*sinApAngle2);
    clear sinApAngle2;
    %Scale so that the total intensity is 1 for a unity uniform
    %illumination
    apertureFieldTransmission=apertureFieldTransmission./sqrt(apertureArea);
    
    pupilFunctionX=apertureFieldTransmission.*pupilFunctorH(U/objectiveSinMaxHalfAngleInSample,V/objectiveSinMaxHalfAngleInSample);
    if (vectorialCalculation)
        pupilFunctionY=apertureFieldTransmission.*pupilFunctorV(U/objectiveSinMaxHalfAngleInSample,V/objectiveSinMaxHalfAngleInSample);
        %Convert pupil function to polar coordinates
        T=atan2(V,U); CT=cos(T); ST=sin(T);
        pupilFunctionR =  CT.*pupilFunctionX+ST.*pupilFunctionY; % Radial component is rotated by the focusing
        pupilFunctionA = -ST.*pupilFunctionX+CT.*pupilFunctionY; % Azimutal component is unaffected by the focusing
        %Calculate the polarization change due to focussing
        pupilFunctionZ = sinApAngle.*pupilFunctionR;
        pupilFunctionR = cosApAngle.*pupilFunctionR;
        %Convert back to carthesian coordinates
        pupilFunctionX = CT.*pupilFunctionR-ST.*pupilFunctionA;
        pupilFunctionY = ST.*pupilFunctionR+CT.*pupilFunctionA;
        clear pupilFunctionR pupilFunctionA CT ST T apertureFieldTransmission sinApAngle cosApAngle;

        pupilFunction2D=cat(3,pupilFunctionX,pupilFunctionY,pupilFunctionZ);
        clear pupilFunctionX pupilFunctionY pupilFunctionZ;
    else
        pupilFunction2D=pupilFunctionX;
        clear pupilFunctionX;
    end
    
    %Calculate the focal plain fields
    [psfField psfIntensities]=czt2andDefocus(pupilFunction2D,objectiveSinMaxHalfAngleInSample,xRange/wavelengthInSample,yRange/wavelengthInSample,zRange/wavelengthInSample, focalLengthInSample/wavelengthInSample, projectionDimensions, highestOrderIntensityRequired);

    % Rename the output
    psf=psfIntensities(:,:,:,1);
    if (size(psfIntensities,4)>2)
        varargout=mat2cell(psfIntensities(:,:,:,2:end),size(psfIntensities,1),size(psfIntensities,2),size(psfIntensities,3),ones(1,size(psfIntensities,4)-1));
    else
        if (size(psfIntensities,4)==2)
            varargout={psfIntensities(:,:,:,2)};
        end
    end
    clear psfIntensities;
    
    if (nargout==0)
%         %Store results
%         logMessage('Writing PSF field and intensity to disk...');
%         save('VectorialPsf.mat','psf','psfField','xRange','yRange','zRange','wavelength','pupilFunctorH','pupilFunctorV','pupilFunction2D','objectiveNumericalAperture','refractiveIndexOfSample');

        close all;
        
        %Display results
        maxNormalization=1./max(abs(psf(:)));
        for (zIdx=1:size(psf,3))
            subplot(2,2,1);
            showImage(psf(:,:,zIdx).'.*maxNormalization,[],xRange*1e6,yRange*1e6);
            title(sprintf('Total intensity for z=%0.3f \\mu m',zRange(zIdx)*1e6));
            xlabel('x [\mu m]'); ylabel('y [\mu m]');
            subplot(2,2,2);
            showImage(psfField(:,:,zIdx,3).',-1,xRange*1e6,yRange*1e6);
            title(sprintf('Ez for z=%0.3f \\mu m',zRange(zIdx)*1e6))
            xlabel('x [\mu m]'); ylabel('y [\mu m]');
            subplot(2,2,3);
            showImage(psfField(:,:,zIdx,1).',-1,xRange*1e6,yRange*1e6);
            title(sprintf('Ex for z=%0.3f \\mu m',zRange(zIdx)*1e6))
            xlabel('x [\mu m]'); ylabel('y [\mu m]');
            subplot(2,2,4);
            showImage(psfField(:,:,zIdx,2).',-1,xRange*1e6,yRange*1e6);
            title(sprintf('Ey for z=%0.3f \\mu m',zRange(zIdx)*1e6))
            xlabel('x [\mu m]'); ylabel('y [\mu m]');
            drawnow();
            
            logMessage('Total intensity = %0.3f%%',100*sum(sum(psf(:,:,zIdx))));
            
            pause(1/30);
        end
           
        clear psf; % Don't litter on the command prompt
    end
end

% Calculate the partial spectrum of x using the chirp z transform.
% This returns the complex field.
%
% x is the pupil and should not be ifftshifted, implicit zero padding to the right!
% objectiveSinMaxHalfAngleInSample: the input matrix must cover this disk exactly.
% the following arguments, also specifyable as a list, are:
%     nxRange and nyRange: the sample points in the units of wavelength in the sample medium.
%     nzRange: the sample points in the z dimension in units of wavelength in the sample.
%     focalLengthInSample: (optional, default infinite) The focal length specified in units of wavelength.
%     projectionDimensions: (optional, default none) The dimension along which an integration is done
%     highestOrderIntensityRequired: (optional, default depends on nargout) If specified, (higher order) intensities upto this number are returned as well.
%
% If more than one output argument is specified, the first and higher order
% intensities will be returned as 3D arrays stacked into a single 4D array.
%
function [f, psfIntensities]=czt2andDefocus(x,objectiveSinMaxHalfAngleInSample,nxRange,nyRange,nzRange, focalLengthInSample, projectionDimensions, highestOrderIntensityRequired)
    if (nargin<6 || isempty(focalLengthInSample))
        focalLengthInSample=Inf; %Assuming that focusLength >> z
    end
    if (nargin<7)
        projectionDimensions=[];
    end
    if (nargin<8)
        highestOrderIntensityRequired=max(0,nargout-1);
    end
    
    %Prepare the output matrix with zeros
    inputSize=size(x);
    %inputSize=x(1)*0+inputSize;%Cast to same class as the x input
    outputSize=[length(nxRange) length(nyRange) length(nzRange), inputSize(3:end)]; % Add dimension
    if (~isempty(projectionDimensions))
        outputSize(projectionDimensions)=1;
    end
    f=zeros(outputSize,class(x));
    psfIntensities=zeros([outputSize(1:3) highestOrderIntensityRequired],class(x));
    
    uRange=objectiveSinMaxHalfAngleInSample*2*[-floor(inputSize(1)/2):floor((inputSize(1)-1)/2)]/inputSize(1);
	vRange=objectiveSinMaxHalfAngleInSample*2*[-floor(inputSize(2)/2):floor((inputSize(2)-1)/2)]/inputSize(2);
    [U,V]=ndgrid(uRange,vRange);
    R2=min(1.0,U.^2+V.^2);
    clear U V;
    if (isinf(focalLengthInSample))
        unityDefocusInRad=2*pi*(sqrt(1-R2)-1); %Assuming that focusLength >> z
    end
    
    %Loop through the z-stack
    for zIdx=1:length(nzRange)
        normalizedZ=nzRange(zIdx);
        
        if (~isinf(focalLengthInSample))
            pupilDefocusInRad=2*pi*(sqrt(focalLengthInSample^2+2*focalLengthInSample*sqrt(1-R2)*normalizedZ+normalizedZ^2)-(focalLengthInSample+normalizedZ));
        else
            pupilDefocusInRad=normalizedZ*unityDefocusInRad; %Assuming that focusLength >> z
        end
        pupil=x.*repmat(exp(1i*pupilDefocusInRad),[1 1 inputSize(3:end)]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        psfSlice=czt2fromRanges(pupil,nxRange*2*objectiveSinMaxHalfAngleInSample,nyRange*2*objectiveSinMaxHalfAngleInSample);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Project the output before continuing to save memory
        psfSliceIntensity=sum(abs(psfSlice).^2,3);
        psfSliceIntensity=repmat(psfSliceIntensity,[1 1 highestOrderIntensityRequired]);
        for (photonNb=1:highestOrderIntensityRequired)
            psfSliceIntensity(:,:,photonNb)=psfSliceIntensity(:,:,photonNb).^photonNb;
        end
        if (~isempty(projectionDimensions))
            for (projIdx=1:size(projectionDimensions,2))
                projectionDimension=projectionDimensions(projIdx);
                if (projectionDimension>=3)
                    projectionDimension=projectionDimension-1;
                end
                psfSlice=sum(psfSlice,projectionDimension);
                psfSliceIntensity=sum(psfSliceIntensity,projectionDimension);
            end
        end
        if (any(projectionDimensions==3))
            f(:,:,1,:)=f(:,:,1,:)+psfSlice;
            psfIntensities(:,:,1,:)=psfIntensities(:,:,1,:)+psfSliceIntensity;
        else
            f(:,:,zIdx,:)=psfSlice;
            psfIntensities(:,:,zIdx,:)=psfSliceIntensity;
        end
    end % of for loop over zIdx
end
