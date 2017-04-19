%
%
%
function [sample xRange yRange zRange]=simulateAiryAndBesselSI()
    alpha=7;
    beta=5/100;
    
    xRange=single([-25:.15:25]*1e-6); % propagation
    yRange=single([-20:.15:20]*1e-6); % swipe
    zRange=single([-21:.15:21]*1e-6); % scan
    
    config=struct();
    config.excitation=struct();
    config.excitation.objective=struct();
    config.excitation.objective.refractiveIndex=1.00;
%     config.excitation.objective.numericalAperture=0.42;
    config.excitation.objective.numericalAperture=0.43; % Fig 1C
    config.excitation.objective.magnification=20;
    config.excitation.objective.tubeLength=0.200;
%     config.excitation.wavelength=532e-9;
    config.excitation.wavelength=488e-9; % Fig 1C
	config.excitation.fractionOfNumericalApertureUsed=1;
%     config.sample.refractiveIndex=1.40;
    config.sample.refractiveIndex=1.33; % Fig 1C
    config.detection=struct();
    config.detection.objective=struct();
    config.detection.objective.refractiveIndex=1.00;
%     config.detection.objective.numericalAperture=0.40;
    config.detection.objective.numericalAperture=0.80; % Fig 1C
    config.detection.objective.magnification=20;
    config.detection.objective.tubeLength=0.160;
    config.detection.tubeLength=0.176;
%     config.detection.wavelength=600e-9; %6.12e-7;
    config.detection.wavelength=509e-9; % eGFP, Fig 1C
	config.detection.fractionOfNumericalApertureUsed=1;
    config.detector=struct();
    config.detector.pixelSize=[1 1]*7.4e-6;
    
    % Load synthetic sample
    sample=getTetraedronDensity(xRange,yRange,zRange);
%     sample=getSpiralDensity(xRange,yRange,zRange);

    % Calculate light sheet
    pupilFunctor=@(U,V) U.^2+V.^2>(1-beta)^2;

    [~,psfField]=calcVectorialPsf(yRange,zRange,xRange,config.excitation.wavelength,pupilFunctor,@(U,V) 1i*pupilFunctor(U,V),...
                    config.excitation.objective.numericalAperture,config.sample.refractiveIndex,...
                    config.excitation.objective.magnification,config.excitation.objective.tubeLength);
	psfField=ipermute(psfField,[2 3 1 4]);
    psfFieldSize=size(psfField);
    psfField(1,end*2,1,1)=0; % Pad
	psfFieldFftY=fft(psfField,[],2);
    clear psfField;
    yRangeExt=[yRange yRange(end)+diff(yRange(1:2))*[1:length(yRange)]];

    nbBeams=7;
	nbPhases=5;
    period=1.02e-6*config.excitation.wavelength/488e-9;
    diffractionPeriod=ceil(10e-6/period)*period;

    % Create diffractive grating of nbBeams beams
    logMessage('Passing beam through diffraction grating...');
    gratingFft=0*yRange;
    gratingFft(end*2)=0;
    for diffBeamPos=diffractionPeriod*([1:nbBeams]-floor(1+nbBeams/2)),
        if abs(diffBeamPos)<-2*yRange(1),
            gratingFft=gratingFft+exp(2i*pi*(([1:2*length(yRange)]-1-length(yRange))/(2*length(yRange)))*diffBeamPos/diff(yRange(1:2)));
        end
    end
    gratingFft=ifftshift(gratingFft)./nbBeams;
    
    psf=ifft(psfFieldFftY.*repmat(gratingFft(:).',[psfFieldSize(1) 1 psfFieldSize(3:end)]),[],2);
    clear psfFieldFftY gratingFft;
    psf=sum(abs(psf).^2,4);
    
%     showImage(squeeze(psf(1+floor(end/2),:,:)),-1,zRange*1e6,yRangeExt*1e6);axis equal;
    
	% Create intensity structured illumination and phase shifts
    logMessage('Patterning light sheet intensity...');
    psfFft=fft(psf,[],2);
    clear psf;
    lightSheetFft=zeros([size(psfFft) nbPhases],'single');
    for phaseIdx=1:nbPhases, % Generate phases
        phaseOffset=(phaseIdx-1)/nbPhases;
        
        pulseFft=0*yRange;
        pulseFft(end*2)=0;
        beamPositions=period*(phaseOffset+[1:round(diffractionPeriod/period)]-1-floor(round(diffractionPeriod/period)/2));
        for beamPos=beamPositions,
            pulseFft=pulseFft+exp(2i*pi*(([1:2*length(yRange)]-1-length(yRange))/(2*length(yRange)))*beamPos/diff(yRange(1:2)));
        end
        pulseFft=ifftshift(pulseFft)./length(beamPositions);
        
        lightSheetFft(:,:,:,phaseIdx)=psfFft.*repmat(pulseFft(:).',[psfFieldSize(1) 1 psfFieldSize(3)]);
    end
    clear psfFft;
    lightSheet=ifft(lightSheetFft,[],2,'symmetric');
    clear lightSheetFft;
        
%     for phaseIdx=repmat(1:nbPhases,[1 10])
%         showImage(squeeze(lightSheet(1+floor(end/2),1:end/2,:,phaseIdx)),-1,zRange*1e6,yRange*1e6);axis equal;
%         drawnow();
%         pause(.1);
%     end
    
    lightSheet=lightSheet(:,1:end/2,:,:); % Crop again

    %
    % Calculate detection PSF
    %
    logMessage('Calculating detection PSF...');
    detectionPsf=calcVectorialPsf(xRange,yRange,zRange,config.detection.wavelength,1/sqrt(2),1/sqrt(2),...
                    config.detection.objective.numericalAperture,config.sample.refractiveIndex,...
                    config.detection.objective.magnification,config.detection.objective.tubeLength);
    detectionPsfFft=fft2(ifftshift(ifftshift(detectionPsf,1),2));
    
    %
    % Do the light sheet imaging
    %
    logMessage('Starting light sheet recording...');
    img=zeros([psfFieldSize(1:3) nbPhases]);
    for zIdx=1:length(zRange),
        logMessage('z=%0.1fum in interval [%0.1fum %0.1fum]',[zRange(zIdx)*1e6,zRange([1 end])*1e6]);
        for phaseIdx=1:nbPhases,
            fluorescence=sample(:,:,max(1,min(end,[1:end]+(zIdx-1)-floor(length(zRange)/2)))).*lightSheet(:,:,:,phaseIdx);
            imgSlice=sum(ifft2(fft2(fluorescence).*detectionPsfFft,'symmetric'),3);
            img(:,:,zIdx,phaseIdx)=imgSlice;
        end
    end
    save('img2.mat');
    
    %
    % Deconvolve slice by slice
    %
    [XOtf,YOtf,fRel]=calcOtfGridFromSampleFrequencies(1./[diff(xRange(1:2)) diff(yRange(1:2))],[length(xRange) length(yRange)],(2*config.detection.objective.numericalAperture)/config.detection.wavelength);
    deconvolvedImg=zeros([size(img,1) size(img,2) size(img,3)],'single');
    for zIdx=1:length(zRange),
        slice=img(:,:,zIdx,:);
        orders=[-floor(nbPhases/2):floor(nbPhases/2)];
        order=zeros([size(slice,1) size(slice,2) 1 length(orders)]);
        for orderIdx=1:length(orders),
            grating=repmat((yRange/period),[1 1 1 nbPhases])+repmat(permute(([1:nbPhases]-1)/nbPhases,[1 4 3 2]),[1 length(yRange) 1 1]);
            grating=exp(2i*pi*grating*orders(orderIdx));
            grating=repmat(grating,[length(xRange) 1 1 1]);
            sliceFft=fft2(slice.*grating); %,[],2);
            sliceFft=sum(sliceFft,4); % Remove other orders
            sliceFft=sliceFft.*repmat(ifftshift(fRel)<=1,[1 1 1 size(sliceFft,4)]); % Remove high frequency noise
            order(:,:,1,orderIdx)=fft2(ifft2(sliceFft).*conj(grating(:,:,1,1))); % Shift back
        end
        deconvolvedSlice=ifft(deconvolvedSliceFft,[],2,'symmetric');
        deconvolvedImg(:,:,zIdx)=deconvolvedSlice;
    end

    %
    % Output
    %
    if (nargout==0)
        fig=figure();
%         for zIdx=1:size(sample,3),
%             imagesc(xRange*1e6,yRange*1e6,sample(:,:,zIdx).');
%             axis equal;
%             title(zRange(zIdx)*1e6);
%             drawnow();
%             pause(.1);
%         end
    
        [X,Y,Z]=ndgrid(xRange,yRange,zRange);
        [X,Y,Z, redSample]=reducevolume(X,Y,Z,smooth3(sample,'gaussian',11),[1 1 1]*2);
        isoVal=double(max(redSample(:))/2);
        fv=isosurface(X*1e6,Y*1e6,Z*1e6,redSample,isoVal);
        fvCaps=isocaps(X*1e6,Y*1e6,Z*1e6,redSample,isoVal);
        rfv=reducepatch(fv,0.1);
        rfvCaps=reducepatch(fvCaps,0.1);
        patch(rfv,'FaceColor',[1 0 0],'EdgeColor','none');
        patch(rfvCaps,'FaceColor',[0.5 0 0],'EdgeColor','none');
        axis equal;
        lighting phong;
%         shading interp;
        camlight(30,30);
        camlight(30+180,30);
        
%         close(fig);
        clear sample;
    end
    
end

function sample=getSpiralDensity(xRange,yRange,zRange)
    diameter=20e-6;
    period=2*diameter;
    tubeDiameter=10e-6;
    thickness=1e-6;
    
    nbSpirals=3;
    
    innerR2=abs(tubeDiameter/2-thickness/2)^2;
    outerR2=abs(tubeDiameter/2+thickness/2)^2;
    
    % Make all horizontal
    xRange=xRange(:).'; yRange=yRange(:).'; zRange=zRange(:).';
    
    outputSize=[numel(xRange) numel(yRange) numel(zRange)];
    
    sample=false(outputSize);
    
    for xIdx=1:outputSize(1),
        thetas=2*pi*(xRange(xIdx)/period+([1:nbSpirals].'-1)./nbSpirals);
        spiralPos=[thetas*0 cos(thetas) sin(thetas)]*diameter/2;
        inAnySpiral=false(outputSize(2:3));
        for spiralIdx=1:size(spiralPos,1),
            R2=repmat(abs(yRange-spiralPos(spiralIdx,2)).^2.',[1 outputSize(3)])+repmat(abs(zRange-spiralPos(spiralIdx,3)).^2,[outputSize(2) 1]);
            inSpiral=R2>=innerR2 & R2<=outerR2;
            inAnySpiral=inAnySpiral|inSpiral;
        end
        sample(xIdx,:,:)=inAnySpiral;
    end
    
end

function sample=getTetraedronDensity(xRange,yRange,zRange)
    diameter=20e-6;
    thickness=1e-6;
    centerOffset=[0 -5e-6 -2.5e-6];
    
    % Make all horizontal
    xRange=xRange(:).'; yRange=yRange(:).'; zRange=zRange(:).';
    
    outputSize=[numel(xRange) numel(yRange) numel(zRange)];
    
    innerR2=abs(diameter/2-thickness/2)^2;
    outerR2=abs(diameter/2+thickness/2)^2;
    
    sample=false(outputSize);
    
    spherePos=[0 -1/sqrt(6) 2/sqrt(3); -1 -1/sqrt(6) -1/sqrt(3); 1 -1/sqrt(6) -1/sqrt(3); 0 sqrt(3/2) 0]*diameter/2;
    spherePos=spherePos+repmat(centerOffset,[size(spherePos,1) 1]);
    
    for xIdx=1:outputSize(1),
        inAnyShell=false(outputSize(2:3));
        for sphereIdx=1:size(spherePos,1),
            R2=abs(xRange(xIdx)-spherePos(sphereIdx,1)).^2+...
                repmat(abs(yRange-spherePos(sphereIdx,2)).^2.',[1 outputSize(3)])+repmat(abs(zRange-spherePos(sphereIdx,3)).^2,[outputSize(2) 1]);
            inShell=R2>=innerR2 & R2<=outerR2;
            inAnyShell=inAnyShell|inShell;
        end
        sample(xIdx,:,:)=inAnyShell;
    end
    
    sample=single(sample);
end