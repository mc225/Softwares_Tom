% 
%
function compareBeamProfiles()
    close all;
    
    config=struct();
    config.excitation=struct();
    config.excitation.objective=struct();
    config.excitation.objective.refractiveIndex=1.00;
    config.excitation.objective.numericalAperture=0.42;
    config.excitation.objective.magnification=20;
    config.excitation.objective.tubeLength=0.200;
    config.excitation.wavelength=532e-9;
	config.excitation.fractionOfNumericalApertureUsed=1;
    config.sample.refractiveIndex=1.40;
    config.detection=struct();
    config.detection.objective=struct();
    config.detection.objective.refractiveIndex=1.00;
    config.detection.objective.numericalAperture=0.40;
    config.detection.objective.magnification=20;
    config.detection.objective.tubeLength=0.160;
    config.detection.tubeLength=0.176;
    config.detection.wavelength=600e-9; %6.12e-7;
	config.detection.fractionOfNumericalApertureUsed=1;
    config.detector=struct();
    config.detector.pixelSize=[1 1]*7.4e-6;
    
%     coreWidthBessel=2*config.excitation.wavelength/(2*config.excitation.objective.numericalAperture);
%     coreWidthBessel=2*config.excitation.wavelength/(2*config.excitation.objective.numericalAperture*0.88); %2*0.720e-6;
    coreWidthBessel=2*2.4048255576957727686216318793264546431242449091459671*config.excitation.wavelength/(2*pi*config.excitation.objective.numericalAperture);
    
    realMagnification=config.detection.objective.magnification*config.detection.tubeLength/config.detection.objective.tubeLength;
    depthOfFocus=config.sample.refractiveIndex*(config.detection.wavelength/(config.detection.objective.numericalAperture^2)+...
        config.detector.pixelSize(1)/(realMagnification*config.detection.objective.numericalAperture));
    spotLength=4*config.sample.refractiveIndex*config.detection.wavelength/(config.detection.objective.numericalAperture^2);
        
    % Gaussian
    resGaussian=plotUsefulPhotonFraction(config,0,1,spotLength);
 
    % Airy
    resAiry=plotUsefulPhotonFraction(config,7,1,spotLength);

    % Bessel 10
    resBessel10=plotUsefulPhotonFraction(config,0,0.10,[coreWidthBessel spotLength]);    
    % Bessel 5
    resBessel5=plotUsefulPhotonFraction(config,0,0.05,[coreWidthBessel spotLength]);
    % Bessel 2
	resBessel2=plotUsefulPhotonFraction(config,0,0.02,[coreWidthBessel spotLength]);      
    % Bessel 1
	resBessel1=plotUsefulPhotonFraction(config,0,0.01,[coreWidthBessel spotLength]);
    
%     % Show two-photon MTF plots
%     figure;
%     plot(resGaussian.fRel2Slice,abs(resGaussian.otf2Slice),'Color',[.5 .5 .5]); hold on;
%     plot(resAiry.fRel2Slice,abs(resAiry.otf2Slice),'Color',[0 .8 0]); hold on;
%     plot(resBessel10(1).fRel2Slice,abs(resBessel10(1).otf2Slice),'Color',[.3 0 0]); hold on;
%     plot(resBessel5(1).fRel2Slice,abs(resBessel5(1).otf2Slice),'Color',[.5 0 0]); hold on;
%     plot(resBessel2(1).fRel2Slice,abs(resBessel2(1).otf2Slice),'Color',[.7 0 0]); hold on;
%     plot(resBessel1(1).fRel2Slice,abs(resBessel1(1).otf2Slice),'Color',[1 0 0]); hold on;
    
    % Display
    fig=figure();
    plot(resGaussian.zRange*1e6,resGaussian.focalPlaneEfficiency(:,1),'Color',[0 0 0],'LineWidth',2);
    hold on;
    plot(resAiry.zRange*1e6,resAiry.focalPlaneEfficiency(:,1),'Color',[0 .5 0],'LineWidth',3);
    plot(resBessel10(1).zRange*1e6,resBessel10(1).focalPlaneEfficiency(:,1),'Color',[1 0 0],'LineWidth',2,'LineStyle','--');
    plot(resBessel5(1).zRange*1e6,resBessel5(1).focalPlaneEfficiency(:,1),'Color',[.75 0 .75],'LineWidth',1,'LineStyle','--');
    plot(resBessel2(1).zRange*1e6,resBessel2(1).focalPlaneEfficiency(:,1),'Color',[0 1 .5],'LineWidth',1,'LineStyle','--');
    plot(resBessel1(1).zRange*1e6,resBessel1(1).focalPlaneEfficiency(:,1),'Color',[0 0 1],'LineWidth',2,'LineStyle','--');
    plot(resBessel10(1).zRange*1e6,resBessel10(1).coreEfficiency(:,1),'Color',[1 0 0],'LineWidth',2,'LineStyle',':');
    plot(resBessel5(1).zRange*1e6,resBessel5(1).coreEfficiency(:,1),'Color',[.75 0 .75],'LineWidth',1,'LineStyle',':');
    plot(resBessel2(1).zRange*1e6,resBessel2(1).coreEfficiency(:,1),'Color',[0 1 .5],'LineWidth',1,'LineStyle',':');
    plot(resBessel1(1).zRange*1e6,resBessel1(1).coreEfficiency(:,1),'Color',[0 0 1],'LineWidth',2,'LineStyle',':');
    plot(resBessel2(1).zRange*1e6,resBessel2(1).focalPlaneEfficiency(:,2),'Color',[0 1 .5],'LineWidth',2,'LineStyle','-.');
    plot(resBessel1(1).zRange*1e6,resBessel1(1).focalPlaneEfficiency(:,2),'Color',[0 0 1],'LineWidth',2,'LineStyle','-.');
    xlim(resBessel1(1).zRange([1 end])*1e6); ylim([0 1]);
    legend({'Gaussian','Airy','Bessel10-SI','Bessel5-SI','Bessel2-SI','Bessel1-SI',...
        'Bessel10-conf','Bessel5-conf','Bessel2-conf','Bessel1-conf','Bessel2-2PE','Bessel1-2PE'});
    
    FOVs=[resGaussian.FOV resAiry.FOV...
        repmat([resBessel10(1).FOV resBessel5(1).FOV resBessel2(1).FOV resBessel1(1).FOV],[1 3])...
        resBessel10(1).FOV2PE resBessel5(1).FOV2PE resBessel2(1).FOV2PE resBessel1(1).FOV2PE];
    
    FWHMs=[resGaussian(1).fullWidthAtHalfMaximumAtWaist resAiry(1).fullWidthAtHalfMaximumAtWaist...
        resBessel10(2).fullWidthAtHalfMaximumAtWaist resBessel5(2).fullWidthAtHalfMaximumAtWaist resBessel2(2).fullWidthAtHalfMaximumAtWaist resBessel1(2).fullWidthAtHalfMaximumAtWaist...
        resBessel10(1).fullWidthAtHalfMaximumAtWaist resBessel5(1).fullWidthAtHalfMaximumAtWaist resBessel2(1).fullWidthAtHalfMaximumAtWaist resBessel1(1).fullWidthAtHalfMaximumAtWaist...
        resBessel10(1).fullWidthAtHalfMaximumAtWaist resBessel5(1).fullWidthAtHalfMaximumAtWaist resBessel2(1).fullWidthAtHalfMaximumAtWaist resBessel1(1).fullWidthAtHalfMaximumAtWaist...
        resBessel10(2).fullWidthAtHalfMaximumAtWaist2PE resBessel5(2).fullWidthAtHalfMaximumAtWaist2PE resBessel2(2).fullWidthAtHalfMaximumAtWaist2PE resBessel1(2).fullWidthAtHalfMaximumAtWaist2PE...
        ];
    
    waistEfficiencies=[resGaussian.focalPlaneEfficiency(1,1) resAiry.focalPlaneEfficiency(1,1)...
        resBessel10(2).focalPlaneEfficiency(1,1) resBessel5(2).focalPlaneEfficiency(1,1) resBessel2(2).focalPlaneEfficiency(1,1) resBessel1(2).focalPlaneEfficiency(1,1)...
        resBessel10(1).focalPlaneEfficiency(1,1) resBessel5(1).focalPlaneEfficiency(1,1) resBessel2(1).focalPlaneEfficiency(1,1) resBessel1(1).focalPlaneEfficiency(1,1)...
        resBessel10(1).coreEfficiency(1,1) resBessel5(1).coreEfficiency(1,1) resBessel2(1).coreEfficiency(1,1) resBessel1(1).coreEfficiency(1,1)...
        resBessel10(2).focalPlaneEfficiency(1,2) resBessel5(2).focalPlaneEfficiency(1,2) resBessel2(2).focalPlaneEfficiency(1,2) resBessel1(2).focalPlaneEfficiency(1,2)];
    
    averageEfficiencies=[resGaussian.averageFocalPlaneEfficiency(1) resAiry.averageFocalPlaneEfficiency(1)...
        resBessel10(2).averageFocalPlaneEfficiency(1) resBessel5(2).averageFocalPlaneEfficiency(1) resBessel2(2).averageFocalPlaneEfficiency(1) resBessel1(2).averageFocalPlaneEfficiency(1)...
        resBessel10(1).averageFocalPlaneEfficiency(1) resBessel5(1).averageFocalPlaneEfficiency(1) resBessel2(1).averageFocalPlaneEfficiency(1) resBessel1(1).averageFocalPlaneEfficiency(1)...
        resBessel10(1).averageCoreEfficiency(1) resBessel5(1).averageCoreEfficiency(1) resBessel2(1).averageCoreEfficiency(1) resBessel1(1).averageCoreEfficiency(1)...
        resBessel10(1).averageCoreEfficiency(2) resBessel5(1).averageCoreEfficiency(2) resBessel2(1).averageCoreEfficiency(2) resBessel1(1).averageCoreEfficiency(2)];
   
    axialResolution=[resGaussian(1).axialResolution resAiry(1).axialResolution...
        resBessel10(2).axialResolution resBessel5(2).axialResolution resBessel2(2).axialResolution resBessel1(2).axialResolution...
        resBessel10(1).axialResolution resBessel5(1).axialResolution resBessel2(1).axialResolution resBessel1(1).axialResolution...
        resBessel10(1).axialResolution resBessel5(1).axialResolution resBessel2(1).axialResolution resBessel1(1).axialResolution...
        resBessel10(2).axialResolution resBessel5(2).axialResolution resBessel2(2).axialResolution resBessel1(2).axialResolution];
    
    numericalAxialResolution=[resGaussian(1).axialResolution1PE resAiry(1).axialResolution1PE...
        resBessel10(2).axialResolution1PE resBessel5(2).axialResolution1PE resBessel2(2).axialResolution1PE resBessel1(2).axialResolution1PE...
        resBessel10(1).axialResolution resBessel5(1).axialResolution resBessel2(1).axialResolution resBessel1(1).axialResolution...
        resBessel10(1).axialResolution resBessel5(1).axialResolution resBessel2(1).axialResolution resBessel1(1).axialResolution...
        resBessel10(2).axialResolution2PE resBessel5(2).axialResolution2PE resBessel2(2).axialResolution2PE resBessel1(2).axialResolution2PE];
    
    peakValues=[resGaussian(1).peakValueAtWaist resAiry(1).peakValueAtWaist...
        resBessel10(2).peakValueAtWaist resBessel5(2).peakValueAtWaist resBessel2(2).peakValueAtWaist resBessel1(2).peakValueAtWaist...
        resBessel10(1).peakValueAtWaist resBessel5(1).peakValueAtWaist resBessel2(1).peakValueAtWaist resBessel1(1).peakValueAtWaist...
        resBessel10(1).peakValueAtWaist resBessel5(1).peakValueAtWaist resBessel2(1).peakValueAtWaist resBessel1(1).peakValueAtWaist...
        resBessel10(2).peakValueAtWaist2PE resBessel5(2).peakValueAtWaist2PE resBessel2(2).peakValueAtWaist2PE resBessel1(2).peakValueAtWaist2PE];
    
    peakValues=peakValues/peakValues(1);
    relativePeakValues=peakValues./waistEfficiencies;
    relativePeakValues=relativePeakValues/relativePeakValues(1);
    
    logMessage(['Light Sheet Type: Gaussian, Airy, Bessel10, Bessel5, Bessel2, Bessel1, Bessel10-SI, Bessel5-SI, Bessel2-SI, Bessel1-SI, Bessel10-conf, Bessel5-conf, Bessel2-conf, Bessel1-conf, Bessel10-2PE, Bessel5-2PE, Bessel2-2PE, Bessel1-2PE']);
    logMessage(['FOV: ', sprintf('%0.0f ',FOVs*1e6),'um']);
    logMessage(['FWHMs: ', sprintf('%0.0f ',FWHMs*1e9),'nm']);
    logMessage(['half-FOV: ', sprintf('%0.0f ',0.5*FOVs*1e6),'um']);
    logMessage(['waistEfficiencies: ', sprintf('%0.1f ',waistEfficiencies*100),'%']);
    logMessage(['averageEfficiencies: ', sprintf('%0.1f ',averageEfficiencies*100),'%']);
    logMessage(['theoretical axial resolution: ', sprintf('%0.0f ',axialResolution*1e9),'nm']);
    logMessage(['numerical axial resolution: ', sprintf('%0.0f ',numericalAxialResolution*1e9),'nm']);
    logMessage(['peak value: ', sprintf('%0.1f ',100*peakValues),'%']);
    logMessage(['relative peak value: ', sprintf('%0.1f ',100*relativePeakValues),'%']);
end
   
function results=plotUsefulPhotonFraction(config,alpha,beta,coreWidths)
    results=struct([]);
    
    result=struct();
    % Calculate the theorectical FOV
    if beta<1
        % Bessel beam
        result.FOV=config.excitation.wavelength/config.sample.refractiveIndex/(2*(1-sqrt(1-(config.excitation.objective.numericalAperture/config.sample.refractiveIndex)^2))*beta);
        axialResolution=config.excitation.wavelength*pi*.05/(config.excitation.objective.numericalAperture*beta);
        axialResolution2PE=config.excitation.wavelength*pi*.05/(config.excitation.objective.numericalAperture*beta);
    elseif alpha~=0
        % Airy
        result.FOV=6*alpha*config.excitation.wavelength/config.sample.refractiveIndex/(1-sqrt(1-(config.excitation.objective.numericalAperture/config.sample.refractiveIndex)^2));
        maxSpFreq=min(0.88,1/(0.05^2*48*alpha))*config.excitation.objective.numericalAperture/(config.excitation.wavelength/2);
        axialResolution=1./maxSpFreq;
    else
        % Gaussian
        result.FOV=4*config.excitation.wavelength*config.sample.refractiveIndex*config.excitation.objective.numericalAperture^-2;
        axialResolution=config.excitation.wavelength/(2*config.excitation.objective.numericalAperture*0.88);
    end
    result.FOV2PE=2*result.FOV;

    result.xRange=single([-75:.05:75]*1e-6);
    result.yRange=result.xRange;
    result.zRange=single([0:.05:1]*result.FOV/2);
    
    config.modulation.alpha=alpha;
    config.modulation.beta=beta; 

    pupilFunctor=@(U,V) (sqrt(U.^2+V.^2)>=(1-config.modulation.beta)).*exp(2i*pi*config.modulation.alpha*(U.^3+V.^3));

    cache=Cache(); % Get the default cache
    key={result.zRange,result.zRange,result.yRange,config,1/sqrt(2),1i/sqrt(2)};
    if (isfield(cache,key))
        logMessage('Loading PSF projection from cache...');
        psf=cache(key);
    else
        psf=calcVectorialPsf(result.xRange,result.yRange,result.zRange,config.excitation.wavelength,@(U,V) pupilFunctor(U,V)/sqrt(2),@(U,V) 1i*pupilFunctor(U,V)/sqrt(2),config.excitation.objective.numericalAperture,config.sample.refractiveIndex,config.excitation.objective.magnification,config.excitation.objective.tubeLength);
        cache(key)=psf;
    end
    
    % Calc PSF projection and OTF slice
    result.psfAtWaist=sum(psf(:,:,1));
    result.fullWidthAtHalfMaximumAtWaist=calcFullWidthAtHalfMaximum(result.xRange,result.psfAtWaist,'Linear');
    result.fullWidthAtHalfMaximumAtWaist2PE=calcFullWidthAtHalfMaximum(result.xRange*2,abs(result.psfAtWaist).^2,'Linear');
    
    cutOffSpatialFrequency=(2*config.excitation.objective.numericalAperture*config.excitation.fractionOfNumericalApertureUsed)/config.excitation.wavelength;
	[XOtf,YOtf,fRel]=calcOtfGridFromSampleFrequencies(1./[diff(result.xRange(1:2)) diff(result.yRange(1:2))],[length(result.xRange) length(result.yRange)],cutOffSpatialFrequency);
    otf=fftshift(fft2(ifftshift(psf(:,:,1))));
    otf2=fftshift(fft2(ifftshift(abs(psf(:,:,1)).^2)));
    otf2=otf2/otf2(1+floor(end/2),1+floor(end/2));
    fRelSel=single(any(fRel.'<=1).')*single(any(fRel<=1));
    fRelSel2=single(any(fRel.'<=2).')*single(any(fRel<=2));
    selSize=[sum(any(fRel.'<=1)) sum(any(fRel<=1))];
    selSize2=[sum(any(fRel.'<=2)) sum(any(fRel<=2))];
    otf=reshape(otf(logical(fRelSel)),selSize);
    otf2=reshape(otf2(logical(fRelSel2)),selSize2);
    XOtf2=0.5*reshape(XOtf(logical(fRelSel2)),selSize2);
    XOtf=reshape(XOtf(logical(fRelSel)),selSize);
    fRel2=0.5*reshape(fRel(logical(fRelSel2)),selSize2);
    fRel=reshape(fRel(logical(fRelSel)),selSize);
    result.otfSlice=conj(otf(1+floor(end/2),1+floor(end/2):-1:1));
    result.otf2Slice=conj(otf2(1+floor(end/2),1+floor(end/2):-1:1));
    result.XOtfSlice=-XOtf(1+floor(end/2),1+floor(end/2):-1:1);
    result.XOtf2Slice=-XOtf2(1+floor(end/2),1+floor(end/2):-1:1);
    result.fRelSlice=fRel(1+floor(end/2),1+floor(end/2):-1:1);
    result.fRel2Slice=fRel2(1+floor(end/2),1+floor(end/2):-1:1);
    
    result.peakValueAtWaist=max(result.psfAtWaist/sum(result.psfAtWaist));
    result.peakValueAtWaist2PE=max(abs(result.psfAtWaist).^2/sum(abs(result.psfAtWaist).^2));
    
    maxSpFreq1PE=result.XOtfSlice(find(abs(result.otfSlice)/max(abs(result.otfSlice))<.05,1,'first')-1);
    result.axialResolution1PE=1/maxSpFreq1PE;
    maxSpFreq2PE=result.XOtf2Slice(find(abs(result.otf2Slice)/max(abs(result.otf2Slice))<.05,1,'first')-1);
    result.axialResolution2PE=1/maxSpFreq2PE;
    
%     figure;
%     plot(result.XOtfSlice*1e-3,abs(result.otfSlice));
%     figure;
%     ssurf(XOtf*1e-3,YOtf*1e-3,abs(result.otf),angle(result.otf));
%     xlim(cutOffSpatialFrequency*[-1 1]*1e-3);
%     ylim(cutOffSpatialFrequency*[-1 1]*1e-3);
%     zlim([0 1]);
 
    [X,Y]=ndgrid(result.xRange,result.yRange);
    for (coreWidthIdx=1:length(coreWidths))
        coreWidth=coreWidths(coreWidthIdx);
       
        result.axialResolution=min(coreWidth,axialResolution);
    
        corePixels=X.^2+Y.^2<=(coreWidth/2)^2;
        corePixels2PE=2*sqrt(X.^2+Y.^2)<=(coreWidth/2);
        focalPlanePixels=abs(X)<=(coreWidth/2);
        focalPlanePixels2PE=2*abs(X)<=(coreWidth/2);
        result.totalPower=[];
        result.corePower=[];
        result.focalPlanePower=[];
        for zIdx=1:length(result.zRange)
            psfSlice=psf(:,:,zIdx);
            result.totalPower(zIdx,1)=sum(psfSlice(:));
            result.totalPower(zIdx,2)=sum(psfSlice(:).^2);
            result.corePower(zIdx,1)=sum(psfSlice(corePixels));
            result.corePower(zIdx,2)=sum(psfSlice(corePixels2PE).^2);
            result.focalPlanePower(zIdx,1)=sum(psfSlice(focalPlanePixels));
            result.focalPlanePower(zIdx,2)=sum(psfSlice(focalPlanePixels2PE).^2);
    %         showImage(cat(3,psfSlice.*corePixels,psfSlice.*~corePixels,psfSlice*0),-1,result.xRange*1e6,result.yRange*1e6);
        end

        % Calculate some derived metrics
        result.powerOutsideCore=result.totalPower-result.corePower;
        result.powerOutsideFocalPlane=result.totalPower-result.focalPlanePower;

        result.powerOutsideCorePerUsefulPower=result.powerOutsideCore./result.corePower;
        result.powerOutsideFocalPlanePerUsefulPower=result.powerOutsideFocalPlane./result.focalPlanePower;
        result.averagePowerOutsideCorePerUsefulPower=sum(result.powerOutsideCore)./sum(result.corePower);
        result.averagePowerOutsideFocalPlanePerUsefulPower=sum(result.powerOutsideFocalPlane)./sum(result.focalPlanePower);

        result.coreEfficiency=result.corePower./result.totalPower;
        result.focalPlaneEfficiency=result.focalPlanePower./result.totalPower;
        result.averageCoreEfficiency=sum(result.corePower)./sum(result.totalPower);
        result.averageFocalPlaneEfficiency=sum(result.focalPlanePower)./sum(result.totalPower);
        
        % Stack all results in an array
        if (isempty(results))
            results=result;
        else
            results(coreWidthIdx)=result;
        end

        if (nargout==0)
            %
            % Display
            %
            fig(coreWidthIdx)=figure();
            subplot(1,2,1);
            plot(zRange*1e6,result.powerOutsideCorePerUsefulPower);
            ylim([0 ceil(1.5*result.powerOutsideCorePerUsefulPower(1)/10)*10]);
            title(sprintf('core of alpha=%0.0f, beta=%0.0f%%',[config.modulation.alpha config.modulation.beta*100]));
            subplot(1,2,2);
            plot(zRange*1e6,result.powerOutsideFocalPlanePerUsefulPower);
            ylim([0 ceil(1.5*result.powerOutsideFocalPlanePerUsefulPower(1)/10)*10]);
            title(sprintf('section of alpha=%0.0f, beta=%0.0f%%',[config.modulation.alpha config.modulation.beta*100]));
        end
    end
    
end

