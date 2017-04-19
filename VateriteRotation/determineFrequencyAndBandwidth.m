% [frequency bandwidth]=determineFrequencyAndBandWidth(signal,samplesPerSecond)
%
%
function [frequency bandwidth]=determineFrequencyAndBandWidth(signal,samplesPerSecond)
    close all;
    if (nargin<1)
        data=load('C:\Users\Tom\Documents\labsoftware\matlab\VateriteRotation\s.10155.mat');
        signalA=double(data.data16a)*2^-15;
        signalB=double(data.data16b)*2^-15; % Why has the b-signal more power in the higher orders than in the first order?
        signal=signalA; %./(signalA+signalB);
        clear data;
    end
    if (nargin<2)
        samplesPerSecond=2^14;
    end
    
    firstOrderFrequencyInitialEstimate=539;
    
    signal=signal(1:floor(4*end/16));
    
    nbSamples=length(signal);
    samplePeriod=1/samplesPerSecond;
    recordingTime=nbSamples*samplePeriod;
    
    blockPeriod=16; % in seconds
    blockStepPeriod=1; % in seconds
    
    subSampling=1;
    
    times=[0:samplePeriod:recordingTime-samplePeriod];
    
%     [Qs,params,resNorms]=determineProfile(times,signal,firstOrderFrequencyInitialEstimate,1);
    [Qs,params,resNorms]=determineProfile(times,signal,firstOrderFrequencyInitialEstimate,1);
    Qs
    
%     % Normalize signal
%     lowPass=bandPass(times,signal.^2,[0 firstOrderFrequencyInitialEstimate/10]);
%     signal=signal-lowPass;
%     envelope=sqrt(bandPass(times,signal.^2,[0 firstOrderFrequencyInitialEstimate/10]));
%     signal=signal./envelope;
    
    order=1;
    singleOrderOfSignal=bandPass(times,signal,order*firstOrderFrequencyInitialEstimate+[-1 1]*firstOrderFrequencyInitialEstimate/2);
        
    nbSelectedSamples=floor(blockPeriod*samplesPerSecond);
    nbSelectedSubSamples=nbSelectedSamples*subSampling;
    blockFrequencies=([1:nbSelectedSubSamples]-1-floor(nbSelectedSubSamples/2))*samplesPerSecond/nbSelectedSubSamples;
    blockStartTimes=[0:blockStepPeriod:(recordingTime-blockPeriod)];
    dominantFrequencies=zeros(1,length(blockStartTimes));
    spectra=zeros(length(blockStartTimes),nbSelectedSubSamples);
    for blockIdx=1:length(blockStartTimes)
        blockStartTime=blockStartTimes(blockIdx);
        [~, firstSample]=min(abs(times-blockStartTime));
        selectedIndexes=firstSample-1+[1:nbSelectedSamples];
        signalSelection=singleOrderOfSignal(selectedIndexes);
        timesSelection=times(selectedIndexes);
        window=hann(timesSelection,blockStartTime+blockPeriod/2,blockPeriod);
        windowedSignal=signalSelection.*window;
        if (subSampling>1)
            windowedSignal(end*subSampling)=0; % zero pad time domain, sinc interpolate spectrum
        end
        spectra(blockIdx,:)=fftshift(fft(windowedSignal));
        [~,dominantFrequencyIndex]=max(abs(spectra(blockIdx,:)));
        dominantFrequency=abs(blockFrequencies(dominantFrequencyIndex));
        dominantFrequencies(blockIdx)=dominantFrequency;
        if (floor(blockStartTime)>floor(blockStartTime-blockStepPeriod))
            logMessage('Processed %0.0f seconds',blockStartTime);
        end
    end
    dominantFrequencies=dominantFrequencies/order;
    
    figure;
    subplot(2,2,1);
    plot(repmat(blockFrequencies.'*1e-3,[1 length(blockStartTimes)]),abs(spectra.'));
    xlim([0 blockFrequencies(end)]*1e-3);
    subplot(2,2,2);
    plot([0:blockStepPeriod:(recordingTime-blockPeriod)],dominantFrequencies);
    subplot(2,2,[3 4]);
    plot(times,singleOrderOfSignal);
    xlim([0 10/firstOrderFrequencyInitialEstimate]);
end

function window=hann(t,center,width)
    window=(abs(t-center)<width/2).*(0.5*(1+cos(2*pi*(t-center)/width)));
end
function window=hamm(t,center,width)
    alpha=0.54;
    window=(abs(t-center)<width/2).*(alpha+(1-alpha)*cos(2*pi*(t-center)/width));
end

function filteredSignal=bandPass(times,signal,band)
    if (length(times)>1)
        samplesPerSecond=1./diff(times(1:2));
    else
        samplesPerSecond=times;
    end
    
    nbSamples=length(signal);
    
    spectrum=fftshift(fft(signal));
    frequencies=([1:nbSamples]-1-floor(nbSamples/2))*samplesPerSecond/nbSamples;
    
    filter=abs(frequencies)>=band(1) & abs(frequencies)<band(2);
    
    filteredSignal=ifft(ifftshift(spectrum.*filter));
end

function [Qs params resNorms]=determineProfile(times,signal,firstOrderFrequencyInitialEstimate,maxNbOfOrders)
    if (length(times)>1)
        samplesPerSecond=1./diff(times(1:2));
    else
        samplesPerSecond=times;
    end
    nbSamples=length(signal);
    frequencies=([1:nbSamples]-1-floor(nbSamples/2))*samplesPerSecond/nbSamples;
    if (nargin<4)
        maxNbOfOrders=ceil(frequencies(end)/firstOrderFrequencyInitialEstimate);
    end
    
    dFreq=diff(frequencies(1:2));
    
    spectrum=fftshift(fft(signal));
    
    centralFrequencies=[1:maxNbOfOrders]*firstOrderFrequencyInitialEstimate;
    centralFrequencies=centralFrequencies(centralFrequencies<=frequencies(end)-0.5*firstOrderFrequencyInitialEstimate);
    params=zeros(length(centralFrequencies),5);
    resNorms=zeros(1,length(centralFrequencies));
    Qs=zeros(1,length(centralFrequencies));
    estimatedCentralFrequencies=zeros(1,length(centralFrequencies));
    gammas=zeros(1,length(centralFrequencies));
    sigmas=zeros(1,length(centralFrequencies));
    for orderIdx=1:length(centralFrequencies)
        centralFrequency=centralFrequencies(orderIdx);
        band=centralFrequency+[-0.5 0.5]*firstOrderFrequencyInitialEstimate;
        freqSel=frequencies>=band(1) & frequencies<band(2);
        absSpectrumSel=abs(spectrum(freqSel));
        quartileDiffInSampleUnits = diff(quantile(absSpectrumSel,[.25 .75]));
%         gamma=0.5*diff(frequencies(1:2))*quartileDiffInSampleUnits;
        gamma= 1/pi/max(absSpectrumSel);
        params0=[dFreq 0 centralFrequency sum(absSpectrumSel)-median(absSpectrumSel)*numel(absSpectrumSel) median(absSpectrumSel)];
        bounds= [0 0          band(1) min(absSpectrumSel)*length(absSpectrumSel) min(absSpectrumSel);  Inf Inf band(2) 2*sum(absSpectrumSel) max(absSpectrumSel)];
%         [yPrime, params(orderIdx,:), resNorms(orderIdx)] = voigtFit(frequencies(freqSel),absSpectrumSel,params0,bounds);
        [yPrime, params(orderIdx,:), resNorms(orderIdx)] = voigtTranslateScale(frequencies(freqSel),absSpectrumSel,params0,bounds);
%         [yPrime, params(orderIdx,:), resNorms(orderIdx)] = lorentzFit(frequencies(freqSel),absSpectrumSel,params0,bounds);
        figure(1);
        plot(frequencies(freqSel),absSpectrumSel); hold on; plot(frequencies(freqSel),yPrime,'r'); hold on; xlim(540+[-1 1]*20);
        message=sprintf('Res. Norm = %0.3f',sqrt(resNorms/numel(absSpectrumSel)));
        logMessage(message);
        title(message);
        
        %A=f0/Q
        %H(f)=A.*f./(f.^2+A.*f+f0^2)=A.*f/((f-f0).^2+(2*f0+A).*f)
        %H(f)=P1./((f - P2).^2 + P3)=P1./(f.^2   -2 P2 f    +P2.^2 +P3)
        estimatedCentralFrequencies(orderIdx)=params(orderIdx,3);
        gammas(orderIdx)=params(orderIdx,1);
        lineWidths=2*gammas;
        sigmas(orderIdx)=params(orderIdx,2);
%         A=params(orderIdx,3)-2*estimatedCentralFrequencies(orderIdx);
%         Qs(orderIdx)=estimatedCentralFrequencies(orderIdx)/A;
        Qs(orderIdx)=estimatedCentralFrequencies(orderIdx)./lineWidths(orderIdx);
    end
    hold off;
    
end

%     params=[gamma,sigma,x0l,integratedValue,offset];
%     bounds=[lowerParam; upperParam]; % where lower/upperParam are row vectors as params
function [fittedY, params, resNorm] = lorentzFit(X,Y,params0,bounds)
    fitoptions=optimset('Display','final');
    fitFunctor=@(p,X) voigt(X,p(1),params0(2),p(3),p(2),0)+p(4);
    [params,resNorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(fitFunctor,params0([1 3:end]),X,Y,bounds(1,[1 3:end]),bounds(2,[1 3:end]),fitoptions);
%     [params,resNorm,exitflag,output] = fminsearch(@(p) sum(abs(fitFunctor(p,X)-Y).^2),params0,fitoptions);
    fittedY=fitFunctor(params,X);
%     %
%     plot(X,Y,X,fittedY); xlim([500 580]); title(sqrt(resNorm./length(X)))

    params=[params(1) params0(2) params(2:end)];
end


%     params=[gamma,sigma,x0l,integratedValue,offset];
%     bounds=[lowerParam; upperParam]; % where lower/upperParam are row vectors as params
function [fittedY, params, resNorm] = voigtTranslateScale(X,Y,params0,bounds)
    fitoptions=optimset('Display','final');
    fitFunctor=@(p,X) voigt(X,params0(1),params0(2),p(1),p(2),0)+p(3);
    [params,resNorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(fitFunctor,params0(3:end),X,Y,bounds(1,3:end),bounds(2,3:end),fitoptions);
    fittedY=fitFunctor(params,X);
%     %
%     plot(X,Y,X,fittedY); xlim([500 580]); title(sqrt(resNorm./length(X)))
    
    params=[params0(1:2) params];
end

%     params=[gamma,sigma,x0l,integratedValue,offset];
%     bounds=[lowerParam; upperParam]; % where lower/upperParam are row vectors as params
function [fittedY, params, resNorm] = voigtFit(X,Y,params0,bounds)
    fitoptions=optimset('Display','final');
    fitFunctor=@(p,X) voigt(X,p(1),p(2),p(4),p(3),0)+p(5);
    [params,resNorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(fitFunctor,params0,X,Y,bounds(1,:),bounds(2,:),fitoptions);
%     [params,resNorm,exitflag,output] = fminsearch(@(p) sum(abs(fitFunctor(p,X)-Y).^2),params0,fitoptions);
    fittedY=fitFunctor(params,X);
%     %
%     plot(X,Y,X,fittedY); xlim([500 580]); title(sqrt(resNorm./length(X)))
end
