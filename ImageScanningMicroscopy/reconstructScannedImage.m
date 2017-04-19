% [imageScanningImage xRangeOutputImg yRangeOutputImg]=reconstructScannedImage(inputImages,pixelPitch,magnification,outputFileName)
%
function [imageScanningImage xRangeOutputImg yRangeOutputImg]=reconstructScannedImage(inputImages,samplePositions,pixelPitch,magnification,outputFileName)
    close all; % TODO: Remove this line when done testing
    
    if (nargin<3 || isempty(pixelPitch))
        pixelPitch=[1 1]*7.4e-6;
    end
    if (length(pixelPitch)==1)
        pixelPitch(2)=pixelPitch;
    end
    if (nargin<4 || isempty(magnification))
        magnification=100;
    end
    pixelPitchInSample=pixelPitch/magnification;
    pixelPitchForOutput=pixelPitchInSample/2;
    regionOfInterest=[];
    if (nargin<1 || isempty(inputImages))
        inputImages='C:\Users\Tom\Dropbox\ScansChris\12 aug scans\scan3\scanselectedROI';
        regionOfInterest=[]; %[219-16 379-16 32 32];
        
%         % Generate synthetic data
%         imgSize=[32 32]/2; % Pixel size of recorded image
%         rangeX=[-0.5:.05:0.5]*1e-6;
%         rangeY=[-0.5:.05:0.5]*1e-6;
%         [stagePositionsX stagePositionsY]=ndgrid(rangeX,rangeY);
%         clear rangeX rangeY;
%         
%         [inputImages samplePositions]=simulateImageScanningMicroscopy(imgSize,stagePositionsX,stagePositionsY,pixelPitchInSample);
%         clear imgSize stagePositionsX stagePositionsY;

%         % Test output
%         inputImages4D=reshape(inputImages,[size(inputImages,1) size(inputImages,2) [1 1]*sqrt(size(inputImages,3))]);
%         for (xIdx=1:size(inputImages4D,3))
%             for (yIdx=1:size(inputImages4D,4))
%                 tmp([1:size(inputImages4D,1)]+size(inputImages4D,1)*(xIdx-1),[1:size(inputImages4D,2)]+size(inputImages4D,2)*(yIdx-1))=inputImages4D(:,:,xIdx,yIdx);
%             end
%         end
%         imagesc(tmp);
    end
    if (nargin==1 || (nargin>0 && isempty(samplePositions)))
        nbImagesXY=[floor(sqrt(nbImages)) ceil(nbImages/floor(sqrt(nbImages)))];
        stepSize=pixelPitch./magnification;
        rangeX=([1:nbImagesXY(1)]-floor(nbImagesXY(1)/2))*stepSize(1);
        rangeY=([1:nbImagesXY(2)]-floor(nbImagesXY(2)/2))*stepSize(2);
        [stagePositionsX stagePositionsY]=ndgrid(rangeX,rangeY);
        clear nbSamples nbImagesXY stepSize rangeX rangeY;
        samplePositions=[stagePositionsX(:) stagePositionsY(:)];
    end
    if (ischar(inputImages))
        if (isdir(inputImages))
            logMessage('Loading input images from folder %s',inputImages);
            [inputImages samplePositions]=loadInputImages(inputImages,regionOfInterest);
        else
            logMessage('%s is not a folder.',inputImages);
            return;
        end
        
    end
    if (nargin<5)
        outputFileName='recontructedImage.fig';
    end
    
    % For deconvolution
    resolutionIllumination=400e-9;
    resolutionDetection=100e-9;
    
    % Process input
    inputSize=[size(inputImages,1) size(inputImages,2)];
    nbSamples=size(inputImages,3);
    
    % Create output grid and empty output image
    xRangeInputImg=([1:inputSize(1)]-floor(inputSize(1)/2))*pixelPitchInSample(1);
    yRangeInputImg=([1:inputSize(2)]-floor(inputSize(2)/2))*pixelPitchInSample(2);
    [imgInX imgInY]=ndgrid(xRangeInputImg,yRangeInputImg);
    
    xRangeOutputImg=[xRangeInputImg(1)+min(samplePositions(:,1)):pixelPitchForOutput(1):xRangeInputImg(end)+max(samplePositions(:,1))];
    yRangeOutputImg=[yRangeInputImg(2)+min(samplePositions(:,2)):pixelPitchForOutput(2):yRangeInputImg(end)+max(samplePositions(:,2))];
    [imgOutX imgOutY]=ndgrid(xRangeOutputImg,yRangeOutputImg);
        
    %Calculate background
    backgroundImage=min(inputImages,[],3);
    inputImages=inputImages-repmat(backgroundImage,[1 1 nbSamples]);
    totalIntensities=squeeze(sum(sum(inputImages)));
    % Confocal pinhole
    pinholeRadius=3; % pixels
    confocalPinHole=imgInX.^2+imgInY.^2<=pinholeRadius^2;
    confocalIntensities=squeeze(sum(sum(inputImages.*repmat(confocalPinHole,[1 1 nbSamples]))));
    
    %
    % Reconstruct
    %
    interpolantLaserScanning=TriScatteredInterp(samplePositions(:,1),samplePositions(:,2),totalIntensities,'linear');
    interpolatedLaserScanningImage = interpolantLaserScanning(imgOutX,imgOutY);
    interpolatedLaserScanningImage(isnan(interpolatedLaserScanningImage))=0;
    interpolantConfocal=TriScatteredInterp(samplePositions(:,1),samplePositions(:,2),confocalIntensities,'linear');
    interpolatedConfocalImage = interpolantConfocal(imgOutX,imgOutY);
    interpolatedConfocalImage(isnan(interpolatedConfocalImage))=0;
    
    %Image scanning reconstruction
    imageScanningImage=0*imgOutX;
    for (imgIdx=1:nbSamples)
        fractionOfReconstruction=imgIdx/nbSamples
        
        samplePos=samplePositions(imgIdx,:);
        
        inputImageCenter=1+floor(inputSize./2);
%         inputImage=zeros(inputSize);
%         inputImage(inputImageCenter(1),inputImageCenter(2))=confocalIntensities(imgIdx); %totalIntensities(imgIdx);
%         newImage=interpn(imgInX+samplePos(1),imgInY+samplePos(2),inputImage,imgOutX,imgOutY,'*nearest',0); % show as squares
%         newImage(newImage<max(newImage(:)))=0;

%         newImage=interpn(imgInX+samplePos(1),imgInY+samplePos(2),inputImage,imgOutX,imgOutY,'*nearest',0);
        
        % Image Scanning Microscopy reconstruction
        newImage=interpn(imgInX./2+samplePos(1),imgInY./2+samplePos(2),inputImages(:,:,imgIdx),imgOutX,imgOutY,'*cubic',0);
        
%         imageScanningImage=imageScanningImage+newImage;
        imageScanningImage=max(imageScanningImage,newImage);
    end
    
    % Deconvolve
    cutOffSpatialFrequency=1/max(resolutionDetection,resolutionIllumination); %Only for the regularization
    [XOtf,YOtf,fRel]=calcOtfGridFromSampleFrequencies(1./[diff(imgOutX(1:2,1)) diff(imgOutY(1,1:2))],size(imgOutX)*2,cutOffSpatialFrequency);
    noiseToSignalLevel=100;
%     noiseToSignalPowerRatio=(fRel./noiseToSignalLevel).^2;
    imageScanningImagePadded=imageScanningImage([1:end, end*ones(1,ceil(end/2)), ones(1,floor(end/2))],[1:end, end*ones(1,ceil(end/2)), ones(1,floor(end/2))]);
    confocalImagePadded=interpolatedConfocalImage([1:end, end*ones(1,ceil(end/2)), ones(1,floor(end/2))],[1:end, end*ones(1,ceil(end/2)), ones(1,floor(end/2))]);
    otfDoubleIllum=max(0,(1-0.5*sqrt(XOtf.^2+YOtf.^2).*resolutionIllumination));
    otfDoubleDetect=max(0,(1-0.5*sqrt(XOtf.^2+YOtf.^2).*resolutionDetection));
    opticalTransferFunction=otfDoubleIllum.*otfDoubleDetect;
    filter=calcWienerFilter(opticalTransferFunction,noiseToSignalLevel^2,otfDoubleDetect);
    imageScanningDeconvolvedImage=ifft2(fft2(imageScanningImagePadded).*ifftshift(filter),'symmetric');
    confocalDeconvolvedImage=ifft2(fft2(confocalImagePadded).*ifftshift(filter),'symmetric');
    clear imageScanningImagePadded filter;
    imageScanningDeconvolvedImage=imageScanningDeconvolvedImage(1:end/2,1:end/2); % Crop padding again
    confocalDeconvolvedImage=confocalDeconvolvedImage(1:end/2,1:end/2); % Crop padding again
    
    % Make first and second coordinate Y and X for consistency with
    % functions such as imagesc
    tmp=xRangeOutputImg;
    xRangeOutputImg=yRangeOutputImg;
    yRangeOutputImg=tmp;
    clear tmp;
    
    %
    % Output
    %
    
    %If no output specified, display result
    if (nargout<1 || ~isempty(outputFileName))
        fig=figure();
        axs(1)=subplot(2,3,1);
        scatter(samplePositions(:,2)*1e6,samplePositions(:,1)*1e6,10,confocalIntensities,'s','filled'); title('confocal blocks');
        set(gca,'Color',[0 0 .75]);
        axis image; xlabel('x [\mum]'); ylabel('y [\mum]');
        
        axs(2)=subplot(2,3,2);
        imagesc(xRangeOutputImg*1e6,yRangeOutputImg*1e6,interpolatedLaserScanningImage); title('laser scan');
        axis image; colormap(hot); colorbar(); xlabel('x [\mum]'); ylabel('y [\mum]');
        
        axs(3)=subplot(2,3,3);
        imagesc(xRangeOutputImg*1e6,yRangeOutputImg*1e6,interpolatedConfocalImage); title('confocal');
        axis image; colormap(hot); colorbar(); xlabel('x [\mum]'); ylabel('y [\mum]');
        
        axs(4)=subplot(2,3,4);
        imagesc(xRangeOutputImg*1e6,yRangeOutputImg*1e6,imageScanningImage); title('image scan');
%         hold on;
%         scatter(samplePositions(:,2)*1e6,samplePositions(:,1)*1e6,10,[0 0 .75],'+');
%         hold off;
        axis image; colormap(hot); colorbar(); xlabel('x [\mum]'); ylabel('y [\mum]');
        
        axs(5)=subplot(2,3,5);
        imagesc(xRangeOutputImg*1e6,yRangeOutputImg*1e6,confocalDeconvolvedImage); title('confocal deconvolved');
        axis image; colormap(hot); colorbar(); xlabel('x [\mum]'); ylabel('y [\mum]');
        
        axs(6)=subplot(2,3,6);
        imagesc(xRangeOutputImg*1e6,yRangeOutputImg*1e6,imageScanningDeconvolvedImage); title('image scan deconvolved');
        axis image; colormap(hot); colorbar(); xlabel('x [\mum]'); ylabel('y [\mum]');
        
        linkaxes(axs);
        
        % Save result
        if (~isempty(outputFileName))
            saveas(fig,outputFileName);
        end
    
        clear imageScanningImage; % Avoid cluttering the command line
    end
    
end

function [recordedImages stagePositions]=simulateImageScanningMicroscopy(imgSize,samplePositionsX,samplePositionsY,pixelPitchInSample)
    stagePositions=[samplePositionsX(:),samplePositionsY(:)];
    
%     sample=getSpokeTarget(256,256,12,1);
%     sample(1:end/2,1:end/2)=0;
%     sample(1:end/2,end/2:end)=sample(1:end/2,end/2:end)/2;
%     sample(end/2:end,1:end/2)=sample(end/2:end,1:end/2)/4;
    %sample=getTestImage('usaf_512x512.png'); sample=sample(1:4:end,1:4:end);
    sample=zeros(64);
%     sample(1+floor(end/2),1+floor(end/2))=1;
    sample(1+floor(end/2)-8+[1:3],1+floor(end/2)+[1:3])=1;
    sample(1+floor(end/2)+[1:3],1+floor(end/2)+[1:3]+8)=1;
    samplePixelPitch=[1 1]*25e-9;
    
    xRangeSample=([1:size(sample,1)]-floor(size(sample,1)/2))*samplePixelPitch(1);
    yRangeSample=([1:size(sample,2)]-floor(size(sample,2)/2))*samplePixelPitch(2);
    figure; imagesc(yRangeSample*1e6,xRangeSample*1e6,sample); axis image; xlabel('x [\mum]'); ylabel('y [\mum]');
    [sampleX sampleY]=ndgrid(xRangeSample,yRangeSample);
    clear xRangeSample yRangeSample;
    
    xRangeImg=([1:imgSize(1)*2]-1-imgSize(1))*pixelPitchInSample(1);
    yRangeImg=([1:imgSize(2)*2]-1-imgSize(2))*pixelPitchInSample(2);
    [imgX imgY]=ndgrid(xRangeImg,yRangeImg);
    sigma=100e-9;
    illumination=exp(-(imgX.^2+imgY.^2)./(2*sigma^2));
    otf=fft2(ifftshift(illumination));
    otf=otf./otf(1);
    
    for posIdx=1:size(stagePositions,1),
        fractionOfInput=posIdx/size(stagePositions,1)
        stagePos=stagePositions(posIdx,:);
        %shift
        illuminatedSample=interpn(sampleX-stagePos(1),sampleY-stagePos(2),sample,imgX,imgY,'*cubic',0);
        %simulate illumination
        illuminatedSample=illuminatedSample.*illumination;
        % simulate detection
        img=ifft2(fft2(illuminatedSample).*otf,'symmetric');
        % Crop padding again
        img=circshift(img,-ceil(size(img)./4));
        img=img(1:end/2,1:end/2);
        %Add noise
        img=img+0.001*randn(size(img));
        % Store in two large matrices
        if (posIdx>1)
            recordedImages(:,:,posIdx)=img;
        else
            recordedImages=img;
        end
%         showImage(img*100);
%         drawnow();
%         pause(.02);
    end
end

function [recordedImages samplePositions]=loadInputImages(inputImageFolder,regionOfInterest)
    files=dir([inputImageFolder,'\*.png']);
    posIdx=1;
    for fileIdx=1:length(files),
        fractionOfInput=fileIdx/length(files)
        fileName=files(fileIdx).name;
        samplePos=regexpi(fileName,'x(\d+)y(\d+)\.png','tokens');
        if (~isempty(samplePos))
            try
                img=imread([inputImageFolder,'/',fileName]);
                img=double(img)./255; % normalize to dynamic range
                img=mean(img,3); % Convert to gray scale
                % Reduce the ROI and create a sub folder if requested
                if (~isempty(regionOfInterest))
                    img=img(regionOfInterest(1)+[1:regionOfInterest(3)]-1,regionOfInterest(2)+[1:regionOfInterest(4)]-1);
                    selectedROIFolder=[inputImageFolder,'selectedROI'];
                    if (~exist(selectedROIFolder,'dir'))
                        mkdir(selectedROIFolder);
                    end
                    imwrite(img,[selectedROIFolder,'/',fileName]);
                end
                samplePos=samplePos{1};
                samplePos=str2double(samplePos)*1e-9/10; % 1/10 nm units
                if (~all(samplePos==0)) % Ignorethe wide field image
                    % Store in two large matrices
                    if (posIdx>1)
                        recordedImages(:,:,posIdx)=img;
                        samplePositions(posIdx,:)=samplePos;
                    else
                        recordedImages=img;
                        samplePositions=samplePos(:).';
                    end
                    posIdx=posIdx+1;
                else
                    logMessage('Skipping wide field image at coordinate (0,0).');
                end
            catch Exc
                logMessage('Error reading file %s:',fileName);
                logMessage(Exc.message)
            end
        else
            logMessage('Skipping file %s, unrecognized file name.',fileName);
        end
    end
end