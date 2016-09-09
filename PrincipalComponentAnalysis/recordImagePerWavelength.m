% recordImagePerWavelength(sampledWavelengths,outputFolder,nbImagesPerWavelength)
%
% Function to record a series of speckle images for various sampledWavelengths
% 
% sampledWavelengths: the range of sampledWavelengths (in meters)
% outputFolder: a string indicating where the data is to be saved. Default:
%               measurement followed by the data and time.
% nbImagesPerWavelength: the number of images to record per wavelength,
%               default: 100
% 
function recordImagePerWavelength(sampledWavelengths,outputFolder,nbImagesPerWavelength)
    if (nargin<1 || isempty(sampledWavelengths))
%         sampledWavelengths=[695:.2:705 706:720]*1e-9;
        sampledWavelengths=[745:805]*1e-9;
        
        sampledWavelengths=unique(round(sampledWavelengths*1e15)*1e-15); %Remove doubles (upto fm)
    end
    if (nargin<2 || isempty(outputFolder))
        outputFolder=strcat('D:\SpeckleWaveMeter\BroadBand\Alumina\measurement_',datestr(now(),'YYYY-mm-dd_HH_MM_ss'));
    end
    if (nargin<3 || isempty(nbImagesPerWavelength))
        nbImagesPerWavelength=100;
    end
    
    logMessage('Measuring for %d wavelengths.',length(sampledWavelengths));
    
    cam=PikeCam();
    cam.regionOfInterest=[0 0 480 320];
    cam.integrationTime=25e-3; % initial integration time only
    cam.gain=100;
    
    mkdir(outputFolder);
    
    images=[];
    for wavelengthIdx=1:length(sampledWavelengths),
        wavelength=sampledWavelengths(wavelengthIdx);
        logMessage('Set laser wavelength to %0.6f nm and hit enter.',wavelength*1e9);
        input('');
        logMessage('ACQUIRING! DON''T MOVE!');
        images=acquireWithoutSaturating(cam,nbImagesPerWavelength);
        images=single(images);
        logMessage('Peak intensity %f, interframe error %f, rel. interframe error %f',[max(images(:)) mean(mean(std(images,0,3))) mean(mean(std(images,0,3)))/mean(images(:))]);
        % Save output
        outputFileName=strcat(outputFolder,sprintf('/imagesForWavelength%0.6fnm.mat',wavelength*1e9));
        recordingTime=datestr(clock(),'dd-mmm-yyyy HH:MM:SS.FFF');
        camGain=cam.gain;
        camIntegrationTime=cam.integrationTime;
        camRegionOfInterest=cam.regionOfInterest;
        save(outputFileName,'sampledWavelengths','images','recordingTime','camGain','camIntegrationTime','camRegionOfInterest'); 
        logMessage('---------------------- %0.1f%% done ------------------------',100*wavelengthIdx/length(sampledWavelengths));
    end
end

function img=acquireWithoutSaturating(cam,nbImagesPerWavelength)
    if (nargin<2 || isempty(nbImagesPerWavelength))
        nbImagesPerWavelength=1;
    end
    
    shortestIntegrationTime=0.01e-3;
    longestIntegrationTime=2000e-3;
    
    cam.background=0;
    img=cam.acquire();
    while ((max(img(:))>=1 && cam.integrationTime>shortestIntegrationTime*2) || (max(img(:))<.5 && cam.integrationTime<longestIntegrationTime/2))
        if (max(img(:))>=1)
            cam.integrationTime=cam.integrationTime/2;
        else
            cam.integrationTime=cam.integrationTime*2;
        end
%         logMessage('Integration time %f ms',cam.integrationTime*1e3);
        img=cam.acquire();
    end
    if (max(img(:))>=1)
        logMessage('Too bright!!!');
    end
    if (max(img(:))<=.5)
        logMessage('Too dark!!!');
    end
    if (nbImagesPerWavelength>1)
        img(:,:,nbImagesPerWavelength)=0;
    end
    for (imgIdx=2:nbImagesPerWavelength)
        img(:,:,imgIdx)=cam.acquire();
    end
    
    img=img/cam.integrationTime;
end