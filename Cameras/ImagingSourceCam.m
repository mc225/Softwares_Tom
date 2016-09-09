%
% Imaging Source camera class
% 
classdef ImagingSourceCam < Cam
    properties (SetAccess = private)
        maxSize;
        minRegionOfInterestSize=[1 1]*2;
        bitsPerPixel;
    end
    properties (Access = private)
        videoInput;
        videoSource;
        canSetIntegrationTime;
        canSetGain;
    end
    properties (Dependent = true)
        integrationTime; % seconds
        gain;
        regionOfInterest;
        defaultNumberOfFramesToAverage;
    end
    
    methods
        function cam=ImagingSourceCam(deviceID)
            cam=cam@Cam();
            if (nargin<1 || isempty(deviceID))
                deviceID=-1;
            end
            
            %Check if there is already a GigE cam configured
            vids=imaqfind('Tag','ImagingSourceVideoinput');
            vid=[];
            vidIdx=0;
            while (vidIdx<length(vids) && isempty(vid))
                vidIdx=vidIdx+1;
                vid=vids{vidIdx};
                if (deviceID>0 && vid.DeviceID~=deviceID)
                    vid=[];
                end
            end
            if (isempty(vid))
                hwInfo=imaqhwinfo('winvideo');
                deviceInfoIdx=find(cell2mat(hwInfo.DeviceIDs)==deviceID);
                if (~isempty(deviceInfoIdx))
                    formats=hwInfo.DeviceInfo(deviceInfoIdx(1)).SupportedFormats;
                    if (any(strcmp(formats,'Y800_1280x960')))
                        vid=videoinput('winvideo',deviceID,'Y800_1280x960');
                        vid.ReturnedColorspace = 'bayer';
                        vid.BayerSensorAlignment = 'grbg';
                    elseif (any(strcmp(formats,'YUY2_640x480')))
                        vid=videoinput('winvideo',deviceID,'YUY2_640x480');
                        set(vid,'ReturnedColorSpace','YCbCr');
                    elseif (any(strcmp(formats,'YUY2_1280x1024')))
                        vid=videoinput('winvideo',deviceID,'YUY2_1280x1024');
                        set(vid,'ReturnedColorSpace','YCbCr');
                    elseif (any(strcmp(formats,'YUY2_1024x768')))
                        vid=videoinput('winvideo',deviceID,'YUY2_1024x768');
                        set(vid,'ReturnedColorSpace','YCbCr');
                    else
                        devInfo=imaqhwinfo('winvideo',deviceID);
                        vid=videoinput('winvideo',deviceID,devInfo.DefaultFormat);
                        set(vid,'ReturnedColorSpace','YCbCr');
                    end
                    vid.Tag='ImagingSourceVideoinput';
                else
                     error('No winvideo camera found.');
                end
            end
            %Retrieve the number of bits per pixel from the format string
            cam.bitsPerPixel=8;
            
            src=getselectedsource(vid);
            
            cam.videoInput=vid;
            cam.videoSource=src;
            
            ensureVideoStopped(cam);
            
            cam.videoInput.ROIPosition=[0 0 cam.videoInput.VideoResolution];
            exposureTime=-4; % [-13 +3] for Imaging Source camera
            try
                cam.videoSource.Exposure=exposureTime;
                
                cam.canSetIntegrationTime=true;
            catch Exc
                logMessage('Could not set exposure, most likely this is not an ImagingSource device or driver.');
                cam.canSetIntegrationTime=false;
            end
            try
                cam.videoSource.GainMode = 'manual';
                cam.videoSource.Gain=720; %260-1023;
                cam.canSetGain=true;
            catch Exc
                logMessage('Could not set gain, most likely this is not an ImagingSource device or driver.');
                cam.canSetGain=false;
            end
            try
                cam.videoSource.ExposureMode = 'manual';
            catch Exc
            end
            try
                cam.videoSource.Gamma=100;
            catch Exc
            end
            try
                frameRateInfo=propinfo(src,'FrameRate');
                frameRateOptions=frameRateInfo.ConstraintValue;
                frameRateOptionsNumeric=str2num(cell2mat(frameRateOptions(:))).';
                [ign maxI]=max(frameRateOptionsNumeric);
                src.FrameRate = frameRateOptions{maxI};
            catch Exc
                logMessage('Could not set framerate to 60Hz, keeping it at %d',src.FrameRate);
            end
                
            set(cam.videoInput,'TriggerRepeat',Inf);
            triggerconfig(cam.videoInput,'manual');
            ensureVideoStarted(cam);
            
            cam.maxSize([2 1])=cam.videoInput.VideoResolution;
            
            cam.regionOfInterest=[0 0 cam.maxSize];
            
            cam.defaultNumberOfFramesToAverage=1;
        end
        function cam=set.integrationTime(cam,newIntegrationTime)
            if (cam.canSetIntegrationTime)
                ensureVideoStopped(cam);
                exposureValue=round(log2(newIntegrationTime));
                propInfo=propinfo(cam.videoSource,'Exposure');
                constraints=propInfo.ConstraintValue;
                exposureValue=max(min(exposureValue,constraints(2)),constraints(1));
                cam.videoSource.Exposure=exposureValue;
                ensureVideoStarted(cam);
            end
        end
        function integrationTime=get.integrationTime(cam)
            if (cam.canSetIntegrationTime)
                integrationTime=2.^(cam.videoSource.Exposure);
            else
                integrationTime=0;
            end
        end
        function cam=set.gain(cam,newGain)
            if (cam.canSetGain)
                ensureVideoStopped(cam);
                registerValue=log10(newGain)*1023/1.8; % 10 bits value linear with gain in dB with maximum of 18dB it seems
                % 36dB accorduing to http://www.theimagingsourceforums.com/archive/index.php/t-324235.html
                cam.videoSource.Gain=max(260,min(1023,round(registerValue))); %260-1023;
                ensureVideoStarted(cam);
            end
        end
        function gain=get.gain(cam)
            if (cam.canSetGain)
                gain=10^(1.8*cam.videoSource.Gain/1023);
            else
                gain=0;
            end
        end
        function cam=set.defaultNumberOfFramesToAverage(cam,newDefaultNumberOfFramesToAverage)
            ensureVideoStopped(cam);
            cam.videoInput.FramesPerTrigger=newDefaultNumberOfFramesToAverage;
            ensureVideoStarted(cam);
        end
        function defaultNumberOfFramesToAverage=get.defaultNumberOfFramesToAverage(cam)
            defaultNumberOfFramesToAverage=cam.videoInput.FramesPerTrigger;
        end
        function cam=set.regionOfInterest(cam,newROI)
            cam.background=0;
            
            if (isempty(newROI))
                newROI=[0 0 cam.maxSize];
            end
            
            ensureVideoStopped(cam);
            %Make sure that the size is at least 128x128
            sizeIncrease=max(0,cam.minRegionOfInterestSize-newROI(3:4));
            newROI(1:2)=newROI(1:2)+floor(-sizeIncrease/2);
            newROI(3:4)=newROI(3:4)+sizeIncrease;
            %Round to even numbers when using a Bayer filter
            newROI=2*floor(newROI./2);
            %Clip the region of interest to the active area of the CCD
            newROI(1:2)=min(max(0,newROI(1:2)),cam.maxSize-1);
            newROI(3:4)=min(newROI(1:2)+max(2,newROI(3:4)),cam.maxSize)-newROI(1:2);
            cam.videoInput.ROIPosition=newROI([2 1 4 3]);
            ensureVideoStarted(cam);
        end
        function roi=get.regionOfInterest(cam)
            roi([2 1 4 3])=cam.videoInput.ROIPosition;
        end
        function cam=acquireBackground(cam,nbFrames)
            if (nargin<2)
                nbFrames=cam.defaultNumberOfFramesToAverage;
            end
            cam=acquireBackground@Cam(cam,nbFrames);
        end
        function img=acquire(cam,nbFrames)
            if (nargin<2)
                nbFrames=cam.defaultNumberOfFramesToAverage;
            end
            img=acquire@Cam(cam,nbFrames);
        end
        function delete(cam)
            delete@Cam(cam);
            
            ensureVideoStopped(cam);
            delete(cam.videoInput);
            clear cam.videoInput;
        end
    end
    methods(Access = protected)
        function img=acquireDirect(cam,nbFrames)
            framesPerTrigger=cam.defaultNumberOfFramesToAverage;
             
            %Round the number of frames up to a multiple of the frames per trigger
            nbFrames=nbFrames*ceil(nbFrames/framesPerTrigger);
            for (triggerIdx=1:nbFrames/framesPerTrigger)
                snapshots=[];
                while (isempty(snapshots))
                    ensureVideoStarted(cam);
                    trigger(cam.videoInput);
                    try
                        snapshots=double(getdata(cam.videoInput,framesPerTrigger));
                    catch Exc
                        logMessage('Something went wrong, reseting the videoinput object...');
                        stop(cam.videoInput);
                        start(cam.videoInput);
                    end
                end
                %Debug info
                if (any(snapshots(:)==(2^cam.bitsPerPixel-1)))
                    logMessage('%u pixels are saturated in %u frames!',[sum(snapshots(:)==(2^cam.bitsPerPixel-1)) framesPerTrigger]);
                end
                %Normalize the maximum graylevel to 1
                snapshots=snapshots./(2^cam.bitsPerPixel-1);
                
                snapshots=snapshots(:,:,1);
                snapshots=sum(snapshots,4);
                
                if (triggerIdx==1)
                    img=snapshots;
                else
                    img=img+snapshots;
                end

            end
            
            if (nbFrames>1)
                img=img./nbFrames;
            end
            
        end
    end
    methods (Access = private)
        function ensureVideoStarted(cam)
            if (~strcmpi(cam.videoInput.Running,'on'))
                start(cam.videoInput);
            end
        end
        function ensureVideoStopped(cam)
            if (strcmpi(cam.videoInput.Running,'on'))
                stop(cam.videoInput);
            end
        end
    end
    
end


