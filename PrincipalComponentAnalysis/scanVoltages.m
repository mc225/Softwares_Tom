% scanVoltages(voltages,outputFolder,nbImagesPerWavelength,laserStabilizationTime)
%
% Function to record a series of speckle images for various wavelengths
% automatically.
% 
% voltages: the range of voltages to drive the Sacher laser
% outputFolder: a string indicating where the data is to be saved. Default:
%               D:\SpeckleWaveMeter\NarrowBand\FabryPerot\withoutDiffuserAnd2F\measurement followed by the data and time.
% nbImagesPerWavelength: the number of images to record per wavelength,
%               default: 50
% laserStabilizationTime: the number of seconds to wait after changing the wavelength
% 
function scanVoltages(voltages,outputFolder,nbImagesPerWavelength,laserStabilizationTime)
    if (nargin<1 || isempty(voltages))
        %voltages=[-5:0.5:-1 -0.9:0.1:0.9 1:.5:5];
        voltages=[-5:0.1:5];
    end
    if (nargin<2 || isempty(outputFolder))
        outputFolder=strcat('D:\SpeckleWaveMeter\NarrowBand\FabryPerot\withoutDiffuserAnd2F\measurement_',datestr(now(),'YYYY-mm-dd_HH_MM_ss'));
    end
    if (nargin<3 || isempty(nbImagesPerWavelength))
        nbImagesPerWavelength=100;
    end
    if (nargin<4 || isempty(laserStabilizationTime))
        laserStabilizationTime=60; % in seconds
    end
    
    emailToInformOfProgress='tv2@st-andrews.ac.uk';
    
    calibrationVoltages=[-5:5];
    calibrationWavelengths=(785+[.067 .096 .139 .186 .239 .296 .357 .414 .473 .535 .596])*1e-9;
    voltageWavelengthModel=fit(calibrationVoltages.',calibrationWavelengths.','smoothingspline');
    sampledWavelengths=voltageWavelengthModel(voltages);
    
    cam=PikeCam();
    cam.regionOfInterest=[0 0 480 320];
    cam.integrationTime=4000e-3; % initial integration time only
    cam.gain=630;
    
    startTime=clock();
    for idx=1:4,
        img=cam.acquire();
    end
    acquisitionTime=etime(clock(),startTime)/4;
    
    logMessage('Measuring %d images at %d wavelengths. This experiment will take approximately %.1f hours',[nbImagesPerWavelength length(sampledWavelengths),length(sampledWavelengths)*(laserStabilizationTime+nbImagesPerWavelength*acquisitionTime)/60/60]);
    
    DAQSession = daq.createSession('ni');
    DAQSession.addAnalogOutputChannel('Dev2',0,'Voltage');
    
    mkdir(outputFolder);
    
    logMessage('Setting the initial wavelength...');
    % Move to beginning of range in one second
    for (fraction=[0:.01:1])
        DAQSession.outputSingleScan(voltages(1)*fraction);
        pause(.01);
    end
    
    %Wait for user to start the experiment
    minutesToWait=input('Enter minutes to wait for before starting the experiment (default 0):');
    try
        if (isempty(minutesToWait))
            minutesToWait=0;
        else
            minutesToWait=minutesToWait(1);
        end
        if (minutesToWait>0)
            logMessage('Waiting for %0.1f minutes.',minutesToWait);
            pause(60*minutesToWait);
        end
        logMessage('Starting now...');

        % Scan all wavelengths
        images=[];
        for wavelengthIdx=1:length(sampledWavelengths),
            wavelength=sampledWavelengths(wavelengthIdx);
            logMessage('Setting the laser wavelength to %0.6f nm and waiting %0.1f seconds...',[wavelength*1e9 laserStabilizationTime]);
            DAQSession.outputSingleScan(voltages(wavelengthIdx));
            pause(laserStabilizationTime);
            logMessage('Recording...');
            for (imageIdx=1:nbImagesPerWavelength)
                img=cam.acquire();
                if (isempty(images))
                    images=zeros([size(img,1) size(img,2) nbImagesPerWavelength],'single');
                end
                images(:,:,imageIdx)=single(img);
            end
            logMessage('Peak intensity %f, interframe error %f, rel. interframe error %f',[max(images(:)) mean(mean(std(images,0,3))) mean(mean(std(images,0,3)))/mean(images(:))]);
            outputFileName=strcat(outputFolder,sprintf('/imagesForWavelength%0.6fnm.mat',wavelength*1e9));
            recordingTime=datestr(clock(),'dd-mmm-yyyy HH:MM:SS.FFF');
            camGain=cam.gain;
            camIntegrationTime=cam.integrationTime;
            camRegionOfInterest=cam.regionOfInterest;
            save(outputFileName,'sampledWavelengths','images','recordingTime','camGain','camIntegrationTime','camRegionOfInterest','voltages');  
            logMessage('---------------------- %0.1f%% done ------------------------',100*wavelengthIdx/length(sampledWavelengths));
            if (mod(wavelengthIdx,10)==1)
                progressMessage=sprintf('%0.1f%% done ---------------------------------',100*wavelengthIdx/length(sampledWavelengths));
                logMessage([progressMessage,'\n    and written to file: %s'],outputFileName,'tv2@st-andrews.ac.uk');
            end
        end

        logMessage('Returning to central wavelength...');
        % Move back to center in one second
        for (fraction=[1:-.01:0])
            DAQSession.outputSingleScan(voltages(end)*fraction);
            pause(.01);
        end

        delete(DAQSession);

        logMessage(['Sacher measurement - all done and written to folder ' strrep(outputFolder,'\','/')],[],emailToInformOfProgress);
    catch Exc
        logMessage('Sacher measurement - error: %s',Exc.message,emailToInformOfProgress);
        rethrow(Exc);
    end
end
