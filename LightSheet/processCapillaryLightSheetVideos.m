% processCapillaryLightSheetVideos(folderNames,reprocessData,centerOffset,scanShear,perspectiveScaling,gaussianIlluminationStd,beamAngle)
%
% Reads the recorded date, converts it the matrices and deconvolves.
%
% Input:
%      folderNames: a string or cell array of strings indicating folders with avi and json files.
%      reprocessData: When false, date which already has a mat file with
%                     a restoredDataCube matrix in it will be skipped.
%      centerOffset: The offset of the light sheet beam waist in meters given
%                    as a two-elent vector containing the vertical (swipe direction, Y) and the
%                    horizontal (propagation direction, X) offset, respectivelly.
%      scanShear: The linear shift of the recording in [x,y]=[vertical,horizontal] when increasing z by one meter
%      perspectiveScaling: The linear change in [x,y]=[vertical,horizontal] magnification when increasing z by one meter
%
%
% Note: projections can be shown using:
%         figure;axs=subplot(1,2,1);imagesc(yRange*1e6,zRange*1e6,squeeze(max(recordedImageStack,[],1)).');axis equal;axs(2)=subplot(1,2,2);imagesc(yRange*1e6,tRange*1e6,squeeze(max(restoredDataCube,[],1)).'); linkaxes(axs);
%         figure;imagesc(yRange*1e6,xRange*1e6,squeeze(max(recordedImageStack,[],3))); axis equal;
% 
function processCapillaryLightSheetVideos(folderNames,reprocessData,centerOffset,scanShear,perspectiveScaling,gaussianIlluminationStd,beamAngle)
    if (nargin<1 || isempty(folderNames))
        folderNames={'Z:\RESULTS\20120501_alpha7_int1000_g1_littlemorepower_betterfocus_beads_leftofprobe\backwardScans'};
%         folderNames={'Z:\RESULTS\20120501_alpha7_int1000_g1_morepower_beads_leftofprobe'};
    end
    if (nargin<2) 
        reprocessData=true;
    end
    if (nargin<3)
        % the first two dimensions of the data cube [swipe propagation], x-y here, y-x in paper
        centerOffset=[0 42e-6]; % Leave empty [] to keep default
    end
    if (nargin<4)
        % the first two dimensions of the data cube [swipe propagation], x-y here, y-x in paper
        scanShear=[0.025 0.023]; % Leave empty [] to keep default
    end
    if (nargin<5)
        % the first two dimensions of the data cube [swipe propagation], x-y here, y-x in paper
        perspectiveScaling=[]; % [373.4900  722.6800]; % Leave empty [] to keep default
    end
    if (nargin<6)
        gaussianIlluminationStd=2/3;
    end
    if (nargin<7)
        beamAngle=0.05;
    end
    
    if (ischar(folderNames))
        folderNames={folderNames};
    end
        
    %Load the default configuration
    functionName=mfilename();
    configPath=mfilename('fullpath');
    configPath=configPath(1:end-length(functionName));
    defaultConfigFileName=strcat(configPath,'capillarySetup.json');
    defaultConfig=loadjson(defaultConfigFileName);
    
    defaultConfig.sample.signalLevel=0.3;

    % Go through all the specified folders
    for (folderName=folderNames(:).')
        folderName=folderName{1};
        logMessage('Checking folder %s for recorded videos.',folderName);
        experimentConfigFileName=fullfile(folderName,'experimentConfig.json');
        if (exist(experimentConfigFileName,'file'))
            experimentConfig=loadjson(experimentConfigFileName);
        else
            experimentConfig=struct();
            experimentConfig.detector=struct();
            experimentConfig.detector.center=[0 0];
            experimentConfig.detector.scanShear=[0 0];
            experimentConfig.detector.perspectiveScaling=[0 0];
            experimentConfig.excitation.illuminationClippingFactors=zeros(2); %0.06*[1 1; 1 1];
        end
        if (~isempty(centerOffset))
            experimentConfig.detector.center=centerOffset;
        end
        if (~isempty(scanShear))
            experimentConfig.detector.scanShear=scanShear;
        end
        if (~isempty(perspectiveScaling))
            experimentConfig.detector.perspectiveScaling=perspectiveScaling;
        end
        if (~isempty(gaussianIlluminationStd))
            experimentConfig.excitation.gaussianIlluminationStd=gaussianIlluminationStd;
        end
        if (~isempty(beamAngle))
            experimentConfig.excitation.beamAngle=beamAngle;
        end
%         experimentConfig.excitation.excitationApertureFilling=.92; % not used?
        savejson([],experimentConfig,experimentConfigFileName);
        experimentConfig=structUnion(defaultConfig,experimentConfig);
        % and process all videos
        videoFileNameList=dir(strcat(folderName,'/*.avi'));
        for (fileName={videoFileNameList.name})
            fileName=fileName{1}(1:end-4);
            filePathAndName=strcat(folderName,'/',fileName);
            logMessage('Processing %s...',filePathAndName);
            
            % Load specific configuration description
            configFile=strcat(filePathAndName,'.json');
            inputFileName=strcat(filePathAndName,'.avi');
            outputFileName=strcat(filePathAndName,'.mat');
            reprocessFile=true;
            if (~reprocessData && exist(outputFileName,'file'))
                storedVariables = whos('-file',outputFileName);
                if (ismember('restoredDataCube', {storedVariables.name}))
                    reprocessFile=false;
                end
            end
            if (reprocessFile)
                if (~exist(configFile,'file'))
                    logMessage('Description file with extension .json is missing, deducing parameters from file name!');
                    specificConfig={}; specificConfig.modulation={};
                    if (~isempty(regexp(configFile,'Airy_[^\\/]*\.json$', 'once')))
                        alpha=regexpi(configFile,'_alpha(\d+)_','tokens');
                        alpha=alpha{end};
                        specificConfig.modulation.alpha=-str2double(alpha);
                    else
                        specificConfig.modulation.alpha=0;
                    end
                    beta=regexp(configFile,'Bessel(\d+)_[^\\/]*\.json$','tokens');
                    if (~isempty(beta))
                        beta=beta{1};
                        specificConfig.modulation.beta=str2double(beta)/100;
                    else
                        specificConfig.modulation.beta=1;
                    end
                    stagePositions=csvread([configFile(1:end-5),'.txt'])*1e-6;
                    stagePositions=stagePositions*experimentConfig.sample.refractiveIndex*2; %Adjust for the capillary movement in air
                    specificConfig.stagePositions.actual=stagePositions;
                    specificConfig.stagePositions.target=median(diff(stagePositions))*([1:length(stagePositions)]-1)+stagePositions(1);
                    savejson([],specificConfig,configFile);
                else
                    specificConfig=loadjson(configFile);
                end
                setupConfig=structUnion(experimentConfig,specificConfig);

                % Load recorded data
                logMessage('Loading %s...',inputFileName);
                recordedImageStack=readDataCubeFromAviFile(inputFileName);
                %If backward scan, flip axis
                if (median(diff(specificConfig.stagePositions.target))<0)
                    specificConfig.stagePositions.target=specificConfig.stagePositions.target(end:-1:1);
                    specificConfig.stagePositions.actual=specificConfig.stagePositions.actual(end:-1:1);
                    recordedImageStack=recordedImageStack(:,:,end:-1:1,:);
                end

                % Store partial results
                logMessage('Saving recorded data to %s...',outputFileName);
                save(outputFileName,'recordedImageStack','setupConfig', '-v7.3');
                
                % Deconvolve
                logMessage('Starting image reconstruction...');
                [recordedImageStack lightSheetDeconvFilter lightSheetOtf ZOtf xRange,yRange,zRange tRange lightSheetPsf]=deconvolveRecordedImageStack(recordedImageStack,setupConfig);
                restoredDataCube=recordedImageStack; clear recordedImageStack; % This operation does not take extra memory in Matlab
                
                % Append the rest of the results
                logMessage('Saving restored data cube to %s...',outputFileName);
                save(outputFileName,'restoredDataCube','xRange','yRange','zRange','tRange','ZOtf','lightSheetPsf','lightSheetOtf','lightSheetDeconvFilter', '-v7.3','-append');
                clear restoredDataCube;
            else
                logMessage('Skipping file %s, already done.',inputFileName);
            end
        end
        
        % Check for subfolders and handles these recursively
        directoryList=dir(folderName);
        for listIdx=1:length(directoryList),
            if directoryList(listIdx).isdir && directoryList(listIdx).name(1)~='.'
                expandedFolderName=strcat(folderName,'/',directoryList(listIdx).name);
                processCapillaryLightSheetVideos(expandedFolderName,reprocessData);
            end
        end
        
    end
end
