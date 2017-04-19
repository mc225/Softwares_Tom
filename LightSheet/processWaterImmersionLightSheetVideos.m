% processWaterImmersionLightSheetVideos(folderNames,reprocessData,centerOffset,scanShear,perspectiveScaling)
%
% Reads the recorded date, converts it the matrices and deconvolves.
%
% Input:
%      folderNames: a string or cell array of strings indicating folders with avi and json files.
%      reprocessData: When false, date which already has a mat file with
%                     a restoredDataCube matrix in it will be skipped.
%      centerOffset: The offset of the light sheet beam waist in meters given
%                    as a two-element vector containing the vertical (swipe direction, Y) and the
%                    horizontal (propagation direction, X) offset, respectivelly.
%      scanShear: The fraction change in x and y (vertical and horizontal) when moving in z (axially)
%      perspectiveScaling: The linear change in [x,y]=[vertical,horizontal] magnification when increasing z by one meter, units m^-1.
%
%
% Note: projections can be shown using:
%         figure;axs=subplot(1,2,1);imagesc(yRange*1e6,zRange*1e6,squeeze(max(recordedImageStack,[],1)).');axis equal;axs(2)=subplot(1,2,2);imagesc(yRange*1e6,tRange*1e6,squeeze(max(restoredDataCube,[],1)).'); linkaxes(axs);
%         figure;imagesc(yRange*1e6,xRange*1e6,squeeze(max(recordedImageStack,[],3))); axis equal;
% 
function processWaterImmersionLightSheetVideos(folderNames,reprocessData,centerOffset,scanShear,perspectiveScaling)
    if (nargin<1 || isempty(folderNames))
%         folderNames={pwd()};
%         folderNames={'20120815_215406_gain300_amph1'};
%         folderNames='E:\RESULTS\TEST\2012-12-05 19_45_35.753';
%         folderNames='E:\RESULTS\TEST\2012-12-06 19_07_05.401';
%         folderNames='E:\RESULTS\TEST\2012-12-07 17_01_37.652';
%         folderNames='E:\RESULTS\TEST\600nmBeadsInAgarose';
%         folderNames='E:\RESULTS\TEST\600nmBeadsInAgarose_2';
%         folderNames='E:\RESULTS\TEST\600nmBeadsInAgarose_3';
%         folderNames='E:\RESULTS\TEST\600nmBeadsInAgarose_Andor';
%         folderNames={'E:\RESULTS\TEST\HEK_Andor','E:\RESULTS\TEST\HEK_Andor2'};
%         folderNames={'E:\RESULTS\TEST\HEK_Andor2','E:\RESULTS\TEST\2012-12-13 16_01_45.688'};
%         folderNames={'E:\RESULTS\TEST\2013-01-14 17_19_36.146_greenblue_2um'};
%         folderNames={'E:\RESULTS\2013_01_15_SmallBeadSample/','E:\RESULTS\2013_01_17_SmallBeadSample/','E:\RESULTS\2013_01_17_SmallBeadSample_2/'};
%         folderNames={'E:\RESULTS\Amphioxus1'};
%         folderNames={'E:\RESULTS\Amphioxus2'};
%         folderNames={'E:\RESULTS\TEST\2013-01-29 11_11_55.143_mixedSmallBeads_inAgarose_on_PDMS_green_blue'};
%         folderNames={'E:\RESULTS\smallBeadSample2'};
%         folderNames={'E:\RESULTS\TEST\DoubleRingSampleHolder\2013-02-12_green_blue'};
%         folderNames={'E:\RESULTS\NanoschellsInMCF7Spheroids'};
%         folderNames={'E:\RESULTS\NonClippedAiry'};
%         folderNames={'E:\RESULTS\NonClippedAiry3'};
%         folderNames={'E:\RESULTS\BaslerTestForPerspectiveError'};
%         folderNames={'E:\RESULTS\BaslerTestForPerspectiveError2\fullAp'};
%         folderNames={'E:\RESULTS\BaslerTestForPerspectiveError2\2013-02-20 11_49_09.543_halfAp_blue_green'};
%         folderNames={'E:\RESULTS\TEST\AciniMarch\NotFixedToSurface'};
%         folderNames={'E:\RESULTS\TEST\AciniMarch7'};
%         folderNames={'Z:\RESULTS\'};
%         folderNames={'Z:\RESULTS\03_27_2013\'};
%         folderNames={'Z:\RESULTS\2013-04-03\'};
%         folderNames={'Z:\RESULTS\03_27_2013\Acini_1\2013-03-27 10_38_45.594'};
%         folderNames={'Z:\RESULTS\03_27_2013\Bead_Test','Z:\RESULTS\03_27_2013\Acini_2\2013-03-27 10_43_10.564'};
%         folderNames={'F:\RESULTS\04_03_2013_bluebeads'};
%         folderNames={'F:\RESULTS\2013-04-04'};
%         folderNames={'F:\RESULTS\March22'};
%         folderNames={'F:\RESULTS\2013-04-05notmovingGaussian'};
%         folderNames={'F:\RESULTS\2013-04-05b','F:\RESULTS\2013-04-05d'};
%         folderNames={'F:\RESULTS\2013-04-08_multiColorBeads'};
%         folderNames={'F:\RESULTS\600nmBeads100micronScanIncBessel_BackApOneThird'};
%         folderNames={'F:\RESULTS\BeadsUnderAgarose100micronScan_BackApOneThird'};
%          folderNames={'F:\RESULTS\mirrorScan'};
%          folderNames={'F:\RESULTS\redbeads2013-04-18\2013-04-18 17_11_51.448'};
%          folderNames={'F:\RESULTS\redbeads2013-04-19halfAp\2013-04-19 10_42_17.196'};
%          folderNames={'F:\RESULTS\redbeads2013-04-19halfAp\2013-04-19 10_42_17.196'};
%          folderNames={'F:\RESULTS\2013-04-22beads\'};
%          folderNames={'F:\RESULTS\2013-04-22beads\2013-04-22 18_38_59.782\','F:\RESULTS\2013-04-22cellsHEK\','F:\RESULTS\2013-04-22cellsHEK2\'};
%         folderNames={'F:\RESULTS\2013-04-23beadsApOneThird'};
%         folderNames={'F:\RESULTS\2013-04-23HEK_ApOneThird'};
%         folderNames={'F:\RESULTS\2013-04-24HEK_ApOneThird','F:\RESULTS\2013-04-24Amphioxus'};
%          folderNames={'F:\RESULTS\redbeads2013-04-19halfAp\2013-04-19 10_42_17.196'};
%          folderNames={'Z:\RESULTS\2013-04-24HEK_ApOneThird\'};
%         folderNames={'Z:\RESULTS\2013-04-24HEK_ApOneThird\test'};
%         folderNames={'Z:\RESULTS\2013-04-26HEK_ApOneThird'};
%         folderNames={'Z:\RESULTS\2013-04-30 16_38_27.946_tiny'};
%         folderNames={'Z:\RESULTS\2013-04-26HEK_ApOneThird'};
%         folderNames={'Z:\RESULTS\2013-05-14_cellspheroidsWGA'};
%         folderNames={'Z:\RESULTS\2013-05-15'};
%         folderNames={'F:\RESULTS\2013-08-01_HEKdense','F:\RESULTS\2013-08-02_HEKdense1','F:\RESULTS\2013-08-02_HEKdense2','F:\RESULTS\2013-08-02_HEKdense3'};
%           folderNames={'F:\RESULTS\2013-08-08','F:\RESULTS\2013-08-09'};  % center 1051 603
%         folderNames={'F:\RESULTS\2013-08-29'};
%         folderNames={'H:\RESULTS\'};
%         folderNames={'F:\RESULTS\2013-08-29\HEKCellsApOneHalf'};
%         folderNames={'H:\RESULTS\2013-09-06'};
%         folderNames={'H:\RESULTS\2013-09-09_600nmRed_ApOneThird'};
%         folderNames={'H:\RESULTS\2013-09-11 13_17_32.106 Bessel in Fluorescene'};
%         folderNames={'H:\RESULTS\2013-09-11_600nmRed_ApOneThird_limitedScans\offsetScan'}; % [0 -63e-6]
%         folderNames={'H:\RESULTS\200nmBeadsApOneThird','H:\RESULTS\600nmBeadsApOneThird'}; % [0 0]
%         folderNames={'H:\RESULTS\2013-09-13_Amphioxus'}; % [0 0]
%         folderNames={'Z:\RESULTS\2013-04-10_Beads200nm'};
%         folderNames={'Z:\RESULTS\2013-04-10_Beads200nm\2013-04-10 16_10_20.097_ApOneHalf'};
%         folderNames={'Z:\RESULTS\2013-10-10_200nmBeadsNearEntranceInSmallAgaroseDroplet'};
%         folderNames={'E:\Stored Files\RESULTS\Amphioxus'};
        folderNames={'H:\RESULTS\2013-11-26_zebrafish_1','H:\RESULTS\2013-11-26_zebrafish_2'...
            ,'H:\RESULTS\2013-11-26_zebrafish_3','H:\RESULTS\2013-11-26_zebrafish_4','H:\RESULTS\2013-11-26_zebrafish_5'};
    end
    if (nargin<2) 
%         reprocessData=false;
        reprocessData=true;
    end
    if (nargin<3)
%         centerOffset=[0 -63e-6]; % [] == [0 0] means: don't change. % vertical horizontal
        centerOffset=[0 0];
    end
    if (nargin<4)
        scanShear=[]; % [0  0.09]; % [fraction] vertical horizontal
    end
    if (nargin<5)
        perspectiveScaling=[]; % [0  500]; % [m^-1] vertical horizontal
    end
    
    if (ischar(folderNames))
        folderNames={folderNames};
    end
    
    %Load the default configuration
    functionName=mfilename();
    configPath=mfilename('fullpath');
    configPath=configPath(1:end-length(functionName));
    defaultConfigFileName=strcat(configPath,'/waterImmersion.json');
    defaultConfig=loadjson(defaultConfigFileName);
    
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
            reprocessThisFile=true;
            if (~reprocessData && exist(outputFileName,'file'))
                storedVariables = whos('-file',outputFileName);
                if (ismember('restoredDataCube', {storedVariables.name}))
                    reprocessThisFile=false;
                    logMessage('Already done %s, skipping it!',outputFileName);
                end
            end
            if (reprocessThisFile)
                if (exist(configFile,'file'))
                    try
                        specificConfig=loadjson(configFile);
        %                 specificConfig.stagePositions.target=specificConfig.stagePositions.target*1e-6;
        %                 specificConfig.stagePositions.target=specificConfig.stagePositions.actual*1e-6;
        %                 savejson([],specificConfig,configFile);
                        setupConfig=structUnion(experimentConfig,specificConfig);
                    catch Exc
                        logMessage('Could not read json config file %s, assuming defaults!',configFile);
                        setupConfig=experimentConfig;
                    end
                else
                    logMessage('Description file with extension .json is missing, assuming defaults!');
                end

                % Load recorded data
                logMessage('Loading %s...',inputFileName);
                try
                    recordedImageStack=readDataCubeFromAviFile(inputFileName);
                catch Exc
                    logMessage('Failed to load data stack from %s!',inputFileName);
                    recordedImageStack=[];
                end

                if (~isempty(recordedImageStack))
                    % Store partial results
                    logMessage('Saving recorded data to %s...',outputFileName);
                    save(outputFileName,'recordedImageStack','setupConfig', '-v7.3');

                    % Deconvolve
                    logMessage('Starting image reconstruction...');
                    [recordedImageStack lightSheetDeconvFilter lightSheetOtf ZOtf xRange,yRange,zRange tRange lightSheetPsf]=deconvolveRecordedImageStack(recordedImageStack,setupConfig);
                    restoredDataCube=recordedImageStack; clear recordedImageStack; % This operation does not take extra memory in Matlab

                    % Append the rest of the results
                    logMessage('Saving restored data cube to %s...',outputFileName);
                    save(outputFileName,'restoredDataCube','xRange','yRange','zRange','tRange','ZOtf','lightSheetPsf','lightSheetOtf','lightSheetDeconvFilter','-append');
                    clear restoredDataCube;
                end
            else
                logMessage('Skipping file %s, already done.',inputFileName);
            end
        end
        
        % Check for subfolders and handles these recursively
        directoryList=dir(folderName);
        for listIdx=1:length(directoryList),
            if directoryList(listIdx).isdir && directoryList(listIdx).name(1)~='.'
                expandedFolderName=strcat(folderName,'/',directoryList(listIdx).name);
                processWaterImmersionLightSheetVideos(expandedFolderName,reprocessData);
            end
        end
        
    end
end
