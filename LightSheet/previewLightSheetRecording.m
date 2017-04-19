% outputFullFileName=previewLightSheetRecording(fileNameGreen,fileNameRed)
%
% A function to show recorded light sheet movies before deconvolution.
% The second argument is optional.
% Instead of two arguments, both file names can be specified as a cell
% array of strings
%
% Returns the full file name path to a video recording.
%
function outputFullFileName=previewLightSheetRecording(fileNameGreen,fileNameRed)
    close all force;
    
    if nargin>1 && ~isempty(fileNameRed)
        fullFileNames={fileNameGreen,fileNameRed};
    else
        fullFileNames=fileNameGreen;
        if ~iscell(fullFileNames)
            fullFileNames={fullFileNames};
        end
    end
    clear fileNameRed fileNameGreen;
    
	if isempty(fullFileNames)
		[fileNames,pathName]=uigetfile({'*.avi;*.tif;*.tiff;recording_*.mat','Image Sequences (*.avi,*.tif,*.tiff,recording_*.mat)';...
            '*.avi','AVI Videos (*.avi)';...
            '*.tif;*.tiff','TIFF Image Sequences (*.tif,*.tiff)';...
            'recording_*.mat','Matlab Image Sequences (recording_*.mat)';...
            '*.*','All Files (*.*)'},...
            'Select Recording(s)','recording_lambda488nm_alpha0_beta100.avi','MultiSelect','on');
        if (~iscell(fileNames) && ~ischar(fileNames) && ~isempty(fileNames))
            logMessage('File open canceled.');
            return;
        end
        if (~iscell(fileNames))
            fileNames={fileNames};
        end
        fullFileNames={fullfile(pathName,fileNames{1})};
        if (length(fileNames)>1)
            fullFileNames(2)=fullfile(pathName,fileNames{2});
        end
    end
    
    try
        wavelength=regexpi(fullFileNames{1},'lambda([\d]*)nm_alpha[^\\/]+','tokens');
        wavelengths(1)=str2double(wavelength{1})*1e-9;
    catch Exc
        jsonFile=[fullFileNames{1}(1:end-4),'.json'];
        if (exist(jsonFile))
            setupConfig=loadjson(jsonFile);
            wavelengths(1)=setupConfig.excitation.wavelength;
        end
    end
    % If multiple files are specified, make sure that the first one is the green fluorescent one (488 excitation)
	if numel(fullFileNames)>=2
        wavelength=regexpi(fullFileNames{2},'lambda([\d]*)nm_alpha[^\\/]+','tokens');
        wavelengths(2)=str2double(wavelength{1})*1e-9;
        if (round(wavelengths(2)*1e9)==488) % Wrong order
            fullFileNames=fullFileNames{[2 1]};
            wavelengths=wavelengths([2 1]);
        end
    end

    nbChannels=length(fullFileNames);
    
    %Open all video files
    nbFrames=Inf;
	for (channelIdx=1:nbChannels)
        fullFileName=fullFileNames{channelIdx};
        logMessage('Reading channel from %s',fullFileName);
        fileExtension=lower(fullFileName(max(1,end-3):end));
        switch (fileExtension)
            case '.avi'
                dataSource=VideoReader(fullFileName);
                imgSize = [dataSource.Height, dataSource.Width];
                newNbFrames=dataSource.NumberOfFrames;
                nbFrames=min(nbFrames,newNbFrames);
            case '.mat'
                matFile=matfile(fullFileName,'Writable',false);
                dataSource=matFile;
                if (isprop(matFile,'restoredDataCube'))
                    imgSize=whos(matFile,'restoredDataCube');
                else
                    if (isprop(matFile,'recordedImageStack'))
                        imgSize=whos(matFile,'recordedImageStack');
                    else
                        logMessage('Did not find deconvolvedDataCube nor recordedImageStack in %s!',fullFileName);
                        return;
                    end
                end
                imgSize=imgSize.size; nbFrames=imgSize(3); imgSize=imgSize(1:2);
            otherwise
                % Must be multi-page image
                info=imfinfo(fullFileName);
                imgSize=[info(1).Height info(1).Width];
                nbFrames=length(info);
                dataSource=fullFileName;
        end
        dataSources{channelIdx}=dataSource;
    end
    
    [xRange yRange zRange alpha beta]=loadSettings(fullFileNames,[imgSize nbFrames]);
    
    nbChannels=length(fullFileNames);

	fig=figure('Visible','off','Units','pixels','CloseRequestFcn',@(obj,event) shutDown(obj),'Renderer','zbuffer');
    getBaseFigureHandle(fig);
	loadStatus();
    
    % Prepare video output.
    [outputPath, outputFileName]=fileparts(fullFileNames{1});
    if (length(fullFileNames)>1)
        [~,outputFileName2]=fileparts(fullFileNames{2});
        outputFileName=strcat(outputFileName,'_and_',outputFileName2);
    end
    outputFullFileName=fullfile(outputPath,strcat(outputFileName,'.mp4'));
    
    videoWriter=VideoWriter(outputFullFileName,'MPEG-4');
    open(videoWriter);
    
    % Open the figure window
    setUserData('drawing',true);
    windowMargins=[1 1]*128;
    initialPosition=get(0,'ScreenSize')+[windowMargins -2*windowMargins];
    set(fig,'Position',getStatus('windowPosition',initialPosition));
    set(fig,'ResizeFcn',@(obj,event) updateStatus('windowPosition',get(obj,'Position')));
    
    axs(1)=subplot(2,1,1,'Visible','off');
    axs(2)=subplot(2,1,2,'Visible','off');
    
    colorMaps{1}=@(values) interpolatedColorMap(values,[0 0 0; 0 1 0; 1 1 1],[0 .5 1]);
    colorMaps{2}=@(values) interpolatedColorMap(values,[0 0 0; 1 0 0; 1 1 1],[0 .5 1]);
    
    logMessage('Estimating SNR for %d channel(s)...',nbChannels);
    for (channelIdx=1:nbChannels)
        img=readDataCubeFromFile(dataSources{channelIdx},[],ceil(ceil(nbFrames/10)):ceil(nbFrames/5):nbFrames,false);
        normalizations(channelIdx)=1.333/max(img(:));
        backgrounds(channelIdx)=median(median(min(img,[],3)));
        logMessage('Channel %d: background=%0.2f%%, normalization=%0.1fx.',[channelIdx 100*backgrounds(channelIdx) normalizations(channelIdx)]);
    end
	
    logMessage('Displaying %d channel(s)...',nbChannels);
    set(fig,'Visible','on');
    projection=zeros(imgSize(2),nbFrames,nbChannels);
    setUserData('closing',false);
    frameIdx=0;
    while (~getUserData('closing') && frameIdx<nbFrames)
        frameIdx=frameIdx+1;
        img=zeros([imgSize nbChannels]);
        for (channelIdx=1:nbChannels)
            img(:,:,channelIdx)=readDataCubeFromFile(dataSources{channelIdx},[],frameIdx,false);
        end
        projectionSelection=1+floor(size(img,1)/2)+[-2:2];
        projection(:,nbFrames-(frameIdx-1),:)=max(img(projectionSelection,:,:),[],1);
                
        if (mod(frameIdx,2)==1 || frameIdx==nbFrames)
            % Create a color image
            if (nbChannels==1)
                wideField=(img-backgrounds(1))*normalizations(1);
                topView=(projection-backgrounds(1))*normalizations(1);
            else
                wideField=zeros([size(img,1) size(img,2) 3]);
                for channelIdx=1:nbChannels
                    wideField=wideField+mapColor((img(:,:,channelIdx)-backgrounds(channelIdx))*normalizations(channelIdx),colorMaps{channelIdx});
                end
                topView=zeros([size(projection,1) size(projection,2) 3]);
                for channelIdx=1:nbChannels
                    topView=topView+mapColor((projection(:,:,channelIdx)-backgrounds(channelIdx))*normalizations(channelIdx),colorMaps{channelIdx});
                end
            end
            
            margins=[75 75; 10 10; 0 0]; % before; center; after
            wideFieldSize=[diff(yRange([1 end])) diff(xRange([1 end]))];
            topViewSize=[diff(yRange([1 end])) diff(zRange([1 end]))];
            windowSize=get(fig,'Position'); windowSize=windowSize(3:4)-sum(margins);
            wideFieldSize=wideFieldSize.*windowSize(1)/wideFieldSize(1);
            topViewSize=topViewSize.*windowSize(1)/topViewSize(1);
            scaling=min(1,windowSize(2)/(wideFieldSize(2)+topViewSize(2)));
            wideFieldSize=wideFieldSize*scaling;
            topViewSize=topViewSize*scaling;
            posOffset=margins(1,:)+(windowSize-[wideFieldSize(1) wideFieldSize(2)+topViewSize(2)])/2;
            
            showImage(wideField,[],yRange*1e6,xRange*1e6,axs(1));
            set(axs(1),'Units','pixels','Position',[posOffset+[0 margins(2,2)+topViewSize(2)] wideFieldSize]);
            set(axs(1),'LineWidth',2,'TickDir','out','TickLength',[.01 .01]);
            set(axs(1),'XTick',[-1000:20:1000],'XTickLabel',[]);
            set(axs(1),'YTick',[-1000:20:1000],'FontSize',18,'FontWeight','bold');
            ylabel(axs(1),'swipe [\mum]','FontSize',20,'FontWeight','bold');
            showImage(permute(topView(:,end:-1:1,:),[2 1 3]),[],yRange*1e6,zRange*1e6,axs(2));
            set(axs(2),'Units','pixels','Position',[posOffset topViewSize]);
            set(axs(2),'LineWidth',2,'TickDir','out','TickLength',[.01 .01]);
            set(axs(2),'XTick',[-1000:20:1000],'FontSize',18,'FontWeight','bold');
            set(axs(2),'YTick',[-1000:20:1000],'FontSize',18,'FontWeight','bold');
            xlabel(axs(2),'propagation [\mum]','FontSize',20,'FontWeight','bold'); ylabel(axs(2),'scan [\mum]','FontSize',20,'FontWeight','bold');
            
            set(axs,'Visible','on');
            linkaxes(axs,'x');
            drawnow();
            
            writeVideo(videoWriter,getFrame(fig));
        end
        %Update the figure window title
        if (frameIdx<nbFrames)
            progressText=sprintf('%0.0f%%%% read - ',100*frameIdx/nbFrames);
        else
            progressText='';
        end
        wavelengthText=[', lambda=',sprintf('%0.0fnm, ',wavelengths*1e9)];
        if (length(wavelengths)<2), wavelengthText=wavelengthText(1:end-2); end
        set(gcf,'NumberTitle','off','Name',sprintf([progressText,'alpha=%0.0f, beta=%0.0f%%',wavelengthText],[alpha beta*100]));
    end
	
    % Close all video files
	for (channelIdx=1:nbChannels)
        dataSource=dataSources{channelIdx};
        if (isa(dataSource,'VideoReader'))
            delete(dataSource);
        end
    end
    
    close(videoWriter);

    setUserData('drawing',false);
    
    if (getUserData('closing'))
        closereq();
    end
end

function shutDown(obj,event)
    setUserData('closing',true);
    if (~getUserData('drawing'))
        closereq();
        % else that routine will take care of closing the window
    end
end


function [xRange yRange zRange alpha beta]=loadSettings(specificConfigFileNames,dataCubeSize)
    %Load the configuration
    functionName=mfilename();
    configPath=mfilename('fullpath');
    configPath=configPath(1:end-length(functionName));
    defaultConfigFileName=strcat(configPath,'/waterImmersion.json');
    defaultConfig=loadjson(defaultConfigFileName);
    if (~iscell(specificConfigFileNames))
        specificConfigFileNames={specificConfigFileNames};
    end
    % Find a json file
    fileIdx=1;
    while(fileIdx<=length(specificConfigFileNames) && ~exist([specificConfigFileNames{fileIdx}(1:end-4),'.json'],'file')),
        fileIdx=fileIdx+1;
    end
    if (fileIdx<=length(specificConfigFileNames))
        specificConfig=loadjson([specificConfigFileNames{fileIdx}(1:end-4),'.json']);
        config=structUnion(defaultConfig,specificConfig);
    else
        logMessage('Description file with extension .json is missing, trying to load .mat file.');
        try
            matFile=matfile(specificConfigFileNames{1},'Writable',false);
            config=matFile.setupConfig;
        catch Exc
            logMessage('Failed... assuming defaults!');
            config=defaultConfig;
        end
    end
    
    if (~isfield(config.detector,'center') || isempty(config.detector.center))
        config.detector.center=[0 0];
    end
    
    % Prepare the results
    stageTranslationStepSize=norm(median(diff(config.stagePositions.target)));    
    realMagnification=config.detection.objective.magnification*config.detection.tubeLength/config.detection.objective.tubeLength;
    xRange=-config.detector.center(1)+config.detector.pixelSize(1)*([1:dataCubeSize(1)]-floor(dataCubeSize(1)/2)-1)/realMagnification; % up/down
    yRange=-config.detector.center(2)+config.detector.pixelSize(2)*([1:dataCubeSize(2)]-floor(dataCubeSize(2)/2)-1)/realMagnification; % left/right
    zRange=stageTranslationStepSize*([1:dataCubeSize(3)]-floor(dataCubeSize(3)/2+1))/config.detector.framesPerSecond; %Translation range (along z-axis)
    
    alpha=config.modulation.alpha;
    beta=config.modulation.beta;
end

%
% Data access functions
%
%
% User data function are 'global' variables linked to this GUI, though not
% persistent between program shutdowns.
%
function value=getUserData(fieldName)
    fig=getBaseFigureHandle();
    userData=get(fig,'UserData');
    if (isfield(userData,fieldName))
        value=userData.(fieldName);
    else
        value=[];
    end
end
function setUserData(fieldName,value)
    if (nargin<2)
        value=fieldName;
        fieldName=inputname(1);
    end
    fig=getBaseFigureHandle();
    userData=get(fig,'UserData');
    userData.(fieldName)=value;
    set(fig,'UserData',userData);
end
function fig=getBaseFigureHandle(fig)
    persistent persistentFig;
    if (nargin>=1)
        if (~isempty(fig))
            persistentFig=fig;
        else
            clear persistentFig;
        end
    else
        fig=persistentFig;
    end
end
%
% 'Status' variables are stored persistently to disk and reloaded next time
%
function loadStatus(pathToStatusFile)
    if (nargin<1 || isempty(pathToStatusFile))
        pathToStatusFile=[mfilename('fullpath'),'.mat'];
    end
    try
        load(pathToStatusFile,'status');
    catch Exc
        status={};
        status.version=-1;
    end
    
    status.pathToStatusFile=pathToStatusFile;
       
    setUserData('status',status);
end
function value=getStatus(fieldName,defaultValue)
    status=getUserData('status');
    if (isfield(status,fieldName))
        value=status.(fieldName);
    else
        if (nargin<2)
            defaultValue=[];
        end
        value=defaultValue;
    end
end
function updateStatus(fieldName,value)
    status=getUserData('status');
    status.(fieldName)=value;
    setUserData('status',status);
    
    saveStatus();
end
function saveStatus()
    status=getUserData('status');
    
    save(status.pathToStatusFile,'status');
end

