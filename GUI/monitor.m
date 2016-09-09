%
% Executes the aberration correction GUI.
%
function varargout = monitor(varargin)
% MONITOR M-file for monitor.fig
%      MONITOR, by itself, creates a new MONITOR or raises the existing
%      singleton*.
%
%      H = MONITOR returns the handle to a new MONITOR or the handle to
%      the existing singleton*.
%
%      MONITOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MONITOR.M with the given input arguments.
%
%      MONITOR('Property','Value',...) creates a new MONITOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before monitor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to monitor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help monitor

    % Last Modified by GUIDE v2.5 24-Jun-2011 19:57:48

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @monitor_OpeningFcn, ...
                       'gui_OutputFcn',  @monitor_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    % End initialization code - DO NOT EDIT
    end

% --- Executes just before monitor is made visible.
function monitor_OpeningFcn(fig, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to monitor (see VARARGIN)

    % Choose default command line output for monitor
    handles.output = fig;

    % Update handles structure
    guidata(fig, handles);

    % UIWAIT makes monitor wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
    
    % Attempt to load the GUI status from disk
    loadStatus();
    updateStatus('version',0.99);
    
    setUserData('isClosing',false);
    
    setUserData('probing',false);
    setUserData('slm',[]);
    setUserData('camRequestedRegionOfInterest',[-1 -1 -1 -1]);
    setUserData('currentImageOnSLM',[]);
    setUserData('initialTimeOfGUI',clock());
    setUserData('frameUpdateTimes',NaN*ones(1,1000));
    setUserData('fullScreenManager',DefaultScreenManager.instance);
    
    fontSize=16;
    
    set(fig,'Units','normalized','Position',getStatus('mainWindowPosition',[0 .1 .8 .8]),'Resize','on','ResizeFcn',@(obj,event) updateStatus('mainWindowPosition',get(obj,'Position')));
    set(fig,'BackingStore','off');
%     set(fig,'DoubleBuffer','off'); %Required?
    cameraPanel=uipanel('Parent',fig,'Units','normalized','Position',[0 0 .5 1],'Title','Camera');
    slmPanel=uipanel('Parent',fig,'Units','normalized','Position',[.5 0 .5 1],'Title','Spatial Light Modulator');
    
    histogramAxes=axes('Parent',cameraPanel);
    set(histogramAxes,'Position',[.15 .875 .8 .1],'Units','normalized');
    set(histogramAxes,'XTickMode','manual','YTickMode','manual','XTick',[0:.1:1],'XTickLabel',[0:10:100]);
    set(histogramAxes,'XLim',[0 1]);
    set(histogramAxes,'FontName','Arial','FontWeight','bold','FontSize',fontSize*.8);
    xlabel(histogramAxes,'intensity [%]','FontName','Arial','FontWeight','bold','FontSize',fontSize);
    ylabel(histogramAxes,'# [%]','FontName','Arial','FontWeight','bold','FontSize',fontSize);
    setUserData('histogramAxes',histogramAxes);
    cameraAxes=axes('Parent',cameraPanel);
    set(cameraAxes,'Position',[.15 .15 .8 .6],'Units','normalized');
    set(cameraAxes,'FontName','Arial','FontWeight','bold','FontSize',fontSize*.8);
    xlabel(cameraAxes,'x [pixels]','FontName','Arial','FontWeight','bold','FontSize',fontSize);
    ylabel(cameraAxes,'y [pixels]','FontName','Arial','FontWeight','bold','FontSize',fontSize);
    setUserData('cameraAxes',cameraAxes);
    pupilAxes=axes('Parent',slmPanel);
    set(pupilAxes,'Position',[.05 .2 .8 .6],'Units','normalized');
    colormap(pupilAxes,gray(256));
    setUserData('pupilAxes',pupilAxes);
    
    set(fig,'CloseRequestFcn',@closeApp);
    cameraControlPanel=uipanel('Parent',cameraPanel,'Title','Camera Control','Units','pixels','Position',[10 10 320 40]);
    changeCamPopupMenu=uicontrol('Parent',cameraControlPanel,'Position',[5 5 80 20],'Tag','changeCam','Style', 'popupmenu', 'Callback', @changeCam);
    detectedCameras=detectCameras();
    setUserData('detectedCameras',detectedCameras);
    %See if we can find the camera of last time again
    defaultCameraType=getStatus('defaultCameraType','DummyCam');
    defaultCameraIndex=getStatus('defaultCameraIndex',1);
    defaultCameraDropDownIndex=find(strcmpi({detectedCameras.type},defaultCameraType));
    switch(length(defaultCameraDropDownIndex))
        case 0
            defaultCameraDropDownIndex=1; % Dummy Cam 1
        case 1
            % Index may have changed, but we don't care
        otherwise
            closerMatchForDefaultCameraDropDownIndex=find(strcmpi({detectedCameras.type},defaultCameraType) & [detectedCameras.index]==defaultCameraIndex);
            if (~isempty(closerMatchForDefaultCameraDropDownIndex))
                defaultCameraDropDownIndex=closerMatchForDefaultCameraDropDownIndex(1);
            else
                defaultCameraDropDownIndex=defaultCameraDropDownIndex(1); % Just pick the first one, even if the index doesn't match
            end
    end
    set(changeCamPopupMenu,'String',{detectedCameras.description},'Value',defaultCameraDropDownIndex);
    uicontrol('Parent',cameraControlPanel,'Position',[90 5 40 20],'Style', 'text','String','int. time:');
    integrationTimeEdit=uicontrol('Parent',cameraControlPanel,'Position',[130 5 40 20],'Tag','integrationTimeEdit','Style', 'edit', 'Callback', @adjustIntegrationTime,'String',getStatus('integrationTime','20'));
    setUserData('integrationTimeEdit',integrationTimeEdit);
    uicontrol('Parent',cameraControlPanel,'Position',[160 5 20 20],'Style', 'text','String','ms');
    uicontrol('Parent',cameraControlPanel,'Position',[180 5 25 20],'Style', 'text','String','gain:');
    gainEdit=uicontrol('Parent',cameraControlPanel,'Position',[205 5 20 20],'Tag','gainEdit','Style', 'edit', 'Callback', @adjustGain, 'String',getStatus('gain','1'));
    setUserData('gainEdit',gainEdit);
    uicontrol('Parent',cameraControlPanel,'Position',[225 5 20 20],'Style', 'text','String','#:');
    numberOfFramesEdit=uicontrol('Parent',cameraControlPanel,'Position',[240 5 20 20],'Tag','numberOfFramesEdit','Style', 'edit', 'Callback', @adjustNumberOfFrames,'String',getStatus('numberOfFrames','1'));
    setUserData('numberOfFramesEdit',numberOfFramesEdit);
    uicontrol('Parent',cameraControlPanel,'Position',[260 5 30 20],'Tag','darkButton','Style', 'togglebutton','String','Dark','FontWeight','bold','Callback',@updateDarkImage);
    pausePlayButton=uicontrol('Parent',cameraControlPanel,'Position',[290 5 20 20],'Tag','pausePlayButton','Style', 'togglebutton','String','','FontWeight','bold','Callback',@pausePlay);
    setUserData('pausePlayButton',pausePlayButton);
    setUserData('recording',false);

    changeCam(changeCamPopupMenu,[]);
    
    recordingPanel=uipanel('Parent',cameraPanel,'Title','Recording Control','Units','pixels','Position',[330 10 215 40]);
    uicontrol('Parent',recordingPanel,'Position',[5 5 60 20],'Style','pushbutton','String','Browse...','Callback',@selectOutputFile);
    outputFileEdit=uicontrol('Parent',recordingPanel,'Position',[65 5 110 20],'Style', 'edit','String','');
    setUserData('outputFileEdit',outputFileEdit);
    uicontrol('Parent',recordingPanel,'Position',[175 5 35 20],'Style', 'togglebutton','String','REC!','Callback',@toggleOutputRecording);
    
    %
    % Spatial Light Modulator
    %
    slmControlPanel=uipanel('Parent',slmPanel,'Title','SLM Control','Units','pixels','Position',[10 10 530 90]);
    changeSLMPopUpMenu=uicontrol('Parent',slmControlPanel,'Position',[5 55 90 20],'Tag','changeSLMPopUpMenu','Style', 'popupmenu', 'Callback', @updateSLMDisplayNumberOrSLM);
    set(changeSLMPopUpMenu,'String',{'phase only SLM','dual head SLM'},'Value',getStatus('slmTypeIndex',1));
    setUserData('changeSLMPopUpMenu',changeSLMPopUpMenu);
    uicontrol('Parent',slmControlPanel,'Position',[95 55 20 20],'Style','text','String','on:');
    slmDisplayNumberPopUpMenu=uicontrol('Parent',slmControlPanel,'Position',[112 55 100 20],'Tag','slmDisplayNumberPopUpMenu','Style', 'popupmenu', 'Callback', @updateSLMDisplayNumberOrSLM);
    populateSLMDisplayNumberPopupMenu(slmDisplayNumberPopUpMenu);
    setUserData('slmDisplayNumberPopUpMenu',slmDisplayNumberPopUpMenu);
    setUserData('slmAxes',[]);
    defaultDisplayIndex=getStatus('slmDisplayNumber',1);
    if (defaultDisplayIndex>length(get(slmDisplayNumberPopUpMenu,'String')))
        defaultDisplayIndex=1; % popup
    end
    set(slmDisplayNumberPopUpMenu,'Value',defaultDisplayIndex);
    updateSLMDisplayNumberOrSLM(slmDisplayNumberPopUpMenu);
    setUserData(slmDisplayNumberPopUpMenu);
    uicontrol('Parent',slmControlPanel,'Position',[215 55 20 20],'Style','text','String','2pi=');
    twoPiEquivalentEdit=uicontrol('Parent',slmControlPanel,'Position',[235 55 25 20],'Tag','twoPiEquivalentEdit','Style', 'edit', 'Callback', @adjustTwoPiEquivalent,'String',getStatus('twoPiEquivalent','100'));
    uicontrol('Parent',slmControlPanel,'Position',[260 55 20 20],'Tag','estimateTwoPiEquivalentButton','Style', 'pushbutton','String','%!','FontWeight','bold','Callback',@estimateTwoPiEquivalentCallback);
    %uicontrol('Parent',slmControlPanel,'Position',[263 55 10 20],'Style','text','String','%');
    setUserData('twoPiEquivalentEdit',twoPiEquivalentEdit);
    uicontrol('Parent',slmControlPanel,'Position',[280 55 95 20],'Style','text','String','deflect(x,y,{w20}):');
    deflectionEdit=uicontrol('Parent',slmControlPanel,'Position',[370 55 70 20],'Tag','deflectionEdit','Style', 'edit', 'Callback', @adjustDeflection,'String',getStatus('deflection','1/10 1/10 0'));
    uicontrol('Parent',slmControlPanel,'Position',[445 55 30 20],'Style','text','String','pix^-1');
    setUserData('deflectionEdit',deflectionEdit);   
    uicontrol('Parent',slmControlPanel,'Position',[5 30 50 20],'Style','text','String','correction:');
    correctionEdit=uicontrol('Parent',slmControlPanel,'Position',[60 30 250 20],'Tag','correctionEdit','Style', 'edit', 'Callback',@adjustCorrection,'String',getStatus('correction',''));
    uicontrol('Parent',slmControlPanel,'Position',[310 30 60 20],'Style','pushbutton','String','Browse...','Callback',@loadCorrection);
    probeMethodToggle=uicontrol('Parent',slmControlPanel,'Position',[370 30 15 20],'Tag','probeMethod','Style','togglebutton','String',getStatus('probeMethodToggle','M'),'FontWeight','bold','Callback',@toggleProbeMethod);
    setUserData('probeMethodToggle',probeMethodToggle);
    uicontrol('Parent',slmControlPanel,'Position',[385 30 50 20],'Tag','estimateCorrectionButton','Style','pushbutton','String','Probe!','FontWeight','bold','Callback',@estimateCorrection);
    setUserData('correctionEdit',correctionEdit);
    uicontrol('Parent',slmControlPanel,'Position',[435 30 25 20],'Style','text','String','size:');
    probeSizeEdit=uicontrol('Parent',slmControlPanel,'Position',[460 30 20 20],'Tag','probeSizeEdit','Style', 'edit', 'Callback',@updateProbeSize,'String',getStatus('probeSize','1'));
    setUserData('probeSizeEdit',probeSizeEdit);
    updateProbeSize();
    uicontrol('Parent',slmControlPanel,'Position',[480 30 25 20],'Style','text','String','amp/');
    amplificationLimitEdit=uicontrol('Parent',slmControlPanel,'Position',[505 30 20 20],'Tag','amplificationLimitEdit','Style', 'edit', 'Callback',@updateAmplificationLimit,'String',getStatus('amplificationLimit','1'));
    setUserData('amplificationLimitEdit',amplificationLimitEdit);
    uicontrol('Parent',slmControlPanel,'Position',[5 5 30 20],'Style','text','String','pupil:');
    maskEdit=uicontrol('Parent',slmControlPanel,'Position',[45 5 280 20],'Tag','maskEdit','Style', 'edit', 'Callback', @updateSLMDisplayAndOutput,'String',getStatus('mask','R<200'));
    uicontrol('Parent',slmControlPanel,'Position',[330 5 30 20],'Style','text','String','pix^-1');
    uicontrol('Parent',slmControlPanel,'Position',[360 5 20 20],'Style','pushbutton','String','T','FontWeight','bold','Callback', @(obj,event) insertPupilRadiusAndUpdateSLM('R<pupilRadius'));
    uicontrol('Parent',slmControlPanel,'Position',[380 5 20 20],'Style','pushbutton','String','B','FontWeight','bold','Callback', @(obj,event) insertPupilRadiusAndUpdateSLM('R>0.9*pupilRadius & R<pupilRadius'));
    uicontrol('Parent',slmControlPanel,'Position',[400 5 20 20],'Style','pushbutton','String','V','FontWeight','bold','Callback', @(obj,event) insertPupilRadiusAndUpdateSLM('(Rho<pupilRadius)*exp(i*Phi)'));
    uicontrol('Parent',slmControlPanel,'Position',[420 5 20 20],'Style','pushbutton','String','A','FontWeight','bold','Callback', @(obj,event) insertPupilRadiusAndUpdateSLM('(R<pupilRadius)*exp(3*2i*pi*((X/pupilRadius)^3+(Y/pupilRadius)^3))'));
    uicontrol('Parent',slmControlPanel,'Position',[440 5 20 20],'Style','pushbutton','String','+','FontWeight','bold','Callback', @(obj,event) insertPupilRadiusAndUpdateSLM('(R<pupilRadius)*(abs(R*cos(P-t*pi/20))<pupilRadius/50 | abs(R*sin(P-t*pi/20))<pupilRadius/50)'));
    setUserData('maskEdit',maskEdit); 
    
%     align(HandleList,'Fixed',5,'Fixed',10);

    changeSLM(changeSLMPopUpMenu,[]);
    
    dragzoom(cameraAxes);
    set(fig,'Name','monitor - SLM aberration correction GUI');
end

% --- Outputs from this function are returned to the command line.
function varargout = monitor_OutputFcn(fig, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
%     varargout{1} = handles.output;
    pausePlay(); %Start recording
    
    varargout={};
end

function updatePupilModulationOnSLM()
    pupilFunction=getUserData('pupilFunction');
    lastPupilFunctionUpdate=getUserData('lastPupilFunctionUpdate');
    initialTimeOfSLM = getUserData('initialTimeOfSLM');
    timeNow=etime(clock(),initialTimeOfSLM);
    
    slm=getUserData('slm');
    regionOfInterestSize=slm.regionOfInterest; regionOfInterestSize=regionOfInterestSize(3:4);
    regionOfInterestCenter=floor(regionOfInterestSize./2)+1;
    
    if (lastPupilFunctionUpdate<0 || getUserData('pupilFunctionTimeDependent'))
        [X,Y]=meshgrid([1:regionOfInterestSize(2)]-regionOfInterestCenter(2),[1:regionOfInterestSize(1)]-regionOfInterestCenter(1));
        if (lastPupilFunctionUpdate<0 || any(any(pupilFunction(X,Y,lastPupilFunctionUpdate)~=pupilFunction(X,Y,timeNow))))
            refreshPupilModulationOnSLM();

            %
            % Display also on screen
            %
            % Remove image if the size is different
            pupilAxes=getUserData('pupilAxes');
            imgObject=get(pupilAxes,'Children');
            if (~isempty(imgObject))
                %Determine current image size
                currentImgSize=size(get(imgObject,'CData'));
                %Check size and remove if needed
                if (any(currentImgSize(1:2)~=regionOfInterestSize))
                    set(imgObject,'CData',[]);
                end
            end
            % Update pupil display
            pupil=pupilFunction(X,Y,timeNow);
            if (max(abs(imag(pupil(:))))<10*eps(class(pupil)))
                pupil=pupil.*exp(1e-6i); %Make complex for display
            end
            showImage(pupil,[],X,Y,pupilAxes);
        end
    end
end
function refreshPupilModulationOnSLM()
    pupilFunction=getUserData('pupilFunction');
    initialTimeOfSLM = getUserData('initialTimeOfSLM');
    timeNow=etime(clock(),initialTimeOfSLM);
    
    slm=getUserData('slm');
    
    % Set the stabilization time temporarily to zero because for this
    % situation we do not care if the camera is out-of-sync.
    origStabilizationTime=slm.stabilizationTime;
    slm.stabilizationTime=0;
    slm.modulate(@(X,Y) pupilFunction(X,Y,timeNow) );
    slm.stabilizationTime=origStabilizationTime;
    setUserData('lastPupilFunctionUpdate',timeNow);
end
function acquireAndDisplayLoop()
    while (~getUserData('isClosing') && getUserData('recording')),
        acquireAndDisplay();
%         pause(0.05);
    end
end
function acquireAndDisplay()
    %Update the SLM if the pupil function is time dependent
    if (~getUserData('probing'))
        updatePupilModulationOnSLM();
    end
    
    camAx=getUserData('cameraAxes');
    xLim=get(camAx,'XLim'); yLim=get(camAx,'YLim');
    newROI=round([yLim(1) xLim(1) diff(yLim) diff(xLim)]);

    %Acquire a new image
    cam=getUserData('cam');
    updateCamROI(newROI);
    img=cam.acquire();
    %Draw image
    displayOnMonitor(img);
end
%Draw image
function displayOnMonitor(img)
    cam=getUserData('cam');
    
    %Check if anybody just changed the camera
    if (any(size(img,1)>cam.regionOfInterest(3:4)))
        img=zeros(cam.regionOfInterest(3:4));
    end
    
    nonSaturatedPixels=img+cam.background<1;
    normalizedImg=img;
    if (all(nonSaturatedPixels(:)))
        maxLevel=max(img(:));
        normalizedImg=img/maxLevel; %Normalize image intensity to maximum
    else
        maxLevel=1;
    end
    normalizedImg=uint8(normalizedImg*256); %map to [0:255] interval
    %Mark saturated pixels as red
    colorImg=repmat(normalizedImg.*uint8(nonSaturatedPixels),[1 1 3]); % Only copy the non-saturated, leaving other pixels black
    colorImg(~nonSaturatedPixels)=255; %Mark as red when
    
    completeImage=zeros([cam.maxSize 3],'uint8');
    try
        completeImage(cam.regionOfInterest(1)+[1:cam.regionOfInterest(3)],cam.regionOfInterest(2)+[1:cam.regionOfInterest(4)],:)=colorImg;
    catch Exc
        Exc
        cam
        cam.regionOfInterest
        size(completeImage)
        size(colorImg)
        rethrow(Exc);
    end
    camAx=getUserData('cameraAxes');
    imgHandle=get(camAx,'Children');
    while(~ishandle(imgHandle(end)) || ~strcmp(get(imgHandle(end),'Type'),'image'))
        imgHandle=imgHandle(1:end-1); %Remove elements that dragzoom added to it
    end
    imgHandle=imgHandle(end);
    set(imgHandle,'EraseMode','none');
    set(imgHandle,'CData',completeImage);
    
    %Plot histogram
    histAx=getUserData('histogramAxes');
    nbGrayLevelsForHistogram=256;
    relPixelCountsPerGrayLevel=histc(img(:),[0:nbGrayLevelsForHistogram]/nbGrayLevelsForHistogram)./numel(img);
    relPixelCountsPerGrayLevel(end-1)=relPixelCountsPerGrayLevel(end-1)+relPixelCountsPerGrayLevel(end); %The last bin has the border cases of value == 1.0
    relPixelCountsPerGrayLevel=relPixelCountsPerGrayLevel(1:end-1).'; % Remove last bin after copying in penultimate
    cla(histAx);
    patch(repmat([0.5:nbGrayLevelsForHistogram]/nbGrayLevelsForHistogram,[4 1])+repmat([-1 -1 1 1].'*.5/nbGrayLevelsForHistogram,[1 nbGrayLevelsForHistogram]),[0 1 1 0].'*relPixelCountsPerGrayLevel,[0:nbGrayLevelsForHistogram-1],'LineStyle','none','Parent',histAx);
    hold(histAx,'on');
    plot((maxLevel+0.01)*[1; 1],[0; 1],'-','Color',[.8 .8 .8],'LineWidth',3,'Parent',histAx);
    yLims=[0 max(1e-4,10^ceil(log10(2*max(relPixelCountsPerGrayLevel(2:end-1)))))];
    set(histAx,'YLim',yLims);
    set(histAx,'YTick',[0 yLims(2)/2 yLims(2)],'YTickLabel',{'0',formatNumberForTickLabel(100*yLims(2)/2,100*yLims(2)/2),formatNumberForTickLabel(100*yLims(2),100*yLims(2)/2)});
    hold(histAx,'off');
    colorMap=gray(256); colorMap(end,2:3)=0; %Add red marker for saturation
    colormap(histAx,colorMap);
    set(histAx,'Color',[0 .333 1]);
    
    %Show the frame rate
    initialTimeOfGUI=getUserData('initialTimeOfGUI');
    frameUpdateTimes=getUserData('frameUpdateTimes');
    frameUpdateTimes=[frameUpdateTimes(2:end) etime(clock(),initialTimeOfGUI)];
    setUserData('frameUpdateTimes',frameUpdateTimes);
    relativeTimes=frameUpdateTimes-frameUpdateTimes(end);
    weights=exp(relativeTimes);
    averageFrameUpdateTime=[0 diff(relativeTimes(~isnan(relativeTimes)))]*weights(~isnan(relativeTimes)).'/sum(weights(~isnan(relativeTimes)));
    title(camAx,sprintf('%0.2f fps',1/averageFrameUpdateTime),'FontName','Arial','FontWeight','bold','FontSize',16);
    
    drawnow();
end
function shutDown()
    %Keep the current window position
    get(getBaseFigureHandle(),'Position');
    updateStatus('mainWindowPosition',get(getBaseFigureHandle(),'Position'));

    %Clean up before exit
    cam=getUserData('cam');
    cam.delete();
    
%     slmAxes=getUserData('slmAxes');
%     if (ishandle(slmAxes))
%         slmFigure=get(slmAxes,'Parent');
%         updateStatus('slmFigurePosition',get(slmFigure,'Position'));
%         close(slmFigure);
%     else
%         getUserData('fullScreenManager').delete();
%     end
    
    closereq();
end

%% Callbacks
function adjustIntegrationTime(obj,event)
    if (nargin<1)
        obj=getUserData('integrationTimeEdit');
    end
    newIntegrationTime=str2num(get(obj,'String'))*1e-3;
    if (~isempty(newIntegrationTime))
        newIntegrationTime=newIntegrationTime(1);
    end
    if (isnumeric(newIntegrationTime) && ~isnan(newIntegrationTime))
        newIntegrationTime=abs(newIntegrationTime);
        cam=getUserData('cam');
        cam.integrationTime=newIntegrationTime;
        setUserData('cam',cam);
        
        set(obj,'String',sprintf('%g',cam.integrationTime*1e3));
    end
    
    updateStatus('integrationTime',get(obj,'String'));
end
function adjustGain(obj,event)
    if (nargin<1)
        obj=getUserData('gainEdit');
    end
    newGain=str2num(get(obj,'String'));
    if (isnumeric(newGain) && ~isnan(newGain))
        newGain=abs(newGain(1));
        cam=getUserData('cam');
        cam.gain=newGain;
        setUserData('cam',cam);
        
        set(obj,'String',sprintf('%g',cam.gain));
    end
    
    updateStatus('gain',get(obj,'String'));
end
function adjustNumberOfFrames(obj,event)
    if (nargin<1)
        obj=getUserData('numberOfFramesEdit');
    end
    newNumberOfFrames=str2num(get(obj,'String'));
    if (isnumeric(newNumberOfFrames) && ~isnan(newNumberOfFrames))
        newIntegrationTime=max(1,round(newNumberOfFrames(1)));
        cam=getUserData('cam');
        cam.defaultNumberOfFramesToAverage=newNumberOfFrames;
        setUserData('cam',cam);
    end
    
    updateStatus('numberOfFrames',get(obj,'String'));
end
function updateDarkImage(obj,event)
    cam=getUserData('cam');
    
    if (get(obj,'Value')>0)
        logMessage('Updating dark image...');

        cam=cam.acquireBackground();
    else
        logMessage('Not using dark image.');
        cam.background=0;
    end
    
    setUserData('cam',cam);
end
function pausePlay(obj,event)
    if (nargin<1)
        obj=getUserData('pausePlayButton');
    end
    
    recording=~getUserData('recording');
        
    setUserData('recording',recording);
    pausePlayButton=getUserData('pausePlayButton');
    set(pausePlayButton,'Value',recording);
    setUserData('pausePlayButton',pausePlayButton);
    
    switch(recording)
        case true
            set(obj,'String','||');
        case false
            set(obj,'String','>>');
    end
    if (recording && ~getUserData('probing'))
        acquireAndDisplayLoop();
    end
    if (getUserData('isClosing'))
        shutDown();
    end
end
% function adjustGain(obj,event)
%     if (nargin<1)
%         obj=getUserData('gainEdit');
%     end
%     newGain=str2num(get(obj,'String'));
%     if (isnumeric(newGain) && ~isnan(newGain))
%         newGain=round(abs(newGain(1)));
%         cam=getUserData('cam');
%         cam.gain=newGain;
%         setUserData('cam',cam);
%     end
% end
function changeCam(obj,event)
    selectedCamIdx=get(obj,'Value');
    detectedCameras=getUserData('detectedCameras');
    if (selectedCamIdx<1 || selectedCamIdx>length(detectedCameras))
        selectedCamIdx=1;
    end
    selectedCamera=detectedCameras(selectedCamIdx);
    
    try
        switch(selectedCamera.type),
            case 'BaslerGigE'
                %Basler GigE
                cam=BaslerGigECam(selectedCamera.index);
            case 'Andor'
                %Andor
                logMessage('Andor camera not implemented');
            case 'ImagingSource'
                %Imaging Source
                cam=ImagingSourceCam(selectedCamera.index);
            case 'DirectShow'
                %Direct Show
                cam=DirectShowCam(selectedCamera.index);
            case 'PikeCam'
                %Pike Show
                cam=PikeCam(selectedCamera.index);
            otherwise
                %DummyCam
                switch (selectedCamera.index)
                    case 2
                        imgSize=[480 640];
                    otherwise
                        imgSize=[1 1]*256;
                end
                cam=DummyCam(imgSize);

                cam.acquireDirectFunctor=@(cam,nbFrames) calcDummyImage(imgSize,cam,nbFrames);
        end
        setUserData('cam',cam);
    
        updateStatus('cam',get(obj,'Value'));
        
        updateStatus('defaultCameraType',selectedCamera.type);
    	updateStatus('defaultCameraIndex',selectedCamera.index);

        cameraAxes=getUserData('cameraAxes');
%         dragzoom(cameraAxes,'off');
        initializeCamROIAndAxes();
        dragzoom(cameraAxes); % Update the default FOV
        adjustIntegrationTime();
    catch Exc
        logMessage(strcat('Could not initialize camera: ',selectedCamera.description,', falling back to Dummy Cam...'));
        set(obj,'Value',1);
        changeCam(obj,event);
    end
end

function updateCamROI(newROI)
    cam=getUserData('cam');
    
    %Clip the new region of interest
    newROI(3:4)=max(min(newROI(3:4),cam.maxSize),1);
    newROI(1:2)=min(max(0,newROI(1:2)),cam.maxSize-newROI(3:4));
    
    if (~all(getUserData('camRequestedRegionOfInterest')==newROI))
        cam.regionOfInterest=newROI;
        setUserData('cam',cam);
        setUserData('camRequestedRegionOfInterest',newROI);
    end
end
function initializeCamROIAndAxes()
    cam=getUserData('cam');
    roi=[0 0 cam.maxSize];
    updateCamROI(roi);
        
    roiSize = roi(3:4);
    ax = getUserData('cameraAxes');
    hImage = image(roiSize,'Parent',ax); axis(ax,'equal');
    set(ax,'XLim',[0 roiSize(2)],'YLim',[0 roiSize(1)]);
end

%% SLM callbacks
function adjustDeflection(obj,event)
    if (nargin<1)
        obj=getUserData('deflectionEdit');
    end
    newDeflection=str2num(get(obj,'String'));
    if (isnumeric(newDeflection) && ~any(isnan(newDeflection)))
        if (length(newDeflection)<2)
            logMessage('The deflection frequency should be two or three scalars, assuming diagonal deflection.');
            newDeflection(2)=newDeflection(1);
        end
        if (length(newDeflection)<3)
            newDeflection(3)=0;
        end
        if (length(newDeflection)>3)
            logMessage('The deflection frequency should be maximum 3 scalars.');
            newDeflection=newDeflection(1:3);
        end
        slm=getUserData('slm');
        slm.referenceDeflectionFrequency=newDeflection([2 1 3]);
        setUserData('slm',slm);
        
        updateSLMDisplayAndOutput();
        
        updateStatus('deflection',get(obj,'String'));
    end
end
function adjustTwoPiEquivalent(obj,event)
    if (nargin<1)
        obj=getUserData('twoPiEquivalentEdit');
    end
    newTwoPiEquivalent=str2num(get(obj,'String'));
    if (isnumeric(newTwoPiEquivalent) && ~any(isnan(newTwoPiEquivalent)))
        newTwoPiEquivalent=newTwoPiEquivalent(1)/100;
        slm=getUserData('slm');
        slm.twoPiEquivalent=newTwoPiEquivalent;
        setUserData('slm',slm);
        
        updateSLMDisplayAndOutput();
        
        updateStatus('twoPiEquivalent',get(obj,'String'));
    end
end
function adjustCorrection(obj,event)
    if (nargin<1)
        obj=getUserData('correctionEdit');
    end
    correctionFunctionFileName=strtrim(get(obj,'String'));
    
    slm=getUserData('slm');
    
    if (~isempty(correctionFunctionFileName))
        if (~strcmp(correctionFunctionFileName(end-3:end),'.mat'))
            correctionFunctionFileName=strcat(correctionFunctionFileName,'.mat');
        end
        try
            correctionFileContents=whos('-file',correctionFunctionFileName);
            correctionFileContents={correctionFileContents.name};
            % Configure the region of interest of the SLM as given in the correction archive
            if (any(strcmp(correctionFileContents,'slmRegionOfInterest')))
                load(correctionFunctionFileName,'slmRegionOfInterest');
                slm.regionOfInterest=slmRegionOfInterest;
            else
                logMessage('slmRegionOfInterest variable not found in the correction file %s, leaving the region-of-interest unchanged.',correctionFunctionFileName);
            end
            % If the amplification limit is specified in the archive,
            % recalculate the correction if the original info is known.
            amplificationLimit=getStatus('amplificationLimit');
            if (isempty(amplificationLimit) && any(strcmp(correctionFileContents,'measuredPupilFunction')))
                load(correctionFunctionFileName,'amplificationLimit');                
                % Update the GUI with the amplification limit read from file
                updateStatus('amplificationLimit',amplificationLimit);
                setUserData('amplificationLimitEdit',num2str(amplificationLimit));
            end
            if (~isempty(amplificationLimit) && any(strcmp(correctionFileContents,'measuredPupilFunction')))
                load(correctionFunctionFileName,'measuredPupilFunction');
                if (any(strcmp(correctionFileContents,'initialCorrection')))
                    load(correctionFunctionFileName,'initialCorrection');
                else
                    initialCorrection=1;
                end
                slm.correctionFunction=calcCorrectionFromPupilFunction(measuredPupilFunction./initialCorrection,amplificationLimit);
                logMessage('Correction function updated using amplification limit %d.',amplificationLimit);
            else
                logMessage('No amplification limit specified, using the correction as specified in the file %s.',correctionFunctionFileName);
                load(correctionFunctionFileName,'pupilFunctionCorrection');
                slm.correctionFunction=pupilFunctionCorrection;
            end
        catch Exc
            logMessage('Couldn''t load %s file, not a valid Matlab archive.',correctionFunctionFileName);
        end
    else
        slm.correctionFunction=1;
        slm.regionOfInterest=[];
        logMessage('No correction loaded.');
    end
    
    setUserData('slm',slm);

    updateStatus('correction',correctionFunctionFileName);
    
    refreshPupilModulationOnSLM();
end
function loadCorrection(obj,event)
    correctionEdit=getUserData('correctionEdit');
    
    [fileName,pathName,filterIndex] = uigetfile({'*.mat'},'Select Correction Function','pupilFunctionCorrection.mat'); 
    
    if (~isempty(fileName))
        set(correctionEdit,'String',strcat(pathName,fileName));
    end
    
    adjustCorrection();
end
function selectOutputFile(obj,event)
    outputFileEdit=getUserData('outputFileEdit');
    
    % Suggest a non-existing file name
    defaultFileName='recording.avi';
    if (exist(defaultFileName,'file'))
        idx=2;
        while (exist(strcat(defaultFileName,'_',num2str(idx)),'file'))
            idx=idx+1;
        end
        defaultFileName=strcat(defaultFileName,'_',idx);
    end
    [fileName,pathName,filterIndex] = uiputfile({'*.mat';'*.avi';'*.png';'*.tif';'*.jpg'},'Save as...',defaultFileName);
    
    if (~isempty(fileName))
        set(outputFileEdit,'String',strcat(pathName,fileName));
    end
end
function updateProbeSize(obj,event)
    if (nargin<1)
        obj=getUserData('probeSizeEdit');
    end
    probeSize=str2num(get(obj,'String')).';
    probeSize=probeSize(1:min(2,end));
    if isempty(probeSize) || any(isnan(probeSize))
        probeSize=1;
    end
    probeSize=max(1,round(probeSize));
    if (length(probeSize)<2)
        probeSize(2)=probeSize(1);
    end
    
    updateStatus('probeSize',probeSize);
    set(obj,'String',sprintf('%d ',probeSize));
    setUserData('probeSize',probeSize);
end
function updateAmplificationLimit(obj,event)
    if (nargin<1)
        obj=getUserData('amplificationLimitEdit');
    end
    amplificationLimit=str2num(get(obj,'String'));
    if (length(amplificationLimit)>1)
        amplificationLimit=amplificationLimit(1);
    end
    if (~isempty(amplificationLimit) && isnan(amplificationLimit))
        amplificationLimit=[];
    end
    amplificationLimit=max(1,amplificationLimit);
    
    updateStatus('amplificationLimit',amplificationLimit);
    setUserData('amplificationLimitEdit',num2str(amplificationLimit));
    
    % If a correction is already selected, update the correction pattern using the new amplification limit.
    adjustCorrection();
end
function updateSLMDisplayNumberOrSLM(obj,event)
    slm=getUserData('slm');
    if (~isempty(slm))
        delete(slm); % Close any current SLM windows
    end
    
    slmDisplayNumberPopUpMenu=getUserData('slmDisplayNumberPopUpMenu');
    slmDisplayNumber=get(slmDisplayNumberPopUpMenu,'Value');
    descriptors=get(slmDisplayNumberPopUpMenu,'String');
    if (slmDisplayNumber<1 || slmDisplayNumber>length(descriptors))
        slmDisplayNumber=1;
    end
    if (strcmpi(descriptors{slmDisplayNumber},'popup'))
        slmAxes=getUserData('slmAxes');
        if (isempty(slmAxes) || ~ishandle(slmAxes))
            mainFigHandle=getBaseFigureHandle();
            slmFigure=figure('Name','Dummy SLM','NumberTitle','off','UserData',struct('mainFigure',mainFigHandle)); %Make sure that we still know where the userData is kept
            slmAxes=axes('Parent',slmFigure);
            set(slmFigure,'Position',getStatus('slmFigurePosition',[0 0 200 150]),'Resize','on','ResizeFcn',@(obj,event) updateStatus('slmFigurePosition',get(obj,'Position')));
            setUserData('slmAxes',slmAxes);
        end
        img=image(zeros(300/2,400/2,3),'Parent',slmAxes);
        slmDisplayNumber=get(img,'Parent');
    else
        if (strcmpi(descriptors{slmDisplayNumber},'popup XL'))
            slmAxes=getUserData('slmAxes');
            close(get(slmAxes,'Parent'));
            slmAxes=[];
            if (isempty(slmAxes) || ~ishandle(slmAxes))
                mainFigHandle=getBaseFigureHandle();
                slmFigure=figure('Name','Dummy SLM XL','NumberTitle','off','UserData',struct('mainFigure',mainFigHandle)); %Make sure that we still know where the userData is kept
                slmAxes=axes('Parent',slmFigure);
                set(slmFigure,'Position',getStatus('slmFigurePosition',[0 0 200 150]),'Resize','on','ResizeFcn',@(obj,event) updateStatus('slmFigurePosition',get(obj,'Position')));
                setUserData('slmAxes',slmAxes);
            end
            img=image(zeros(300,400,3),'Parent',slmAxes);
            slmDisplayNumber=get(img,'Parent');
        else
            slmDisplayNumber=sscanf(descriptors{slmDisplayNumber},'%u'); % A full-screen, not a pop-up
        end
    end
    setUserData('slmDisplayNumber',slmDisplayNumber);
    
    populateSLMDisplayNumberPopupMenu(); %Update the list in case more were connected
    
    if (~isempty(slm))
        changeSLM();
    end
    
    updateStatus('slmDisplayNumber',get(slmDisplayNumberPopUpMenu,'Value'));
end
function populateSLMDisplayNumberPopupMenu(slmDisplayNumberPopUpMenu)
    if (nargin<1)
        slmDisplayNumberPopUpMenu=getUserData('slmDisplayNumberPopUpMenu');
    end 
    % Find out the number of displays and some properties
    fullScreenManager=getUserData('fullScreenManager');
    optionsList={'popup',fullScreenManager.screens.description}; %,'popup XL'

    set(slmDisplayNumberPopUpMenu,'String',optionsList);
    setUserData('slmDisplayNumberPopUpMenu',slmDisplayNumberPopUpMenu);
end
function changeSLM(obj,event)
    if (nargin<1)
        obj=getUserData('changeSLMPopUpMenu');
    end 
        
    slmDisplayNumber=getUserData('slmDisplayNumber');
%     slmAxes=getUserData('slmAxes');
%     if (~isempty(slmAxes) && slmDisplayNumber~=slmAxes && ishandle(slmAxes))
%         slmFigure=get(slmAxes,'Parent');
%         updateStatus('slmFigurePosition',get(slmFigure,'Position'));
%         close(slmFigure);
%         setUserData('slmAxes',[]);
%     else
%         getUserData('fullScreenManager').close('all');
%     end
    
    slmTypeIndex=get(obj,'Value');
    switch (slmTypeIndex)
        case 1
            slm=PhaseSLM(slmDisplayNumber);
        case 2
            slm=DualHeadSLM(slmDisplayNumber);
    end
    updateStatus('slmTypeIndex',slmTypeIndex);
    if (ishandle(slmDisplayNumber) && strcmpi(get(slmDisplayNumber,'Type'),'axes'))
        slm.stabilizationTime=0.01; % For testing with the popup SLM
    else
        slm.stabilizationTime=0.20;
    end
    % Adapt the SLM so we can know what it is doing when using the DummyCam
    slm.modulatePostFunctor=@(complexModulation,currentImageOnSLM) setUserData('currentImageOnSLM',currentImageOnSLM);
    setUserData('slm',slm);
    
    adjustDeflection();
    adjustTwoPiEquivalent();
    adjustCorrection();
    updateSLMDisplayAndOutput();
end
% Called by the pupil function text input
function updateSLMDisplayAndOutput(obj,event)
    slm=getUserData('slm');
    maskEdit=getUserData('maskEdit');
    
    pupilEquation=get(maskEdit,'String');
    
    [pupilEquationFunctor argumentsUsed]=parsePupilEquation(pupilEquation);
        
    setUserData('pupilFunction',pupilEquationFunctor);
    setUserData('pupilFunctionTimeDependent',argumentsUsed(3))
    setUserData('lastPupilFunctionUpdate',-1);
    setUserData('initialTimeOfSLM',clock()); %Reset clock
    
    updateStatus('mask',get(maskEdit,'String'));
end
% Called by the buttons
function filledInString = insertPupilRadiusAndUpdateSLM(stringWithPupilRadius)
    slm=getUserData('slm');
    pupilRadius=min(slm.regionOfInterest(3:4))/2;
    filledInString=regexprep(stringWithPupilRadius,'(pupil|max)R(ad(ius)?)?',sprintf('%0.0f',pupilRadius));
    
    maskEdit=getUserData('maskEdit'); 
    set(maskEdit,'String',filledInString);
    setUserData('maskEdit',maskEdit); 
    
    updateSLMDisplayAndOutput();
end


function img=calcDummyImage(imgSize,cam,nbFrames)
    wellDepth=40e3;
    darkPhotonElectrons=100;
    photoElectronsPerGraylevel=wellDepth/(2^cam.bitsPerPixel);
    
    slm=getUserData('slm');
    if (isempty(getUserData('currentImageOnSLM')))
        pupilFunction=getUserData('pupilFunction');
        lastPupilFunctionUpdate=getUserData('lastPupilFunctionUpdate');
        [X,Y]=ndgrid([1:slm.regionOfInterest(3)]-floor(slm.regionOfInterest(3)/2)-1,[1:slm.regionOfInterest(4)]-floor(slm.regionOfInterest(4)/2)-1);
        pupil=pupilFunction(Y,X,lastPupilFunctionUpdate);
        pupil=pupil.*slm.referenceDeflection.*slm.correctionFunction;
    else
        currentImageOnSLM=getUserData('currentImageOnSLM');
        pupil=currentImageOnSLM(:,:,3).*exp(2i*pi*currentImageOnSLM(:,:,2)); %Should work for both phase and amplitude SLMs
    end
    % if the SLM ROI is really small:
    if (any(size(pupil)<[256 256]))
        %pupil=pupil(floor(1:.5:(end+.5)),floor(1:.5:(end+.5))); %assume square pixels and interpolate
        pupil(256,256)=0; % assume dark outside
    end
    %Resize to the imgSize
    if (any(imgSize>size(pupil)))
        pupil(imgSize(1),imgSize(2))=0;
    end
    %Suppose the camera pixels are square, then the pupil must be square
    if (size(pupil,1)<size(pupil,2))
        pupil(size(pupil,2),1)=0;
    elseif (size(pupil,2)<size(pupil,1))
        pupil(1,size(pupil,1))=0;
    end
    
    % Calculate the simulated image with a Fourier transform
    frm=numel(pupil)*abs(fftshift(ifft2(ifftshift(pupil(end:-1:1,end:-1:1)))).^2)/numel(pupil);
    
    %Crop if bigger than imgSize
    frm=frm(floor((size(frm,1)-imgSize(1))/2)+[1:imgSize(1)],floor((size(frm,2)-imgSize(2))/2)+[1:imgSize(2)]);
    
    %Pad if ROI larger than the calculated frame
    sizeDifference=cam.regionOfInterest(3:4)-size(pupil);
    if (any(sizeDifference>0))
        frm(cam.regionOfInterest(3),cam.regionOfInterest(4))=0;
        frm=circshift(frm,floor(max(0,1+sizeDifference)./2));
    end
    % Restrict to ROI
    frm=frm(cam.regionOfInterest(1)+[1:cam.regionOfInterest(3)],cam.regionOfInterest(2)+[1:cam.regionOfInterest(4)],:);
    
    %Simulate noise
    frm=frm*wellDepth; %Convert to photo electrons
%     frm=frm+darkPhotonElectrons; %Add dark noise
    frm=frm*cam.gain;
    integrationTimeUnits=nbFrames*cam.integrationTime/20e-3; %Assume that the laser power is adjusted for an integration time of 20ms
    frm=frm*integrationTimeUnits;
%     frm=frm+sqrt(frm*integrationTimeUnits).*randn(size(frm));
%     frm=floor(frm/photoElectronsPerGraylevel);

    img=frm/nbFrames;

    %Normalize the maximum graylevel to 1
    img=min(1,img./(2^cam.bitsPerPixel-1));
end

function closeApp(obj,event)
    setUserData('isClosing',true);
    if (~getUserData('recording'))
        shutDown();
    end %else wait for the loop to exit and clean itself up
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
    value=userData.(fieldName);
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
function fig=getBaseFigureHandle()
    fig=gcbf();
    if (isempty(fig))
        fig=gcf;
    end
    userData=get(fig,'UserData');
    if (isfield(userData,'mainFigure'))
        fig=userData.mainFigure;
    end
end
%
% 'Status' variables are stored persistently to disk and reloaded next time
%
function loadStatus()
    fullPath=mfilename('fullpath');
    try
        load(strcat(fullPath,'.mat'),'status');
    catch Exc
    end
    if (~exist('status','var'))
        status={};
        status.version=-1;
    end
    
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
    
    fullPath=mfilename('fullpath');
    save(strcat(fullPath,'.mat'),'status');
end

function estimateCorrection(obj,event)
    oldValue=get(obj,'String');
    set(obj,'String','Stop!');
    % Locally defined functions to be passed into the aberration correction algorithm
    amplificationLimit=getStatus('amplificationLimit');
    cam=getUserData('cam');
    slm=getUserData('slm');
    if (~getUserData('probing'))
        [cam centerPos]=selectRegionOfInterestAroundPeakIntensity(cam,slm,min(cam.regionOfInterest(3:4),[1 1]*32),[]);
    else
        centerPos=[]; % We are going to stop the probing
        set(obj,'String','Wait....');
    end
    camRegionOfInterest=cam.regionOfInterest;
    function value=probeFunctor() %combinedDeflection)
        probeSize=getUserData('probeSize');
        
        %Make sure that nobody changes the region of interest of the camera!
        if (~all(cam.regionOfInterest~=camRegionOfInterest))
            cam.regionOfInterest=camRegionOfInterest;
        end
        
        img=cam.acquire();
        
        if (size(img,1)<probeSize(1))
            img(probeSize(1),:)=mean(img,1);
        end
        if (size(img,2)<probeSize(2))
            img(:,probeSize(2))=mean(img,2);
        end
        vals=img(centerPos(1)-cam.regionOfInterest(1)+[-floor(probeSize(1)/2):floor((probeSize(1)-1)/2)],centerPos(2)-cam.regionOfInterest(2)+[-floor(probeSize(1)/2):floor((probeSize(1)-1)/2)]);
        value=mean(vals(:));
    end
    prevFractionDone=0;
    function cont=progressFunctor(fractionDone,currentPupilFunctionEstimate)
        cont=true;
        if (floor(fractionDone*100)>floor(prevFractionDone*100))
            logMessage('%0.0f%% done.',100*fractionDone);
            prevFractionDone=fractionDone;
            currentCorrectionEstimate=calcCorrectionFromPupilFunction(currentPupilFunctionEstimate,amplificationLimit);
            slm.modulate(currentCorrectionEstimate);
            showImage(currentPupilFunctionEstimate+.00001i,-1,[],[],getUserData('pupilAxes'));
            if (getUserData('recording'))
                acquireAndDisplay();
            end
            if (getUserData('isClosing') || getUserData('probingInterupted'))
                cont=false;
            end
        end
    end
    if (~getUserData('probing'))
        setUserData('probingInterupted',false);
        setUserData('probing',true);

        %
        % Determine the file where to store the data
        %
        correctionEdit=getUserData('correctionEdit');
        calibrationFileName=get(correctionEdit,'String');
        if (isempty(strtrim(calibrationFileName)))
            calibrationFileName=strcat(pwd(),'/calibrateSetup_',datestr(now(),'YYYY-mm-DD'),'.mat');
            set(correctionEdit,'String',calibrationFileName);
        end

        %
        % Determine the aberration and the correction
        %
        slm=getUserData('slm');
        code=get(getUserData('probeMethodToggle'),'String'); code=upper(code(1));
        switch(code)
            case 'C'
                probeGridSize=[12 16 8]; % [15 20 8]
                logMessage('Using Tom Cizmar''s aberration measurement method with %dx%d probes and %d phases.\nThe region of interest of the spatial light modulator is %dx%d and the probe size is %dx%d.',[probeGridSize, slm.regionOfInterest(3:4) floor(slm.regionOfInterest(3:4)./probeGridSize(1:2))]);
                measuredPupilFunction=aberrationMeasurementCizmarMethod(slm,probeGridSize,@probeFunctor,@progressFunctor);
            case 'Z'
                zernikeCoefficientIndexes=[2 3  4  5 6  7 8  9 10  11];
                logMessage('Using Zernike wavefront measurement method with %d coefficients.',[size(zernikeCoefficientIndexes,2)]);
                measuredPupilFunction=aberrationMeasurementZernikeWavefront(slm,zernikeCoefficientIndexes,@probeFunctor,@progressFunctor);
            otherwise
                probeGridSize=[25 25 3];
                %probeGridSize=[12 12 3]; %test with smaller grid size; Mingzhou
                logMessage('Using Michael Mazilu''s aberration measurement method with a maximum of %dx%d probes and %d phases.',probeGridSize);
                measuredPupilFunction=aberrationMeasurement(slm,probeGridSize,@probeFunctor,@progressFunctor);
        end
        logMessage('Calculating correction function for a maximum amplitude reduction of %0.3f.',amplificationLimit);
        initialCorrection=slm.correctionFunction;
        pupilFunctionCorrection=calcCorrectionFromPupilFunction(measuredPupilFunction.*conj(initialCorrection),amplificationLimit);
        
        %
        % Store
        %
        %Additional information to be stored in the aberration measurement file
        referenceDeflectionFrequency=slm.referenceDeflectionFrequency;
        slmRegionOfInterest=slm.regionOfInterest;
        twoPiEquivalent=slm.twoPiEquivalent;

        save(calibrationFileName,'referenceDeflectionFrequency','slmRegionOfInterest','twoPiEquivalent','initialCorrection','measuredPupilFunction','pupilFunctionCorrection','amplificationLimit','centerPos');
        set(obj,'String',oldValue);
        setUserData('probing',false);

        adjustCorrection();
    else
        setUserData('probingInterupted',true);
    end
    
end
function toggleProbeMethod(obj,event)
    toggleButton=getUserData('probeMethodToggle');
    code=get(toggleButton,'String'); code=upper(code(1));
    switch(code)
        case 'C'
            set(toggleButton,'String','Z');
        case 'Z'
            set(toggleButton,'String','M');
        otherwise
            set(toggleButton,'String','C');
    end
    setUserData('probeMethodToggle',toggleButton);
    updateStatus('probeMethodToggle',get(toggleButton,'String'));
end

function estimateTwoPiEquivalentCallback(obj,event)
%     function cont=progressFunctor(fractionDone,img)
%         cont=true;
%         logMessage('%0.1f%% done.',100*fractionDone);
%         if (floor(fractionDone*100)~=floor((fractionDone*100-1)))
%             if (getUserData('recording'))
%                 displayOnMonitor(img);
%                 if (getUserData('isClosing'))
%                     cont=false;
%                 end
%             end
%         end
%     end
    if (~getUserData('probing'))
        setUserData('probing',true);

        %
        % Determine the file where to store the data
        %
        correctionEdit=getUserData('correctionEdit');
        calibrationFileName=get(correctionEdit,'String');
        if (isempty(strtrim(calibrationFileName)))
            calibrationFileName=strcat(pwd(),'/calibrateSetup_',datestr(now(),'YYYY-mm-DD'),'.mat');
            set(correctionEdit,'String',calibrationFileName);
        end

        %
        % Determine the two-pi-phase equivalent
        %
        cam=getUserData('cam');
        slm=getUserData('slm');
        
        [cam centerPos]=selectRegionOfInterestAroundPeakIntensity(cam,slm,min(cam.regionOfInterest(3:4),[1 1]*128),[]);

        try
            % Determine the phase shift with gray level
            %[phases,graylevels,values,slm,cam]=calibratePhase(slm,cam,@progressFunctor);
            [twoPiEquivalent slm]=estimateTwoPiEquivalent(slm,cam,centerPos,[1 1]*5);
            if (exist(calibrationFileName,'file'))
                save(calibrationFileName,'twoPiEquivalent','-append');
            else
                save(calibrationFileName,'twoPiEquivalent');
            end

            set(getUserData('twoPiEquivalentEdit'),'String',sprintf('%1.0f',twoPiEquivalent*100));

            setUserData('slm',slm);
            setUserData('cam',cam);
        catch Exc
            % Probably user interupted, just stop
        end

        refreshPupilModulationOnSLM(); %Reset SLM
        setUserData('probing',false);
    else
        logMessage('Already probing, wait until correct measurement is done.');
    end
end

function toggleOutputRecording(obj,event)
	outputFileEdit=getUserData('outputFileEdit');
    outputFileName=get(outputFileEdit,'String');
    outputFileName=strtrim(outputFileName);
    
    % Local function, only used to record movies to .mat files
    function appendFrame(img)
        recordingObj=getUserData('recordingObj');
        if (~isempty(recordingObj))
            recordingObj(:,:,end+1)=img;
        else
            recordingObj=img;
        end
        setUserData('recordingObj',recordingObj);
    end

    if (~isempty(outputFileName))
        cam=getUserData('cam');

        fileExtension=lower(outputFileName(end-3:end));
        switch (fileExtension)
            case {'.mat','.avi'}
                if (get(obj,'Value')>0)
                    logMessage('Starting recording...');
                    if (strcmp(fileExtension,'.avi'))
                        try
                            recordingObj = VideoWriter(outputFileName,'Uncompressed AVI'); %'Motion JPEG AVI'); %'Uncompressed AVI');
                            recordingObj.FrameRate=25;
        %                     recordingObj.Quality=75;
                            open(recordingObj);

                            cam.acquisitionFunctor=@(img) writeVideo(recordingObj,min(1,max(0,img))); 
                        catch Exc
                            recordingObj = VideoWriter(outputFileName,'FrameRate',25);

                            cam.acquisitionFunctor=@(img) addframe(recordingObj,min(1,max(0,img)));
                        end
                    else
                        recordingObj=[];
                        cam.acquisitionFunctor=@appendFrame;
                    end
    
                    setUserData('recordingObj',recordingObj);
                else
                    cam.acquisitionFunctor=[];
                    recordingObj=getUserData('recordingObj');
                    if (~isnumeric(recordingObj))
                        close(recordingObj);
                    else
                        save(outputFileName,'recordingObj');
                    end
                    logMessage('Stopped the recording.');
                end
            otherwise
                %Add an extension if not provided
                if (~strcmp(fileExtension,'.png') && ~strcmp(fileExtension,'.tif') && ~strcmp(fileExtension,'.bmp'))
                    outputFileName=strcat(outputFileName,'.png');
                end
                
                imwrite(min(1,max(0,cam.acquire())),outputFileName);
                logMessage('Took snapshot.');
                
                set(obj,'Value',0);
        end

        setUserData('cam',cam);
    end
end

function cameras=detectCameras()
    cameras=[];
    cameras(1).type='DummyCam';
    cameras(1).index=1;
    cameras(1).description='Dummy Cam';
    cameras(2).type='DummyCam';
    cameras(2).index=2;
    cameras(2).description='Dummy Cam HD';
    
    hwInfo=imaqhwinfo();
    if (any(strcmpi(hwInfo.InstalledAdaptors,'gige')))
    	hwInfoGigE=imaqhwinfo('gige');
        for (camIdx=1:length(hwInfoGigE.DeviceIDs))
            cameras(end+1).type='BaslerGigE';
            cameras(end).index=hwInfoGigE.DeviceIDs{camIdx};
            cameras(end).description='Basler GigE';
            if (length(hwInfoGigE.DeviceIDs)>1)
                cameras(end).description=strcat(cameras(end).description,sprintf(' %.0f',hwInfoGigE.DeviceIDs{camIdx}));
            end
        end
    end
    if (any(strcmpi(hwInfo.InstalledAdaptors,'andor')))
    	hwInfoAndor=imaqhwinfo('andor');
        for (camIdx=1:length(hwInfoAndor.DeviceIDs))
            cameras(end+1).type='Andor';
            cameras(end).index=hwInfoAndor.DeviceIDs{camIdx};
            cameras(end).description='Andor';
            if (length(hwInfoGigE.DeviceIDs)>1)
                cameras(end).description=strcat(cameras(end).description,sprintf(' %.0f',hwInfoAndor.DeviceIDs{camIdx}));
            end
        end
    end
    if (any(strcmpi(hwInfo.InstalledAdaptors,'winvideo')))
    	hwInfoIC=imaqhwinfo('winvideo');
        for (camIdx=1:length(hwInfoIC.DeviceIDs))
            cameras(end+1).type='ImagingSource';
            cameras(end).index=hwInfoIC.DeviceIDs{camIdx};
            cameras(end).description='Imaging Source';
            if (length(hwInfoIC.DeviceIDs)>1)
                cameras(end).description=strcat(cameras(end).description,sprintf(' %.0f',hwInfoIC.DeviceIDs{camIdx}));
            end
        end
    end
    if (any(strcmpi(hwInfo.InstalledAdaptors,'avtmatlabadaptor64_r2009b')))
    	hwInfoPC=imaqhwinfo('avtmatlabadaptor64_r2009b');
        for (camIdx=1:length(hwInfoPC.DeviceIDs))
            cameras(end+1).type='PikeCam';
            cameras(end).index=hwInfoPC.DeviceIDs{camIdx};
            cameras(end).description='PikeCam';
            if (length(hwInfoPC.DeviceIDs)>1)
                cameras(end).description=strcat(cameras(end).description,sprintf(' %.0f',hwInfoPC.DeviceIDs{camIdx}));
            end
        end
    end
    if (exist('vcapg2'))
        numberOfCards=vcapg2();
        for cardIdx=1:numberOfCards,
            cameras(end+1).type='DirectShow';
            cameras(end).index=cardIdx;
            cameras(end).description='Win Direct Show';
            if (numberOfCards>1)
                cameras(end).description=strcat(cameras(end).description,sprintf(' %.0f',cardIdx));
            end
        end
    end
end

function str=formatNumberForTickLabel(number,numberForPrec)
    digitsAfterDot=max(0,ceil(-log10(numberForPrec)));
    str=sprintf(sprintf('%%0.%if',digitsAfterDot),number);
end
