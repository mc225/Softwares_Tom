function displayDataCubeSlices(folderNames)
    close all;
    if (nargin<1 || isempty(folderNames))
%         folderNames='Z:\RESULTS\2013-05-14_cellspheroidsWGA\2013-05-14_cellspheroidsWGA_nofilter_ap1-3\cropped';  % Cell cluster section at -27.5um
% %         folderNames='Z:\RESULTS\2013-05-14_cellspheroidsWGA\2013-05-14_cellspheroidsWGA_filter_ap1-3'; %,'Z:\RESULTS\2013-05-14_cellspheroidsWGA_nofilter_ap1-3'}; % Cell cluster section at -32um
%         folderNames='Z:\RESULTS\2013-05-14_cellspheroidsWGA\2013-05-14_cellspheroidsWGA_filter_ap1-2B'; % 5-cell cluster

        %folderNames={'Z:\RESULTS\2013-04-26HEK_ApOneThird\2013-04-26 13_03_52.075_blue1','Z:\RESULTS\2013-04-26HEK_ApOneThird\2013-04-26 12_49_28.119_green1'};
        %folderNames={'Z:\RESULTS\2013-05-09_HEK_transientTransfected_apOneThird\2013-05-09 16_15_24.228_filter_nodichr','Z:\RESULTS\2013-05-09_HEK_transientTransfected_apOneThird\2013-05-09 16_37_05.814_nofilter_dichroic'};
%         folderNames='Z:\RESULTS\2013-05-09_HEK_transientTransfected_apOneThird\2013-05-09 16_37_05.814_nofilter_dichroic';
%         folderNames='Z:\RESULTS\2013-05-14_cellspheroidsWGA\2013-05-14_cellspheroidsWGA_filter_ap1-3'; %,'Z:\RESULTS\2013-05-14_cellspheroidsWGA_nofilter_ap1-3'}; % Cell cluster section at -35um
%         folderNames='Z:\RESULTS\2013-05-14_cellspheroidsWGA\2013-05-14_cellspheroidsWGA_filter_ap1-2B'; % 5-cell cluster
%         folderNames='Z:\RESULTS\2013-05-14_cellspheroidsWGA\2013-05-14_cellspheroidsWGA_filter_ap1-3B'; % 5 cell cluster (not used?)
%         folderNames={'Z:\RESULTS\2013-05-15\2013-05-15_cellspheroidsWGA_apOneHalf\2013-05-15 17_31_34.698_A','Z:\RESULTS\2013-05-15\2013-05-15_cellspheroidsWGA_apOneHalf\2013-05-15 17_45_03.053_A'};
%         folderNames={'Z:\RESULTS\2013-05-15\2013-05-15_HEKsWGA_apOneHalf\2013-05-15 16_52_04.156','Z:\RESULTS\2013-05-15\2013-05-15_HEKsWGA_apOneHalf\2013-05-15 17_11_39.629'};
%         folderNames={'Z:\RESULTS\2013-05-15\2013-05-15_HEKWGA_apOneHalf\2013-05-15 19_16_41.859_A','Z:\RESULTS\2013-05-15\2013-05-15_HEKWGA_apOneHalf\2013-05-15 19_26_28.113_A'};
%         folderNames={'Z:\RESULTS\2013-05-15\2013-05-15_HEKWGA_apOneHalf\2013-05-15 19_43_10.533_B','Z:\RESULTS\2013-05-15\2013-05-15_HEKWGA_apOneHalf\2013-05-15 19_52_04.788_B_changedLambda'};
%         folderNames={'Z:\RESULTS\2013-05-15\2013-05-15_HEKWGA_apOneThird\2013-05-15 18_54_03.883','Z:\RESULTS\2013-05-15\2013-05-15_HEKWGA_apOneThird\2013-05-15 19_04_11.913'};
%         folderNames={'Z:\RESULTS\2013-05-15\2013-05-15_cellspheroidsWGA_apOneThird\2013-05-15 18_12_45.370','Z:\RESULTS\2013-05-15\2013-05-15_cellspheroidsWGA_apOneThird\2013-05-15 18_23_14.789_changedColor'}; % Dual color

%         folderNames={'C:\Users\Tom\Dropbox\AirySheetData\DataCubes\2013-05-14_cellspheroidsWGA\2013-05-14_cellspheroidsWGA_filter_ap1-2B'}; % 5 cells

%         folderNames={'E:\Stored Files\RESULTS\2013-08-08\HEK2\2013-08-08 15_57_16.566\'};
%         folderNames={'Z:\RESULTS\2013-09-09\2013-09-09 17_20_19.909 at0um_sliced'};
%         folderNames={'Z:\RESULTS\2013-09-09\2013-09-09 17_28_37.355 at_50um_sliced'};
%         folderNames={'Z:\RESULTS\2013-09-09\2013-09-09 17_38_06.531 at_85um_sliced'};
%         folderNames={'Z:\RESULTS\2013-09-13_Amphioxus\Run_1\2013-09-13 12_31_28.686'}; % 488 nm 
%         folderNames={'Z:\RESULTS\2013-09-13_Amphioxus\Run_1\2013-09-13 12_34_55.726'}; % 532 nm
%         folderNames={'Z:\RESULTS\2013-09-13_Amphioxus\Run_2\2013-09-13 12_44_38.570'}; % 488 nm 
%         folderNames={'Z:\RESULTS\2013-09-13_Amphioxus\Run_2\2013-09-13 12_47_53.135'}; % 532 nm
%         folderNames={'Z:\RESULTS\2013-10-10_200nmBeadsNearEntranceInSmallAgaroseDroplet\sample5ApLarge\2013-10-10 18_41_27.129\'};
%         folderNames={'Z:\RESULTS\2013-10-10_200nmBeadsNearEntranceInSmallAgaroseDroplet\sample5ApSmall\2013-10-10 18_16_34.014\'};
%         folderNames={'Z:\RESULTS\600nmBeadsApOneThird\B\at0\'};
%         folderNames={'Z:\RESULTS\2013-01-17 16_59_54.447_Green_2W_Andor_0.05um_(1)\'};
%         folderNames={'Z:\RESULTS\2013-01-17 17_07_40.377_Blue_30A_Andor_0.2um_(1)'};
%         folderNames={'Z:\RESULTS\2013-04-10_Beads200nm\2013-04-10 16_27_33.222_ApOneQuarter'};
%         folderNames={'Z:\RESULTS\2013-04-10_Beads200nm\2013-04-10 16_10_20.097_ApOneHalf'};
%         folderNames={'Z:\RESULTS\20120501_alpha7_int1000_g1_morepower_beads_leftofprobe\'};
%         folderNames={'Z:\RESULTS\20120501_alpha7_int1000_g1_littlemorepower_betterfocus_beads_leftofprobe\forwardScans\'};
%        folderNames={'Z:\RESULTS\20120501_alpha7_int1000_g1_littlemorepower_betterfocus_beads_leftofprobe\forwardScans\'};
       folderNames={'H:\RESULTS\2013-09-13_Amphioxus\Run_1\2013-09-13 12_31_28.686'}; % 488 nm
%         folderNames={'Z:\RESULTS\2013-09-13_Amphioxus_colorCombined\Run_2'};
    end
    
    if (~iscell(folderNames))
        folderNames={folderNames};
    end
    
    % Load the default values
    fig=figure('Visible','off','CloseRequestFcn',@(obj,event) shutDown(obj));
    getBaseFigureHandle(fig); % Set permanent reference
    loadStatus([folderNames{1},'/displayDataCubeSlices.mat']);
    set(fig,'ResizeFcn',@(obj,event) updateStatus('windowPosition',get(obj,'Position')));
    updateStatus('version',0.10);
    
    dataCubeNames={};
    for (folderName=folderNames)
        folderName=folderName{1};
        namesInStructs=dir([folderName,'/*.mat']);
        if (~isempty(namesInStructs))
            for fileName={namesInStructs.name}
                fileName=fileName{1};
                dataCubeNames{end+1}=[folderName,'/',fileName];
            end
        else
            logMessage('No files found in folder %s.',folderName);
        end
    end
    
    % Make a list of available data sorted by alpha value
    fileDescriptors=struct(); % Use a struct as a hash table
    for (fullFileName=dataCubeNames)
        fullFileName=fullFileName{1};
        try
            % Determine the wavelength
            wavelength=regexpi(fullFileName,'lambda([\d]*)nm_alpha[^\\/]+','tokens');
            if (~isempty(wavelength))
                wavelength=str2double(wavelength{1})*1e-9;
                % Determine the light sheet type
                alpha=regexpi(fullFileName,'alpha([\d]*)_[^\\/]+','tokens');
                alpha=str2double(alpha{1});
                beta=regexpi(fullFileName,'beta([\d]*)[\D\s]','tokens');
                beta=str2double(beta{1})/100;
            else
                configFile=strcat(fullFileName(1:end-4),'.json');
                specificConfig=loadjson(configFile);
                alpha=specificConfig.modulation.alpha;
                beta=specificConfig.modulation.beta;
                wavelength=532e-9;
            end
            fileDescriptorTag=sprintf('alpha%06.0f_beta%03.0f',[abs(alpha)*1000 beta*1000]);
            if (isfield(fileDescriptors,fileDescriptorTag))
                fileDescriptor=getfield(fileDescriptors,fileDescriptorTag);
            else
                fileDescriptor={};
                fileDescriptor.fullFileNameGreen=[];
                fileDescriptor.fullFileNameRed=[];
            end
            if (abs(wavelength-488e-9)<20e-9)
                fileDescriptor.fullFileNameGreen=fullFileName;
            else
                fileDescriptor.fullFileNameRed=fullFileName;
            end
            fileDescriptor.alpha=alpha;
            fileDescriptor.beta=beta;
            if (fileDescriptor.alpha==0)
                if (fileDescriptor.beta==1)
                    lightSheetType='Gaussian';
                else
                    lightSheetType=sprintf('Bessel %0.0f%%',fileDescriptor.beta*100);
                end
            else
                if (fileDescriptor.beta==1)
                    lightSheetType=sprintf('Airy %0.1f',fileDescriptor.alpha);
                else
                    lightSheetType=sprintf('a=%0.1f,b=%0.0f%%',[fileDescriptor.alpha fileDescriptor.beta*100]);
                end
            end
            fileDescriptor.lightSheetType=lightSheetType;

            fileDescriptors=setfield(fileDescriptors,fileDescriptorTag,fileDescriptor);
        catch Exc
            logMessage('Skipping %s.',fullFileName);
        end
    end
    % Convert hash table to a sorted array of structs
    [~, sortI]=sort(fieldnames(fileDescriptors));
    fileDescriptors=cell2mat(struct2cell(fileDescriptors));
    fileDescriptors=fileDescriptors(sortI);
    
    % Create the dynamic figure
    drawFigure(fig,fileDescriptors,folderNames{1});
end

function drawFigure(fig,fileDescriptors,description)
    nbPanels=2;
    % Open the figure window
    windowMargins=[1 1]*128;
    initialPosition=get(0,'ScreenSize')+[windowMargins -2*windowMargins];
    set(fig,'Name',description,'Position',getStatus('windowPosition',initialPosition),'Visible','on');
    panelSpacing=0*[1/30 1/30]; % Of complete box
    axesSpacing=[1/20 1/7]; % Of complete box
    projectionsPanel=uipanel('Parent',fig,'BackgroundColor',[1 1 1],'Units','normalized','Position',[0.15 0 0.85 1]);
    lightSheetTypes={'<select>', fileDescriptors.lightSheetType};
    setUserData('fileDescriptors',fileDescriptors);
    panelInfos=[];
    for panelIdx=1:nbPanels,
        panelInfo=struct();
        
        panelInfo.panel=uipanel('Parent',projectionsPanel,'BackgroundColor',[1 1 1],'BorderWidth',0,'Units','normalized','Position',[panelSpacing(1) panelSpacing(2)+(nbPanels-panelIdx)/nbPanels 1-panelSpacing(1) 1/nbPanels-panelSpacing(2)]);
        panelInfo.createMovieButton=uicontrol('Parent',panelInfo.panel,'Position',[10 130 80 20], 'Style','pushbutton', 'String','Create Movie', 'Callback', {@createMovie,fig,panelIdx});
        panelInfo.rawDeconvPopup=uicontrol('Parent',panelInfo.panel,'Position',[10 100 80 20], 'Style','popupmenu', 'String', {'raw','deconvolved'}, 'Value',getStatus(sprintf('rawDeconvPopupValue%i',panelIdx),2), 'Callback', {@reloadData,fig});
%         beamTypePopupValue=find(strcmp(lightSheetTypes,getStatus(sprintf('beamTypePopupString%i',panelIdx),'<select>')),1,'first');
        beamTypePopupValue=1; % Don't preselect
        if (isempty(beamTypePopupValue))
            beamTypePopupValue=1;
        end
        panelInfo.beamTypePopup=uicontrol('Parent',panelInfo.panel,'Position',[10 70 80 20], 'Style','popupmenu', 'String',lightSheetTypes, 'Value',beamTypePopupValue, 'Callback', {@reloadData,fig});
        panelInfo.axes=[];
        panelInfo.axes(1)=axes('Parent',panelInfo.panel,'Units','normalized','Position',[1/6+axesSpacing(1) axesSpacing(2) 2/6-axesSpacing(1) 1-axesSpacing(2)]);
        panelInfo.axes(2)=axes('Parent',panelInfo.panel,'Units','normalized','Position',[3/6+axesSpacing(1) axesSpacing(2) 2/6-axesSpacing(1) 1-axesSpacing(2)]);
        panelInfo.axes(3)=axes('Parent',panelInfo.panel,'Units','normalized','Position',[5/6+axesSpacing(1) axesSpacing(2) 1/6-axesSpacing(1) 1-axesSpacing(2)]);
        
        if (isempty(panelInfos))
            panelInfos=panelInfo;
        else
            panelInfos(end+1)=panelInfo;
        end
    end
    setUserData('panelInfos',panelInfos);
    linkaxes(arrayfun(@(x) x.axes(1),panelInfos));
    linkaxes(arrayfun(@(x) x.axes(2),panelInfos));
    linkaxes(arrayfun(@(x) x.axes(3),panelInfos));
        
    % Draw the controls
    controlsPanel=uipanel('Parent',fig,'Units','normalized','Position',[0 0 0.15 1],'Title','Controls');
    uicontrol('Parent',controlsPanel,'Style','text','String','Colormap:','Position',[10,240,80,20]);
    setUserData('colorMapPopup',uicontrol('Parent',controlsPanel,'Position',[100 240 80 20], 'Style','popupmenu', ...
        'String', {'Screen','Print'}, 'Value',getStatus('colorMapPopupValue',1), 'Callback', {@updateProjections,fig}));
        
    uicontrol('Parent',controlsPanel,'Style','text','String','Background red:','Position',[10,210,100,20]);
    setUserData('redBackgroundEdit',uicontrol('Style','edit','String',getStatus('redBackgroundString','0'),'Position',...
        [100,210,45,20],'Callback',{@updateProjections,fig}));
    uicontrol('Parent',controlsPanel,'Style','text','String','%, green:','Position',[140,210,50,20]);
    setUserData('greenBackgroundEdit',uicontrol('Style','edit','String',getStatus('greenBackgroundString','0'),'Position',...
        [190,210,45,20],'Callback',{@updateProjections,fig}));
    uicontrol('Parent',controlsPanel,'Style','text','String','%','Position',[235,210,20,20]);
    
    uicontrol('Parent',controlsPanel,'Style','text','String','Rel. Int. red:','Position',[10,180,100,20]);
    setUserData('redRelativeDisplayIntensityEdit',uicontrol('Style','edit','String',getStatus('redRelativeDisplayIntensityString','1000'),'Position',...
        [100,180,45,20],'Callback',{@updateProjections,fig}));
    uicontrol('Parent',controlsPanel,'Style','text','String','%, green:','Position',[140,180,50,20]);
    setUserData('greenRelativeDisplayIntensityEdit',uicontrol('Style','edit','String',getStatus('greenRelativeDisplayIntensityString','1000'),'Position',...
        [190,180,45,20],'Callback',{@updateProjections,fig}));
    uicontrol('Parent',controlsPanel,'Style','text','String','%','Position',[235,180,20,20]);
    
    uicontrol('Parent',controlsPanel,'Style','text','String','Green X-shift','Position',[10,150,80,20]);
    setUserData('yShiftEdit',uicontrol('Style','edit','String',getStatus('yShiftString','0'),'Position',...
        [100,150,40,20],'Callback',{@updateDataCubeSection,fig}));
    uicontrol('Parent',controlsPanel,'Style','text','String','Green Y-shift','Position',[10,130,80,20]);
    setUserData('xShiftEdit',uicontrol('Style','edit','String',getStatus('xShiftString','0'),'Position',...
        [100,130,40,20],'Callback',{@updateDataCubeSection,fig}));
    uicontrol('Parent',controlsPanel,'Style','text','String','Green Z-shift','Position',[10,110,80,20]);
    setUserData('zShiftEdit',uicontrol('Style','edit','String',getStatus('zShiftString','0'),'Position',...
        [100,110,40,20],'Callback',{@updateDataCubeSection,fig}));
    % sectioning controls
    uicontrol('Parent',controlsPanel,'Style','text','String','X-stack position','Position',[10,50,80,20]);
    setUserData('yRangeSelEdit',uicontrol('Style','edit','String',getStatus('yRangeSelString',''),'Position',[100,50,90,20],'Callback',{@updateDataCubeSection,fig}));
    uicontrol('Parent',controlsPanel,'Style','text','String','Y-stack position','Position',[10,30,80,20]);
    setUserData('xRangeSelEdit',uicontrol('Style','edit','String',getStatus('xRangeSelString',''),'Position',[100,30,90,20],'Callback',{@updateDataCubeSection,fig}));
    uicontrol('Parent',controlsPanel,'Style','text','String','Z-stack position','Position',[10,10,80,20]);
    setUserData('zRangeSelEdit',uicontrol('Style','edit','String',getStatus('zRangeSelString',''),'Position',[100,10,90,20],'Callback',{@updateDataCubeSection,fig}));
        
    set(fig,'Toolbar','figure');
    drawnow();
    
    reloadData([],[],fig);
end

function createMovie(obj,event,fig,panelIdx)
    
end

function reloadData(obj,event,fig)
    flipPropagationAxis=false;
        
    fig=getBaseFigureHandle();
    fileDescriptors=getUserData('fileDescriptors');
    panelInfos=getUserData('panelInfos');
    
    experimentalSettings=[];
    for panelIdx=1:numel(panelInfos),
        panelInfo=panelInfos(panelIdx);
        beamTypeValue=get(panelInfo.beamTypePopup,'Value')-1; % skip <select>
        if (beamTypeValue>0) % skip the <select>
            beamTypes=get(panelInfo.beamTypePopup,'String');
            if (beamTypeValue<=length(beamTypes))
                beamType=beamTypes{beamTypeValue};
                updateStatus(sprintf('beamTypePopupString%i',panelIdx),beamType);
                set(panelInfo.panel,'Title',beamType);
                beamTypeValue=find(strcmp(beamTypes,beamType),1);
                fileDescriptor=fileDescriptors(beamTypeValue);
                if (isempty(experimentalSettings))
                    % Load the experimental settings first
                    if (~isempty(fileDescriptor.fullFileNameGreen))
                        logMessage(['Loading settings from ',fileDescriptor.fullFileNameGreen,'.']);
                        experimentalSettings=load(fileDescriptor.fullFileNameGreen,'xRange','yRange','zRange','tRange','setupConfig');
                    else
                        logMessage(['Loading settings from ',fileDescriptor.fullFileNameRed,'.']);
                        experimentalSettings=load(fileDescriptor.fullFileNameRed,'xRange','yRange','zRange','tRange','setupConfig');
                    end
                    if (flipPropagationAxis)
                        experimentalSettings.yRange=-experimentalSettings.yRange(end:-1:1);
                    end
                    stageTranslationStepSize=norm(median(diff(experimentalSettings.setupConfig.stagePositions.target)));
                    if isempty(experimentalSettings.setupConfig.detector.center)
                        experimentalSettings.setupConfig.detector.center=[0 0];
                    end
%                     if (~isfield(experimentalSettings,'xRange'))
%                         realMagnification=experimentalSettings.setupConfig.detection.objective.magnification*experimentalSettings.setupConfig.detection.tubeLength/experimentalSettings.setupConfig.detection.objective.tubeLength;
%                         experimentalSettings.xRange=-experimentalSettings.setupConfig.detector.center(1)+experimentalSettings.setupConfig.detector.pixelSize(1)*([1:size(recordedImageStack,1)]-floor(size(recordedImageStack,1)/2)-1)/realMagnification; % up/down
%                         experimentalSettings.yRange=-experimentalSettings.setupConfig.detector.center(2)+experimentalSettings.setupConfig.detector.pixelSize(2)*([1:size(recordedImageStack,2)]-floor(size(recordedImageStack,2)/2)-1)/realMagnification; % left/right
%                         experimentalSettings.tRange=stageTranslationStepSize*([1:size(recordedImageStack,3)]-floor(size(recordedImageStack,3)/2+1))/experimentalSettings.setupConfig.detector.framesPerSecond; %Translation range (along z-axis)
%                         experimentalSettings.zRange=tRange; % detection axis, zero centered
%                     end
                    setUserData('xRange',experimentalSettings.xRange);
                    setUserData('yRange',experimentalSettings.yRange);
                    setUserData('zRange',experimentalSettings.zRange);
%                     if (isfield(experimentalSettings,'tRange'))
                    setUserData('tRange',experimentalSettings.tRange);
%                     else
%                         setUserData('tRange',experimentalSettings.zRange);
%                     end
                end
                dropDownValue=get(panelInfo.rawDeconvPopup,'Value');
                updateStatus(sprintf('rawDeconvPopupValue%i',panelIdx),dropDownValue);
                restored=dropDownValue>1;
                if (restored)
                    dataCubeField='restoredDataCube';
                else
                    dataCubeField='recordedImageStack';
                end

                % Free memory
                panelUserData.dataCube=[]; set(panelInfo.panel,'UserData',panelUserData);
                if (~isempty(fileDescriptor.fullFileNameRed))
                    logMessage(['Loading RED channel from ',fileDescriptor.fullFileNameRed,'...']);
                    loadedData=load(fileDescriptor.fullFileNameRed,dataCubeField);
                    panelUserData.dataCube=getfield(loadedData,dataCubeField);
                    clear loadedData;
                    panelUserData.channels={'red'};
                    if (~isempty(fileDescriptor.fullFileNameGreen))
                        logMessage(['Loading GREEN channel from ',fileDescriptor.fullFileNameGreen,'...']);
                        loadedData=load(fileDescriptor.fullFileNameGreen,dataCubeField);
                        panelUserData.dataCube(:,:,:,2)=getfield(loadedData,dataCubeField);
                        clear loadedData;
                        panelUserData.channels{end+1}='green';
                    end
                else
                    logMessage(['Loading SINGLE GREEN channel from ',fileDescriptor.fullFileNameGreen,'...']);
                    loadedData=load(fileDescriptor.fullFileNameGreen,dataCubeField);
                    panelUserData.dataCube=getfield(loadedData,dataCubeField);
                    clear loadedData;
                    panelUserData.channels={'green'};
                end
                if (flipPropagationAxis)
                    for (colIdx=1:floor(size(panelUserData.dataCube,2)/2))
                        panelUserData.dataCube(:,[colIdx end-(colIdx-1)],:,:)=panelUserData.dataCube(:,[end-(colIdx-1) colIdx],:,:);
                    end
                end
            else
                logMessage('No beams found on disk.');
            end
        else
            panelUserData={};
            panelUserData.dataCube=[];
            panelUserData.channels={'black'};
        end % of if beamTypeValue>0
        set(panelInfo.panel,'UserData',panelUserData);
    end % of for
    
    updateDataCubeSection([],[],fig);
end


function fig=updateDataCubeSection(hObj,event,fig)
%     projector=@mean;
    projector=@(X,dim) max(X,[],dim); % Maximum intensity projection
       
    % Interpret the sectioning controls
    xRange=getUserData('xRange');
    yRange=getUserData('yRange');
    zRange=getUserData('zRange');
    tRange=getUserData('tRange');
    if (~isempty(xRange))
        xRangeSel=readDefaultAndConvertToIndexRange(getUserData('xRangeSelEdit'),xRange,xRange([1 end]));
        yRangeSel=readDefaultAndConvertToIndexRange(getUserData('yRangeSelEdit'),yRange,yRange([1 end]));
        zRangeSel=readDefaultAndConvertToIndexRange(getUserData('zRangeSelEdit'),zRange,zRange([1 end]));
        tRangeSel=zRangeSel;
        xShift=readDefaultAndConvertToIndices(getUserData('xShiftEdit'),xRange,0);
        yShift=readDefaultAndConvertToIndices(getUserData('yShiftEdit'),yRange,0);
        zShift=readDefaultAndConvertToIndices(getUserData('zShiftEdit'),zRange,0);
        [~, dataCubeCenter(1)]=min(abs(xRange));
        [~, dataCubeCenter(2)]=min(abs(yRange));
        [~, dataCubeCenter(3)]=min(abs(zRange));
        greenShiftInPixels=[xShift yShift zShift]-dataCubeCenter;

        % Save the settings for next time we run
        updateStatus('xRangeSelString',get(getUserData('xRangeSelEdit'),'String'));
        updateStatus('yRangeSelString',get(getUserData('yRangeSelEdit'),'String'));
        updateStatus('zRangeSelString',get(getUserData('zRangeSelEdit'),'String'));
        updateStatus('xShiftString',get(getUserData('xShiftEdit'),'String'));
        updateStatus('yShiftString',get(getUserData('yShiftEdit'),'String'));
        updateStatus('zShiftString',get(getUserData('zShiftEdit'),'String'));
        
        % Shift the Green channel
        xRangeSelGreen=min(length(xRange),max(1,xRangeSel-greenShiftInPixels(1)));
        yRangeSelGreen=min(length(yRange),max(1,yRangeSel-greenShiftInPixels(2)));
        zRangeSelGreen=min(length(zRange),max(1,zRangeSel-greenShiftInPixels(3)));
            
        for (panelInfo=getUserData('panelInfos'))
            panelUserData=get(panelInfo.panel,'UserData');
            if (~isempty(panelUserData.dataCube))
                % Calculate the projections
                panelUserData.projZY=projector(panelUserData.dataCube(xRangeSel,yRangeSel,zRangeSel,1),1);
                panelUserData.projYX=projector(panelUserData.dataCube(xRangeSel,yRangeSel,zRangeSel,1),3);
                panelUserData.projZX=projector(panelUserData.dataCube(xRangeSel,yRangeSel,zRangeSel,1),2);
                if (ndims(panelUserData.dataCube)>=4)
                    panelUserData.projZY(:,:,:,2)=projector(panelUserData.dataCube(xRangeSelGreen,yRangeSelGreen,zRangeSelGreen,2),1);
                    panelUserData.projYX(:,:,:,2)=projector(panelUserData.dataCube(xRangeSelGreen,yRangeSelGreen,zRangeSelGreen,2),3);
                    panelUserData.projZX(:,:,:,2)=projector(panelUserData.dataCube(xRangeSelGreen,yRangeSelGreen,zRangeSelGreen,2),2);
                end
                panelUserData.projZY=shiftdim(panelUserData.projZY,1); %Remove new singleton dimension
                panelUserData.projZY=permute(panelUserData.projZY,[2 1 3:ndims(panelUserData.projZY)]);
                panelUserData.projYX=permute(panelUserData.projYX,[1 2 4:ndims(panelUserData.projYX) 3]); %Remove new singleton dimension
                panelUserData.projZX=permute(panelUserData.projZX,[1 3:ndims(panelUserData.projZX) 2]);  %Remove new singleton dimension
                panelUserData.projZX=permute(panelUserData.projZX,[2 1 3:ndims(panelUserData.projZX)]);
            else
                % Fake the projections
                panelUserData.projZY=zeros(length(yRangeSel),length(zRangeSel),length(panelUserData.channels));
                panelUserData.projZY=permute(panelUserData.projZY,[2 1 3:ndims(panelUserData.projZY)]);
                panelUserData.projYX=zeros(length(xRangeSel),length(yRangeSel),length(panelUserData.channels));
                panelUserData.projZX=zeros(length(xRangeSel),length(zRangeSel),length(panelUserData.channels));
                panelUserData.projZX=permute(panelUserData.projZX,[2 1 3:ndims(panelUserData.projZX)]);
            end
            
            set(panelInfo.panel,'UserData',panelUserData);
        end
        
        fig=updateProjections(hObj,event,fig);
    end
end
function fig=updateProjections(obj,event,fig)
    colorMapPopupValue=get(getUserData('colorMapPopup'),'Value');
    updateStatus('colorMapPopupValue',colorMapPopupValue);
    redBackground=str2double(get(getUserData('redBackgroundEdit'),'String'));
    if (isempty(redBackground) || isnan(redBackground))
        redBackground=0;
    end
    redBackground=redBackground(1)/100;
    set(getUserData('redBackgroundEdit'),'String',sprintf('%0.3f',redBackground*100));
    updateStatus('redBackgroundString',sprintf('%0.3f',redBackground*100));
    greenBackground=str2double(get(getUserData('greenBackgroundEdit'),'String'));
    if (isempty(greenBackground) || isnan(greenBackground))
        greenBackground=0;
    end
    greenBackground=greenBackground(1)/100;
    set(getUserData('greenBackgroundEdit'),'String',sprintf('%0.3f',greenBackground*100));
    updateStatus('greenBackgroundString',sprintf('%0.3f',greenBackground*100));

    colorsForPrint=colorMapPopupValue==2;
    colorMapGreen=@(values) interpolatedColorMap(values,[0 0 0; 0 1 0; 1 1 1],[0 .5 1]);
    colorMapRed=@(values) interpolatedColorMap(values,[0 0 0; 1 0 0; 1 1 1],[0 .5 1]);
%         colorMapImageSize=[1024 64];
%         colorMapSingleColorImage=mapColor(repmat([1:colorMapImageSize(1)].',[1 colorMapImageSize(2)]),colorMapSingleColor);
%         imwrite('colorMapSingleColorImage.png',colorMapSingleColorImage);

    % Normalize
    redRelativeDisplayIntensity=str2double(get(getUserData('redRelativeDisplayIntensityEdit'),'String'));
    if (isempty(redRelativeDisplayIntensity))
        redRelativeDisplayIntensity=100;
        set(getUserData('redRelativeDisplayIntensityEdit'),'String',sprintf('%0.0f',redRelativeDisplayIntensityString));
    end
    redRelativeDisplayIntensity=redRelativeDisplayIntensity/100;
    updateStatus('redRelativeDisplayIntensityString',get(getUserData('redRelativeDisplayIntensityEdit'),'String'));
    greenRelativeDisplayIntensity=str2double(get(getUserData('greenRelativeDisplayIntensityEdit'),'String'));
    if (isempty(greenRelativeDisplayIntensity))
        greenRelativeDisplayIntensity=100;
        set(getUserData('greenRelativeDisplayIntensityEdit'),'String',sprintf('%0.0f',greenRelativeDisplayIntensityString));
    end
    greenRelativeDisplayIntensity=greenRelativeDisplayIntensity/100;
    updateStatus('greenRelativeDisplayIntensityString',get(getUserData('greenRelativeDisplayIntensityEdit'),'String'));
        
    panelIdx=1;
    for (panelInfo=getUserData('panelInfos'))
        panelUserData=get(panelInfo.panel,'UserData');
        projZY=panelUserData.projZY;
        projYX=panelUserData.projYX;
        projZX=panelUserData.projZX;
        
        % Remove background and normalize
        normalization=ones(1,1,length(panelUserData.channels));
        for channelIdx=1:length(panelUserData.channels)
            if (strcmp(panelUserData.channels{channelIdx},'red'))
                backgroundOffset=redBackground;
                relativeDisplayIntensity=redRelativeDisplayIntensity;
            else
                backgroundOffset=greenBackground;
                relativeDisplayIntensity=greenRelativeDisplayIntensity;
            end
            projZY(:,:,channelIdx)=projZY(:,:,channelIdx)-backgroundOffset;
            projYX(:,:,channelIdx)=projYX(:,:,channelIdx)-backgroundOffset;
            projZX(:,:,channelIdx)=projZX(:,:,channelIdx)-backgroundOffset;
            normalization(channelIdx)=relativeDisplayIntensity;
        end
%         normalization=max([max(max(projZY)), max(max(projYX)), max(max(projZX))]);
%         if (length(normalization)>=2)
%             normalization(2)=normalization(2)/greenRelativeDisplayIntensity;
%         end
        projZY=projZY.*repmat(normalization,[size(projZY,1) size(projZY,2)]);
        projYX=projYX.*repmat(normalization,[size(projYX,1) size(projYX,2)]);
        projZX=projZX.*repmat(normalization,[size(projZX,1) size(projZX,2)]);

        % Add blue channel, or user colormap in case only one channel is present
        if (ndims(projZY)>2)
            projZY=colorMapChannels(projZY,{colorMapRed,colorMapGreen},colorsForPrint);
            projYX=colorMapChannels(projYX,{colorMapRed,colorMapGreen},colorsForPrint);
            projZX=colorMapChannels(projZX,{colorMapRed,colorMapGreen},colorsForPrint);
        else
            if (strcmp(panelUserData.channels{1},'red'))
                colorMapSingleColor=colorMapRed;
            else
                colorMapSingleColor=colorMapGreen;
            end
            projZY=colorMapChannels(projZY,colorMapSingleColor,colorsForPrint);
            projYX=colorMapChannels(projYX,colorMapSingleColor,colorsForPrint);
            projZX=colorMapChannels(projZX,colorMapSingleColor,colorsForPrint);
        end

        % Draw the projections
        xRange=getUserData('xRange');
        yRange=getUserData('yRange');
        zRange=getUserData('zRange');
        tRange=getUserData('tRange');
        xRangeSel=readDefaultAndConvertToIndexRange(getUserData('xRangeSelEdit'),xRange,xRange([1 end]));
        yRangeSel=readDefaultAndConvertToIndexRange(getUserData('yRangeSelEdit'),yRange,yRange([1 end]));
        zRangeSel=readDefaultAndConvertToIndexRange(getUserData('zRangeSelEdit'),zRange,zRange([1 end]));
        tRangeSel=zRangeSel;
        axesColor=[0 0 0];
        
        showImage(projYX,[],yRange(yRangeSel)*1e6,xRange(xRangeSel)*1e6,panelInfo.axes(1));
        axis equal;
        set(panelInfo.axes(1),'TickDir','out','TickLength',[0.01 0.025]*3,'Color',axesColor,'XColor',axesColor,'YColor',axesColor,'LineWidth',3,'FontSize',16,'FontWeight','bold');
        xlabel('Parent',panelInfo.axes(1),'X (propagation) [\mum]','LineWidth',3,'FontSize',22,'FontWeight','bold');
        ylabel('Parent',panelInfo.axes(1),'Y (swipe) [\mum]','LineWidth',3,'FontSize',22,'FontWeight','bold');
        
        showImage(projZY,[],yRange(yRangeSel)*1e6,zRange(zRangeSel)*1e6,panelInfo.axes(2));
        axis equal;
        set(panelInfo.axes(2),'TickDir','out','TickLength',[0.01 0.025]*3,'Color',axesColor,'XColor',axesColor,'YColor',axesColor,'LineWidth',3,'FontSize',16,'FontWeight','bold');
        xlabel('Parent',panelInfo.axes(2),'X (propagation) [\mum]','LineWidth',3,'FontSize',22,'FontWeight','bold');
        ylabel('Parent',panelInfo.axes(2),'Z (axial) [\mum]','LineWidth',3,'FontSize',22,'FontWeight','bold');
        
        showImage(projZX,[],xRange(xRangeSel)*1e6,tRange(tRangeSel)*1e6,panelInfo.axes(3));
        axis equal;
        set(panelInfo.axes(3),'TickDir','out','TickLength',[0.01 0.025]*3,'Color',axesColor,'XColor',axesColor,'YColor',axesColor,'LineWidth',3,'FontSize',16,'FontWeight','bold');
        xlabel('Parent',panelInfo.axes(3),'Y (swipe) [\mum]','LineWidth',3,'FontSize',22,'FontWeight','bold');
        ylabel('Parent',panelInfo.axes(3),'Z (axial) [\mum]','LineWidth',3,'FontSize',22,'FontWeight','bold');
        
        imwrite(projYX,sprintf('projXY%0.0f.png',panelIdx));
        imwrite(projZY,sprintf('projZX%0.0f.png',panelIdx));
        imwrite(projZX,sprintf('projZY%0.0f.png',panelIdx));
        panelIdx=panelIdx+1;
    end
end
function [rangeSelStart rangeSelEnd]=readDefaultAndConvertToIndices(rangeSelEdit,range,defaultValues)
    rangeSel=str2num(get(rangeSelEdit,'String'))*1e-6;
    if (isempty(rangeSel))
        rangeSel=defaultValues(1);
    end
    if (length(rangeSel)<2 && length(defaultValues)>=2)
        rangeSel(2)=defaultValues(end);
    end
    set(rangeSelEdit,'String',sprintf('%0.3f ',rangeSel*1e6));
    % Convert metric to indexes
    [~, rangeSelStart]=min(abs(range-min(rangeSel)));
    [~, rangeSelEnd]=min(abs(range-max(rangeSel)));
end
function rangeSel=readDefaultAndConvertToIndexRange(rangeSelEdit,range,defaultValues)
    [rangeSelStart rangeSelEnd]=readDefaultAndConvertToIndices(rangeSelEdit,range,defaultValues);
    rangeSel=[rangeSelStart:rangeSelEnd];
end

function RGB=colorMapChannels(imageChannels,colorMaps,colorsForPrint)
    if ~iscell(colorMaps)
        colorMaps={colorMaps};
    end
    if (nargin<3)
        colorsForPrint=false;
    end
    
    inputSize=size(imageChannels);
    if (length(inputSize)<3)
        inputSize(3)=1;
    end
    RGB=zeros([inputSize(1:2) 3]);
    for channelIdx=1:inputSize(3),
        RGB=RGB+mapColor(imageChannels(:,:,channelIdx),colorMaps{channelIdx});
    end
    if (colorsForPrint),
        HSV=rgb2hsv(RGB);
        HSV(:,:,3)=1-HSV(:,:,3); % Swap the intensity
        RGB=hsv2rgb(HSV);
    end
end

function shutDown(obj,event)
    getBaseFigureHandle([]); % Clear the permanent reference
    closereq();
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