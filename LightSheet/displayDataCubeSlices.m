function displayDataCubeSlices(folderName,alphasToShow)
%     close all;
    if (nargin<1 || isempty(folderName))
        rootFolder='F:\RESULTS\';
%         folderName=[rootFolder,'03_27_2013\Acini_1\2013-03-27 10_34_27.693'];
%         folderName=[rootFolder,'03_27_2013\Acini_1\2013-03-27 10_38_45.594']; % quite good
%         folderName=[rootFolder,'03_27_2013\Acini_2\2013-03-27 10_43_10.564']; % quite good
%         folderName=[rootFolder,'03_27_2013\Acini_2\2013-03-27 10_47_08.663']; % A bit blurry
%         folderName=[rootFolder,'03_27_2013\Acini_3\2013-03-27 11_04_52.996']; % Bleached beam
%         folderName=[rootFolder,'03_27_2013\Acini_4\2013-03-27 11_10_20.581']; % out of FOV
%         folderName=[rootFolder,'03_27_2013\Acini_4\2013-03-27 11_14_06.638']; % no good
%         folderName=[rootFolder,'03_27_2013\Bead_Test\2013-03-27 10_01_06.201_FullNA']; % in NA=0.8, out NA about 0.4
%         folderName=[rootFolder,'03_27_2013\Bead_Test\2013-03-27 10_06_29.359']; % quite nice 488nm (both beads visible)
%         folderName=[rootFolder,'03_27_2013\Bead_Test\2013-03-27 10_14_31.638']; % good, 532nm (only large beads)

%          folderName='F:\RESULTS\redbeads2013-04-19halfAp\2013-04-19 10_42_17.196';
%          folderName=[rootFolder,'redbeads2013-04-18\2013-04-18 17_11_51.448'];
%          folderName=[rootFolder,'redbeads2013-04-19halfAp\2013-04-19 10_42_17.196'];
%         folderName=[rootFolder,'2013-04-22beads\2013-04-22 18_38_59.782\'];
%         folderName=[rootFolder,'2013-04-22cellsHEK\2013-04-22 19_14_42.392_SecondCellNearWaist\'];
%         folderName=[rootFolder,'2013-04-22cellsHEK\2013-04-22 19_15_51.708_SecondCelMoved\'];
%         folderName=['F:\RESULTS\2013-04-23beadsApOneThird\2013-04-23 16_46_12.337\'];
%         folderName='F:\RESULTS\2013-04-23HEK_ApOneThird';
%          folderName='F:\RESULTS\redbeads2013-04-18\2013-04-18 17_11_51.448';
%          folderName='F:\RESULTS\redbeads2013-04-19halfAp\2013-04-19 10_42_17.196';
         folderName='Z:\RESULTS\2013-04-24HEK_ApOneThird\2013-04-24 12_52_17.623';
    end
    if (nargin<2 || isempty(alphasToShow))
        alphasToShow=[];
    end
    
    dataCubeNames=dir([folderName,'/*.mat']);
    
    if (~isempty({dataCubeNames.name}))
        figs=[];
        for (fileName={dataCubeNames.name})
            fileName=fileName{1};
            fullFileName=[folderName,'/',fileName];

            % Determine the light sheet type
            alpha=regexpi(fullFileName,'alpha([\d]*)_[^\\/]+','tokens');
            alpha=str2double(alpha{1});
            if (isempty(alphasToShow) || any(alpha==alphasToShow))
                if (alpha>0)
                    lightSheetType=sprintf('Airy%0.0f',alpha);
                else
                    lightSheetType='Gaussian';
                end
                figs(end+1)=updateFigure([],[],[],fullFileName,lightSheetType);
            else
                logMessage('Skipping the recording with alpha=%0.1f.',alpha);
            end
        end
    else
        logMessage('No files found in folder %s.',folderName);
    end
end

function fig=updateFigure(hObj,event,fig,fullFileName,lightSheetType)
    if (nargin<4 || isempty(fullFileName))
        userData=get(fig,'UserData');
        fullFileName=userData.fullFileName;
    end
    if (nargin<5 || isempty(lightSheetType))
        lightSheetType=userData.lightSheetType;
    end
    
%     projector=@mean;
    projector=@(X,dim) max(X,[],dim);
    

    if (isempty(fig));
        % Load the data first
        wavelength=regexpi(fullFileName,'lambda([\d]*)nm[^\\/]+','tokens');
        wavelength=str2double(wavelength{1})*1e-9;
        logMessage(['Loading ',fullFileName,'...']);
        loadedData=load(fullFileName,'xRange','yRange','zRange','tRange','setupConfig');

        % Open the figure window
        fig=figure('Name',[lightSheetType,sprintf(' at %0.0fnm',[wavelength*1e9]),' - ',fullFileName]);

        userData={};
        userData.lightSheetType=lightSheetType;
        userData.fullFileName=fullFileName;
        userData.xRange=loadedData.xRange;
        userData.yRange=loadedData.yRange;
        userData.zRange=loadedData.zRange;
        userData.tRange=loadedData.tRange;
        loadedData=load(fullFileName,'recordedImageStack');
        userData.recordedImageStack=loadedData.recordedImageStack;
        loadedData=load(fullFileName,'restoredDataCube');
        userData.restoredDataCube=loadedData.restoredDataCube;
        clear loadedData;
        userData.colorMap=hot(1024);

        % Default to complete data cube
        xRangeSel=[1:length(userData.xRange)];
        yRangeSel=[1:length(userData.yRange)];
        zRangeSel=[1:length(userData.zRange)];
    else
        xRangeSel=str2num(get(userData.xRangeSelEdit,'String'))*1e-6;
        if(length(xRangeSel)==0)
            xRangeSel=[1:length(userData.xRange)];
        else
            % Convert metric to indexes
            [ign xRangeSelStart]=min(abs(userData.xRange-xRangeSel(1)));
            [ign xRangeSelEnd]=min(abs(userData.xRange-xRangeSel(end)));
            xRangeSel=[xRangeSelStart:xRangeSelEnd];
        end
        yRangeSel=str2num(get(userData.yRangeSelEdit,'String'))*1e-6;
        if(length(yRangeSel)==0)
            yRangeSel=[1:length(userData.yRange)];
        else
            % Convert metric to indexes
            [ign yRangeSelStart]=min(abs(userData.yRange-yRangeSel(1)));
            [ign yRangeSelEnd]=min(abs(userData.yRange-yRangeSel(end)));
            yRangeSel=[yRangeSelStart:yRangeSelEnd];
        end
        zRangeSel=str2num(get(userData.zRangeSelEdit,'String'))*1e-6;
        if(length(zRangeSel)==0)
            zRangeSel=[1:length(userData.zRange)];
        else
            % Convert metric to indexes
            [ign zRangeSelStart]=min(abs(userData.zRange-zRangeSel(1)));
            [ign zRangeSelEnd]=min(abs(userData.zRange-zRangeSel(end)));
            zRangeSel=[zRangeSelStart:zRangeSelEnd];
        end
    end
    tRangeSel=zRangeSel;

    % Draw the controls
    uicontrol('Style','text','String','X-stack position','Position',[10,50,100,20]);
    userData.yRangeSelEdit=uicontrol('Style','edit','String',sprintf('%0.0f %0.0f',userData.yRange(yRangeSel([1 end]))*1e6),'Position',...
        [130,50,60,20],'Callback',{@updateFigure,fig});
    uicontrol('Style','text','String','Y-stack position','Position',[10,30,100,20]);
    userData.xRangeSelEdit=uicontrol('Style','edit','String',sprintf('%0.0f %0.0f',userData.xRange(xRangeSel([1 end]))*1e6),'Position',...
        [130,30,60,20],'Callback',{@updateFigure,fig});
    uicontrol('Style','text','String','Z-stack position','Position',[10,10,100,20]);
    userData.zRangeSelEdit=uicontrol('Style','edit','String',sprintf('%0.0f %0.0f',userData.zRange(zRangeSel([1 end]))*1e6),'Position',...
        [130,10,60,20],'Callback',{@updateFigure,fig});
    set(gcf,'Toolbar','figure');
    drawnow();

    % Calculate the projections
    recordedXZ=squeeze(projector(userData.recordedImageStack(xRangeSel,yRangeSel,zRangeSel),1)).';
    restoredXZ=squeeze(projector(userData.restoredDataCube(xRangeSel,yRangeSel,tRangeSel),1)).';
    recordedXY=squeeze(projector(userData.recordedImageStack(xRangeSel,yRangeSel,zRangeSel),3));
    restoredXY=squeeze(projector(userData.restoredDataCube(xRangeSel,yRangeSel,zRangeSel),3));
    normalizationRaw=max([max(recordedXZ(:)), max(recordedXY(:))]);
    normalizationDeconvolved=max([max(restoredXZ(:)), max(restoredXY(:))]);
    recordedXZ=recordedXZ./normalizationRaw;
    restoredXZ=restoredXZ./normalizationDeconvolved;
    recordedXY=recordedXY./normalizationRaw;
    restoredXY=restoredXY./normalizationDeconvolved;

    % Draw the projections
    axs=subplot(2,2,1);
    showImage(mapColor(recordedXZ,userData.colorMap),[],userData.yRange(yRangeSel)*1e6,userData.zRange(zRangeSel)*1e6);
    axis equal;
    set(gca,'LineWidth',3,'FontSize',16,'FontWeight','bold');
    title(['recorded ',userData.lightSheetType],'FontSize',22,'FontWeight','bold');
    xlabel('X (propagation) [\mum]','LineWidth',3,'FontSize',22,'FontWeight','bold');
    ylabel('Z (axial) [\mum]','LineWidth',3,'FontSize',22,'FontWeight','bold');
    axs(2)=subplot(2,2,2);
    showImage(mapColor(restoredXZ,userData.colorMap),[],userData.yRange(yRangeSel)*1e6,userData.tRange(tRangeSel)*1e6);
    axis equal;
    set(gca,'LineWidth',3,'FontSize',16,'FontWeight','bold');
    linkaxes(axs(1:2));
    title(['deconvolved ',userData.lightSheetType],'FontSize',22,'FontWeight','bold');
    xlabel('X (propagation) [\mum]','LineWidth',3,'FontSize',22,'FontWeight','bold');
    ylabel('Z (axial) [\mum]','LineWidth',3,'FontSize',22,'FontWeight','bold');

    axs(3)=subplot(2,2,3);
    showImage(mapColor(recordedXY,userData.colorMap),[],userData.yRange(yRangeSel)*1e6,userData.xRange(xRangeSel)*1e6);
    axis equal;
    set(gca,'LineWidth',3,'FontSize',16,'FontWeight','bold');
    title(['recorded ',lightSheetType],'FontSize',22,'FontWeight','bold');
    xlabel('X (propagation) [\mum]','LineWidth',3,'FontSize',22,'FontWeight','bold');
    ylabel('Y (beam spread) [\mum]','LineWidth',3,'FontSize',22,'FontWeight','bold');
    axs(4)=subplot(2,2,4);
    showImage(mapColor(restoredXY,userData.colorMap),[],userData.yRange(yRangeSel)*1e6,userData.xRange(xRangeSel)*1e6);
    axis equal;
    set(gca,'LineWidth',3,'FontSize',16,'FontWeight','bold');
    linkaxes(axs(3:4));
    title(['deconvolved ',userData.lightSheetType],'FontSize',22,'FontWeight','bold');
    xlabel('X (propagation) [\mum]','LineWidth',3,'FontSize',22,'FontWeight','bold');
    ylabel('Y (beam spread) [\mum]','LineWidth',3,'FontSize',22,'FontWeight','bold');

    set(fig,'UserData',userData);
end