function createCellRotationMovie(inputFolder,outputFileName)
    close all;
    if (nargin<1 || isempty(inputFolder))
        inputFolder='Z:\RESULTS\2013-05-14_cellspheroidsWGA\2013-05-14_cellspheroidsWGA_filter_ap1-2B';
    end
    if (nargin<2 || isempty(outputFileName))
        outputFileName='rotatingCells.avi';
    end
    
    gaussianFileName=[inputFolder, '/recording_lambda488nm_alpha0_beta100_cropped.mat'];
    airyFileName=[inputFolder, '/recording_lambda488nm_alpha2_beta100_cropped.mat'];

    %Load the data
    [GaussData AiryData xRange yRange zRange]=loadData(gaussianFileName,airyFileName);

    %Normalize the data
    GaussData=.1*GaussData/findPercentileCutOffValue(GaussData(:),.9);
    AiryData=.1*AiryData/findPercentileCutOffValue(AiryData(:),.9);
    
    writerObj = VideoWriter(outputFileName,'Uncompressed AVI');
    set(writerObj,'FrameRate',25);
    open(writerObj);
    
    colorMap=interpolatedColorMap(1024,[0 0 0; 0 .8 0; 1 1 1],[0 .5 1]);
    
    fig=figure('Position',[50 50 640 480],'Color',[0 0 0]);
%     set(fig,'Renderer','zbuffer'); % Turn of Windows 64 bit Aero!!
    [writerObj fig]=constructDataCube(fig,GaussData,AiryData,xRange,yRange,zRange,1,colorMap,writerObj);
    [writerObj fig]=rotateCamera(fig,20*[0 1 0]*pi/180,5,writerObj);
    [writerObj fig]=rotateCamera(fig,2*pi*[1 0 0],25,writerObj,5);
    
    close(writerObj);
end

function cutOffValue=findPercentileCutOffValue(values,percentile)
    values=sort(values);
    cutOffIndex=sum(cumsum(values)<percentile*sum(values));
    cutOffValue=values(max(1,cutOffIndex));
end
    
function [writerObj fig]=constructDataCube(fig,GaussData,AiryData,xRange,yRange,zRange,nbSteps,colorMap,writerObj)
    isovalGauss=0.1;
    isovalAiry=0.1;

    [X,Y,Z]=meshgrid(xRange,yRange,zRange);
    
    % Built up of data cube
    userData={};
    userData.subplots=subplot(1,2,1,'Parent',fig,'Color','none');
    userData.subplots(2)=subplot(1,2,2,'Parent',fig,'Color','none');
    for fractionToShow=[1:nbSteps]./nbSteps,
        userData.subplots(1)=drawPartialVolume(userData.subplots(1),GaussData,isovalGauss,colorMap,'Gaussian',X,Y,Z,fractionToShow);
        userData.subplots(2)=drawPartialVolume(userData.subplots(2),AiryData,isovalAiry,colorMap,'Airy',X,Y,Z,fractionToShow);
        
        drawnow();
        writeVideo(writerObj,getframe(fig));

        if(fractionToShow<1)
            cla(userData.subplots(1));
            cla(userData.subplots(2));
        end
    end

    set(fig,'UserData',userData);
%     linkaxes(userData.subplots);
end

function ax=drawPartialVolume(ax,dataCube,isoval,colorMap,description,X,Y,Z,fractionToShow)
    alphaValue=0.50;
    cameraPosition=[100 -100 500];
    cameraUp=[0 1 0];
    
    fractionSel=size(dataCube,2)-ceil(fractionToShow*size(dataCube,2)-1):size(dataCube,2);
    dataCube=dataCube(:,fractionSel,:);
    
    GaussPatch1=patch(isosurface(X(:,fractionSel,:),Y(:,fractionSel,:),Z(:,fractionSel,:),...
        dataCube,isoval),'Parent',ax,'FaceColor','green',...
        'FaceLighting','phong','AmbientStrength',0.5,...
        'EdgeColor','none','FaceAlpha',alphaValue);
    patch(isocaps(X(:,fractionSel,:),Y(:,fractionSel,:),Z(:,fractionSel,:),...
        dataCube,isoval),'Parent',ax,'FaceColor','interp',...
        'FaceLighting','phong','AmbientStrength',0.5,...
        'EdgeColor','none','FaceAlpha',alphaValue);
    isonormals(X(:,fractionSel,:),Y(:,fractionSel,:),Z(:,fractionSel,:),dataCube,GaussPatch1);
    colormap(ax,colorMap);
    grid(ax,'on'); box(ax,'on');
    axis(ax,'vis3d','image');
    xlim(ax,[X(1) X(end)]);
    ylim(ax,[Y(1) Y(end)]);
    zlim(ax,[Z(1) Z(end)]);
    xlabel(ax,'x [\mu m]','FontSize',14,'FontWeight','bold');
    ylabel(ax,'y [\mu m]','FontSize',14,'FontWeight','bold');
    zlabel(ax,'z [\mu m]','FontSize',14,'FontWeight','bold');
    title(ax,description,'FontSize',20,'FontWeight','bold');
    axisColor=[1 1 1];
    set(ax,'FontSize',16,'LineWidth',1.5,'XColor',axisColor,'YColor',axisColor,'ZColor',axisColor);
%     view(ax,315,30);
    set(ax,'CameraPositionMode','manual','CameraPosition',cameraPosition);
    set(ax,'CameraUpVectorMode','manual','CameraUpVector',cameraUp);
    set(ax,'CameraTargetMode','manual','CameraTarget',[mean(X(:)) mean(Y(:)) mean(Z(:))]);
%     set(ax,'Projection','perspective');
    set(ax,'CameraViewAngleMode','manual','CameraViewAngle',10);
    
    setLights(ax);
end

function [writerObj fig]=rotateCamera(fig,rotationVector,nbSteps,writerObj,nbCycles)
    if (nargin<5 || isempty(nbCycles))
        nbCycles=1;
    end    
    
    userData=get(fig,'UserData');
    initialCameraPosition=get(userData.subplots(1),'CameraPosition');
    initialCameraUpVector=get(userData.subplots(1),'CameraUpVector');
    
    nRAx=rotationVector./norm(rotationVector);
    
    % Rotate to position
    rotationFrames=struct('cdata',{},'colormap',[]);
    for fractionOfTranslation=[1:nbSteps]./nbSteps,
        ct=cos(-fractionOfTranslation*norm(rotationVector));
        st=sin(-fractionOfTranslation*norm(rotationVector));
        rotationMatrix=[ct+nRAx(1)^2*(1-ct)          nRAx(1)*nRAx(2)*(1-ct)-nRAx(3)*st   nRAx(1)*nRAx(3)*(1-ct)+nRAx(2)*st;...
                        nRAx(2)*nRAx(1)*(1-ct)+nRAx(3)*st  ct+nRAx(2)^2*(1-ct)           nRAx(2)*nRAx(3)*(1-ct)-nRAx(1)*st;...
                        nRAx(3)*nRAx(1)*(1-ct)+nRAx(2)*st  nRAx(3)*+nRAx(2)*(1-ct)+nRAx(1)*st  ct+nRAx(3)^2*(1-ct)         ];
        cameraPosition=initialCameraPosition*rotationMatrix;
        cameraUpVector=initialCameraUpVector*rotationMatrix;
        
        set(userData.subplots(1),'CameraPositionMode','manual','CameraPosition',cameraPosition);
        set(userData.subplots(1),'CameraUpVectorMode','manual','CameraUpVector',cameraUpVector);
        setLights(userData.subplots(1));
        set(userData.subplots(2),'CameraPositionMode','manual','CameraPosition',cameraPosition);
        set(userData.subplots(2),'CameraUpVectorMode','manual','CameraUpVector',cameraUpVector);
        setLights(userData.subplots(2));
        
        drawnow();
        rotationFrames(end+1)=getframe(fig);
    end
    
    % Do multiple rotations if requested
    for cycleIdx=1:nbCycles,
        for (rotationFrame=rotationFrames)
            writeVideo(writerObj,rotationFrame);
        end
    end
end

function setLights(ax)
    camH1=camlight(30,30); set(camH1,'Parent',ax);
%     camH2=camlight('left'); set(camH2,'Parent',ax);
    lighting(ax,'phong');
end
    
function [GaussData AiryData xRange yRange zRange]=loadFakeData(gaussianFileName,airyFileName)
    xRange=[-42.8:2:12.8];
    yRange=[-40.4:2:15.4];
    zRange=[-24.4:1:31.4];
    
    rand('seed',0)
    GaussData = rand(length(xRange),length(yRange),length(zRange));
    GaussData = smooth3(GaussData,'box',15);
    AiryData = rand(length(xRange),length(yRange),length(zRange));
    AiryData = smooth3(AiryData,'box',5);
    
    AiryData=swapXY(AiryData);
end

function [GaussData AiryData xRange yRange zRange]=loadData(gaussianFileName,airyFileName)
    load(gaussianFileName,'recordedImageStack','xRange','yRange','zRange');
    GaussData=recordedImageStack;
    load(airyFileName,'restoredDataCube');
    AiryData=restoredDataCube;
    xRange=xRange*1e6;
    yRange=yRange*1e6;
    zRange=zRange*1e6;
    
    [GaussData,xRange,yRange]=swapXY(GaussData,xRange,yRange);
    AiryData=swapXY(AiryData);
    
%     % Subsample
%     step=4;
%     xRange=xRange(1:step:end);
%     yRange=yRange(1:step:end);
%     zRange=zRange(1:step:end);
%     GaussData=GaussData(1:step:end,1:step:end,1:step:end);
%     AiryData=AiryData(1:step:end,1:step:end,1:step:end);
end
   
function [dataCube xRange yRange]=swapXY(dataCube,xRange,yRange)
%     if(nargin>=3)
%         tmp=xRange;
%         xRange=yRange;
%         yRange=tmp;
%     end
    dataCube=permute(dataCube,[2 1 3]);
end