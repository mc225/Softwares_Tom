function createResolutionPlot(inputFolder)
%     close all;
    
%     beamTypes={'Gaussian','Airy'};
    
    if nargin<1 || isempty(inputFolder)
%         inputFolder='C:\Users\Tom\Dropbox\AirySheetData\FWHMs\forwardscans2';
%         inputFolder='C:\Users\Tom\Dropbox\AirySheetData\FWHMs\forwardscans';
%         inputFolder='C:\Users\Tom\Dropbox\AirySheetData\FWHMs\backwardscans';
%         inputFolder='C:\Users\Tom\Dropbox\AirySheetData\FWHMs\morePowerBackward';
        beamTypes={'Bessel10','Bessel5'};
        inputFolder='C:\Users\Tom\Dropbox\AirySheetData\FWHMs\morePowerForward';
    end
    
    description='';
    for beamTypeIdx=1:length(beamTypes)
        description=[description,'_',beamTypes{beamTypeIdx}];
    end
    
    outputFilePath=['C:\Users\Tom\Documents\nanoscope\docs\LightSheet\paper\figures\resolutionVsFOV\fwhmZOfBeads',description,'.svg'];
    
    samples{1}=load(fullfile(inputFolder,['images',beamTypes{1},'.mat']));
    samples{1}=samples{1}.samples;
    samples{2}=load(fullfile(inputFolder,['images',beamTypes{2},'.mat']));
    samples{2}=samples{2}.samples;
    
    %Keep only those points that fall within the plot
    positions=cell2mat({samples{1}.centroid}.');
    inPlot=positions(:,1)>=-20e-6;
    samples{1}=samples{1}(inPlot);
    samples{2}=samples{2}(inPlot);
    positions=positions(inPlot,:);
    
    positions=cell2mat({samples{1}.centroid}.');
    fwhm{1}=cell2mat({samples{1}.fwhm}.');
    fwhm{2}=cell2mat({samples{2}.fwhm}.');
%     distanceFromWaist=-positions(:,2);
    distanceFromWaist=abs(positions(:,2));
    
    fig=figure('Position',[50 50 800 300],'NumberTitle','off','Name',description);
    ax=axes();
    scatterGaussian=scatterV6(distanceFromWaist*1e6,fwhm{1}(:,3)*1e6,'x','MarkerEdgeColor',[1 0 0],'LineWidth',2);
    hold on;
    scatterAiry=scatterV6(distanceFromWaist*1e6,fwhm{2}(:,3)*1e6,50,[0 0.5 0],'filled');
    hold off;
    xlim([0 150]); ylim([0 25]);
    labs(1)=xlabel('x [\mum]'); labs(2)=ylabel('FWHM [\mum]');
    set(ax,'XTick',[0:25:1000]); set(ax,'YTick',[0:5:100]);
    set(ax,'LineWidth',2);
    set([ax labs],'FontSize',22,'FontWeight','bold','FontName','Arial');
    set(labs,'FontSize',24);
    set(ax,'Box','on');
    drawnow();
    
    %attach Callback
%     set([scatterGaussian scatterAiry],'HitTest','off');
    sampleMarkersGaussian=scatterGaussian; %get(scatterGaussian,'Children');
    sampleMarkersAiry=scatterAiry; %get(scatterAiry,'Children');
    set([sampleMarkersGaussian sampleMarkersAiry],'HitTest','on');
    for sampleIdx=1:length(sampleMarkersGaussian)
        bothMarkers=[sampleMarkersGaussian(sampleIdx) sampleMarkersAiry(sampleIdx)];
        set(sampleMarkersGaussian(sampleIdx),'ButtonDownFcn',@(obj,evt) inspectSample(bothMarkers,'gaussian',sampleIdx,samples{1}(sampleIdx),samples{2}(sampleIdx)));
        set(sampleMarkersAiry(sampleIdx),'ButtonDownFcn',@(obj,evt) inspectSample(bothMarkers,'airy',sampleIdx,samples{1}(sampleIdx),samples{2}(sampleIdx)));
    end
    
    uicontrol(fig,'Style','pushbutton','String','Save','Callback',@saveFigure,'UserData',outputFilePath);
    
    set(fig,'ToolBar','figure');
    zoom(fig);
    
end

function markers=scatterV6(X,Y,S,varargin)
    nbMarkers=length(X);
    maxMarkers=50;
    markers=zeros(size(X));
    for markerIdx=1:maxMarkers:nbMarkers
        opts=varargin;
        selI=markerIdx:min(nbMarkers,markerIdx+maxMarkers-1);
        if all(ischar(S)) || isempty(opts) || all(isscalar(opts{1}))
            scatterGroup=scatter(X(selI),Y(selI),S);
        else
            if (length(opts)>1 && strcmpi(opts{2},'filled'))
                scatterGroup=scatter(X(selI),Y(selI),S,opts{1:2});
                opts={opts{3:end}};
            else
                scatterGroup=scatter(X(selI),Y(selI),S,opts{1});
                opts={opts{2:end}};
            end
        end
        for (argIdx=1:2:length(opts)-1)
            set(scatterGroup,opts{argIdx},opts{argIdx+1});
        end
        hold on;
        newMarkers=get(scatterGroup,'Children');
        markers(markerIdx-1+[1:length(newMarkers)])=newMarkers(end:-1:1); % Assume that the order is swapped!
    end
    hold off;
end

function saveFigure(obj,evt)
    fig=get(obj,'Parent');
    outputFilePath=get(obj,'UserData');
    set(obj,'Visible','off');
    logMessage('Writing figure to %s.',outputFilePath);
    plot2svg(outputFilePath,fig);
    set(obj,'Visible','on');
end

function inspectSample(sampleMarkers,sampleSet,sampleIdx,sampleG,sampleA)
    description=sprintf([sampleSet ' sample %d, centroid=(%0.1f,%0.1f,%0.1f)um, FWHM=%0.0f and %0.0f nm'],[sampleIdx,sampleA.centroid*1e6 [sampleG.fwhm(3) sampleA.fwhm(3)]*1e9]);
    logMessage(description);
    
    colorMap=hot(1024);
    
    inspectFig=figure('Position',[200 200 640 480]);
    subplot(2,3,1);showImage(mapColor(sampleG.test.proj3./max(sampleG.test.proj3(:)),colorMap),[],sampleG.yRange*1e6,sampleG.xRange*1e6);axis equal;
    xlabel('x (propagation) [um]'); ylabel('y (swipe) [um]');
    subplot(2,3,2);showImage(mapColor(sampleG.test.proj2./max(sampleG.test.proj2(:)),colorMap),[],sampleG.yRange*1e6,sampleG.zRange*1e6); axis equal;
    xlabel('x (propagation) [um]'); ylabel('z (scan) [um]');
    subplot(2,3,3);showImage(mapColor(sampleG.test.proj1./max(sampleG.test.proj1(:)),colorMap),[],sampleG.xRange*1e6,sampleG.zRange*1e6); axis equal;
    xlabel('y (swipe) [um]'); ylabel('z (scan) [um]');
    subplot(2,3,4);showImage(mapColor(sampleA.test.proj3./max(sampleA.test.proj3(:)),colorMap),[],sampleA.yRange*1e6,sampleA.xRange*1e6);axis equal;
    xlabel('x (propagation) [um]'); ylabel('y (swipe) [um]');
    subplot(2,3,5);showImage(mapColor(sampleA.test.proj2./max(sampleA.test.proj2(:)),colorMap),[],sampleA.yRange*1e6,sampleA.zRange*1e6); axis equal;
    xlabel('x (propagation) [um]'); ylabel('z (scan) [um]');
    subplot(2,3,6);showImage(mapColor(sampleA.test.proj1./max(sampleA.test.proj1(:)),colorMap),[],sampleA.xRange*1e6,sampleA.zRange*1e6); axis equal;
    xlabel('y (swipe) [um]'); ylabel('z (scan) [um]');
    
    set(inspectFig,'NumberTitle','off','Name',description);
    set(inspectFig,'ToolBar','figure');
    zoom(inspectFig);
    
    drawnow();

    keepAction=@(obj,evt) close(inspectFig);
    
    uicontrol(inspectFig,'Style','pushbutton','String','Keep','Position',[10 10 50 20],'Callback',keepAction);
    uicontrol(inspectFig,'Style','pushbutton','String','Remove','Position',[80 10 50 20],'Callback',@(obj,evt) removeMarkers(sampleMarkers,inspectFig));
    
end
function removeMarkers(sampleMarkers,inspectFig)
    set(sampleMarkers,'Visible','off');
    close(inspectFig);
end
    