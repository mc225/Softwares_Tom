function lightSheet3DShapeSimulation()
    close all;
    
    outputFolderName=[pwd(),'/'];
    
    stepSize=[1 1 1]*.1e-6;
        
    %
    % Create figures
    %
    colorMap=jet(1024);
%     colorMap=interpolatedColorMap(1024,[.5 1 .5; .5 .5 1; 0 0 1; 1 1 0; 1 0 0],[0 0.25 .5 .75 1]);
%     colorMap=hot(1024);
    figures=createFigure(0,1,true,9,true,15,@max,[0.2 0.4],stepSize,.001); colormap(colorMap);
    figures(2)=createFigure(0,0.01,true,9,true,15,@plus,[0.2 0.4],stepSize,.001); colormap(colorMap);
    figures(3)=createFigure(5,1,false,.12,false,0,@max,[0.05 0.2],stepSize,-1); colormap(colorMap);
    
    % Output
    saveWithTransparency(figures(1),strcat(outputFolderName,'GaussianBeamSwipedIntoLightSheet.png'));
    saveWithTransparency(figures(2),strcat(outputFolderName,'BesselBeamSwipedIntoLightSheet.png'));
    saveWithTransparency(figures(3),strcat(outputFolderName,'AiryBeamSwipedIntoLightSheet.png'));
end

function fig=createFigure(alpha,beta,circularSymmetric,isoValues,logScale,colorMapOffset,projector,opaquenessLightSheet,stepSize,noiseLevel)
    % Generates the light sheets
    [psf lightSheet X Y Z]=generateLightSheet(alpha,beta,circularSymmetric,projector,stepSize,noiseLevel);
    
    if (logScale)
        psf=log(max(psf,eps(1)));
        lightSheet=log(max(lightSheet,eps(1)));
    end
    psf=psf+colorMapOffset;
    lightSheet=lightSheet+colorMapOffset;
    
    % Draw the figures
    fig=figure('Position',[50 50 1024 768],'Color','white');
 
    % Draw beam
    if (circularSymmetric)
        beamCenterIdx=1+floor(size(psf,1)/2);
        psfSlice=squeeze(psf(1+floor(end/2):end,beamCenterIdx,:));
        pRange=2*pi*[0:.01:1];
        rRange=Y(1+floor(end/2):end,1,1);
        zRange=squeeze(Z(1,1,:));
        [P,R,Zpolar]=meshgrid(pRange,rRange,zRange);
        [Xpolar,Ypolar,Zpolar]=pol2cart(P,R,Zpolar);
        psfSlice=repmat(permute(psfSlice,[1 3 2]),[1 length(pRange) 1]);
        fv=isosurface(Xpolar,Ypolar,Zpolar,psfSlice,isoValues(1));
        if (~isempty(fv.faces))
            p1=patch(fv,'FaceColor','interp','EdgeColor','none');
            patch(isocaps(Xpolar,Ypolar,Zpolar*1.005,psfSlice,isoValues(1)),'FaceColor','interp','EdgeColor','none');
            faceNormals=cross(fv.vertices(fv.faces(:,1),:)-fv.vertices(fv.faces(:,2),:),fv.vertices(fv.faces(:,1),:)-fv.vertices(fv.faces(:,3),:));
            vertexNormals=zeros(size(fv.vertices));
            for (faceIdx=1:size(faceNormals,1))
                vertexNormals(fv.faces(faceIdx,:),:)=vertexNormals(fv.faces(faceIdx,:),:)+repmat(faceNormals(faceIdx,:),[3 1]);
            end
            vertexNormals=-vertexNormals./repmat(sum(abs(vertexNormals).^2,2),[1 3]);
            set(p1,'VertexNormals',vertexNormals);
        end
        clear psfSlice Xpolar Ypolar Zpolar;
    else
        p1=patch(isosurface(X,Y,Z,psf,isoValues(1)),'FaceColor','white','EdgeColor','none');
        patch(isocaps(X,Y,Z*1.005,psf,isoValues(1)),'FaceColor','interp','EdgeColor','none');
        isonormals(X,Y,Z,psf,p1);
    end
    set(p1,'FaceColor','flat','CData',isoValues(1),'CDataMapping','scaled');
    
    %Draw light sheet
    p2=patch(isosurface(X,Y,Z,lightSheet,isoValues(end)),'FaceColor','white','EdgeColor','none','FaceAlpha',opaquenessLightSheet(1));
    patch(isocaps(X,Y,Z,lightSheet,isoValues(end)),'FaceColor','interp','EdgeColor','none','FaceAlpha',opaquenessLightSheet(2));
    axis equal;
    camtarget([mean(X(:)) mean(Y(:)) mean(Z(:))]);
    camup([0 1 0]);
    campos([40 15 -25]);
    camva(33);
    camproj('orthographic');
    camlight(60,30); camlight; lighting gouraud;
    isonormals(X,Y,Z,lightSheet,p2);
    set(p2,'FaceColor','flat','CData',isoValues(end),'CDataMapping','scaled');
    set(gca,'Visible','Off');
    xlim(X([1 end]));ylim(Y([1 end]));zlim(Z([1 end]));
    drawnow();
end

function [psf lightSheet X Y Z]=generateLightSheet(alpha,beta,circularSymmetric,projector,stepSize,noiseLevel)
    xRange=[-10e-6:stepSize(1):25e-6]; yRange=[-10e-6:stepSize(2):10e-6]; zRange=[-10e-6:stepSize(3):10e-6];
    xRangeShort=xRange(1:length(yRange));
    [X,Y,Z]=ndgrid(xRange,yRange,zRange);
        
    wavelength=532e-9;
    
    pupilFunctor=@(U,V) (sqrt(U.^2+V.^2)<=1 & sqrt(U.^2+V.^2)>(1-beta)).*exp(2i*pi*alpha*(U.^3-V.^3-2/3*(U-V)));
    
    objectiveNumericalAperture=0.42;
    refractiveIndex=1.4;
    magnification=20;
    tubeLength=Inf; % Infinite focal length lens simulated, more typically this should be 200e-3;
    psf=calcVectorialPsf(xRangeShort,yRange,zRange,wavelength,...
                              @(U,V) pupilFunctor(U,V)/sqrt(2),@(U,V) 1i*pupilFunctor(U,V)/sqrt(2), ... % Circ polariz
                              objectiveNumericalAperture,refractiveIndex,magnification,tubeLength);
    
    psf=psf./max(psf(:));
    
    psf(psf<noiseLevel)=0;
    if (circularSymmetric)
        [XShort YShort]=ndgrid(xRangeShort,yRange,zRange);
        R=sqrt(XShort.^2+YShort.^2);
        maxRadius=R(1+floor(end/2),find(psf(1+floor(end/2),:,1)>0,1,'first'));
        psf(R>maxRadius)=0;
        clear XShort YShort;
    end
        
    % Calculate the light-sheet
    beamSize=size(psf);
    swipeLength=beamSize(1);
    lightSheet=psf;
    for swipeIdx=2:swipeLength,
        lightSheet(swipeIdx:end,:,:)=projector(lightSheet(swipeIdx:end,:,:),psf(1:(end-(swipeIdx-1)),:,:));
    end
    extension=length(xRange)-length(xRangeShort);
    lightSheet(end+[1:extension],:,:)=repmat(lightSheet(end,:,:),[extension 1 1]);
    
    lightSheet=lightSheet/max(lightSheet(:));
    
    psf(end+extension,1,1)=0;
    
    psf=permute(psf,[2 1 3]);
    lightSheet=permute(lightSheet,[2 1 3]);
    X=permute(X,[2 1 3])*1e6;
    Y=permute(Y,[2 1 3])*1e6;
    Z=permute(Z,[2 1 3])*1e6;
end