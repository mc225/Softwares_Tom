% Takes a subvolume out of corresponding red and green datacubes, makes a
%  maximum intensity projection of the 2 subvolumes, applies the
%  corresponding colormaps to the data and combines in a 3-colour image.
%  Can also output the normalisation factors used to allow the same
%  normalisation to be used across multiple images.

function [TwoColorImage,normalisation]=extractSubVolumeThenProject_2Colors(dataType,alpha,beta,subVolume_x,subVolume_y,subVolume_z,transpose,normalisation)

    fileName_g=strcat('recording_lambda488nm_alpha',num2str(alpha),'_beta',num2str(beta));
    fileName_r=strcat('recording_lambda532nm_alpha',num2str(alpha),'_beta',num2str(beta));
    matfile_g=matfile(fileName_g);
    matfile_r=matfile(fileName_r);
    
    %guess the projection direction - minimum sized dimension
    subVolumeSize=[length(subVolume_x),length(subVolume_y),length(subVolume_z)];
    if min(subVolumeSize)==length(subVolume_x)
        ProjectionGuess=1;
    elseif min(subVolumeSize)==length(subVolume_y)
        ProjectionGuess=2;
    elseif min(subVolumeSize)==length(subVolume_z)
        ProjectionGuess=3;
    end
    
    if strcmp(dataType,'deconv')==1
    image_g=squeeze(max(matfile_g.restoredDataCube(subVolume_x,subVolume_y,subVolume_z),[],ProjectionGuess));
    image_r=squeeze(max(matfile_r.restoredDataCube(subVolume_x,subVolume_y,subVolume_z),[],ProjectionGuess));
    elseif strcmp(dataType,'record')==1
    image_g=squeeze(max(matfile_g.recordedImageStack(subVolume_x,subVolume_y,subVolume_z),[],ProjectionGuess));
    image_r=squeeze(max(matfile_r.recordedImageStack(subVolume_x,subVolume_y,subVolume_z),[],ProjectionGuess));
    end
    
    if strcmp(transpose,'Y_transp')==1
        image_g=image_g.';
        image_r=image_r.';
    end
    
    image_g=image_g.*(image_g>0);
    image_r=image_r.*(image_r>0);
    
    if size(normalisation~=[1 2]);
        normalisation=zeros(1,2);
        normalisation(1)=max(image_g(:));
        normalisation(2)=max(image_r(:));
    end
    
    image_g=image_g./normalisation(1);
    image_r=image_r./normalisation(2);
    
    colorMap=interpolatedColorMap(1024,[0 0 0; 0 .8 0; 1 1 1],[0 .5 1]);
    image_g_RGB=mapColor(image_g,colorMap);
    colorMap=interpolatedColorMap(1024,[0 0 0; .8 0 0; 1 1 1],[0 .5 1]);
    image_r_RGB=mapColor(image_r,colorMap);
    image_RGB=image_r_RGB+image_g_RGB;
    TwoColorImage=min(1,image_RGB);

end
