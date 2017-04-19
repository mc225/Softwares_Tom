% [dataCube maxValue]=readDataCubeFromTiffFile(fileName,projectionDimensionOrSubCube,frameIndexes,normalize)
%
% Inputs:
%     fileName: string representing the tiff file to read as a data cube
%     projectionDimensionOrSubCube: integer (default []). If not empty, the indicated dimension will be returned integrated.
%     frameIndexes: optional list of indexes to read, use -1 to select all (default: -1: all)
%     normalize: boolean indicating if the output data should be normalized to the dynamic range of the input, default: true
%
% Outputs:
%     dataCube: the 3D matrix of values
%     maxValue: the maximum value that could have been stored in this file.
%
% Returns a 3D matrix with the frames of a tiff file in single precision and normalized to 1;
%
function [dataCube maxValue]=readDataCubeFromTiffFile(fileName,projectionDimensionOrSubCube,frameIndexes,normalize)
    if (nargin<2)
        projectionDimensionOrSubCube=[];
    end
    if (nargin<3)
        frameIndexes=-1;
    end
    frameIndexes=sort(frameIndexes);
    if (nargin<4 || isempty(normalize))
        normalize=true;
    end
    
    info=imfinfo(fileName);
    imgSize=[info(1).Height info(1).Width];
    nbFrames=length(info);
    maxValue=2^(info(1).BitDepth)-1;
    
    if (length(frameIndexes)==1 && frameIndexes(1)==-1)
        frameIndexes=[1:nbFrames];
    else
        frameIndexes=intersect([1:nbFrames],frameIndexes);
    end
    
    dataCube=zeros([imgSize length(frameIndexes)],'single');
    for frameIndexIdx = 1:length(frameIndexes)
        frameIndex=frameIndexes(frameIndexIdx);
        img=single(imread(fileName,'Index',frameIndex));
        % (project and) store
        if (isempty(projectionDimensionOrSubCube))
            %Complete data cube
            dataCube(:,:,frameIndexIdx)=img;
        else
            %Project or crop the data cube
            if (max(size(projectionDimensionOrSubCube))==1)
                %Project the full cube along one dimension specified by projectionDimensionOrSubCube
                if (any(projectionDimensionOrSubCube==1))
                    img=max(img,[],1);
                end
                if (any(projectionDimensionOrSubCube==2))
                    img=max(img,[],2);
                end
                if (~any(projectionDimensionOrSubCube==3))
                    dataCube(:,:,frameIndexIdx)=img;
                else
                    dataCube(:,:,1)=dataCube(:,:,1)+img;
                end
            else
                %Crop to a subset of the data cube given by the matrix projectionDimensionOrSubCube.
                if(frameIdx>=projectionDimensionOrSubCube(3,1) && frameIdx<=projectionDimensionOrSubCube(3,2))
                     img=img(projectionDimensionOrSubCube(1,1):projectionDimensionOrSubCube(1,2),projectionDimensionOrSubCube(2,1):projectionDimensionOrSubCube(2,2));
                     dataCube(:,:,frameIndexIdx-projectionDimensionOrSubCube(3,1)+1)=img;
                end
            end
        end
    end
    if (normalize)
        dataCube=dataCube./maxValue;
    end
end

