% [dataCube maxValue]=readDataCubeFromAviFile(fileName,projectionDimensionOrSubCube,frameIndexes,normalize)
%
% Inputs:
%     fileName: string representing the avi file to read as a data cube, or a VideoReader object.
%     projectionDimensionOrSubCube: integer (default []). If not empty, the indicated dimension will be returned integrated.
%     frameIndexes: optional list of indexes to read, use -1 to select all (default: -1: all)
%     normalize: boolean indicating if the output data should be normalized to the dynamic range of the input, default: true
%
% Outputs:
%     dataCube: Returns a 3D matrix with the frames of an avi file in single precision and normalized to 1 unless the normalize argument is false;
%     maxValue: the maximum value that could have been stored in this file.
%
function [dataCube maxValue]=readDataCubeFromAviFile(fileName,projectionDimensionOrSubCube,frameIndexes,normalize)
    if (nargin<2)
        projectionDimensionOrSubCube=[];
    end
    if (nargin<3)
        frameIndexes=-1;
    end
    frameIndexes=sort(frameIndexes);
    if (nargin<4)
        normalize=true;
    end
    
	if (exist('VideoReader','class'))
        if (ischar(fileName))
            vidReader = VideoReader(fileName, 'tag', 'vidReader1');
        else
            vidReader=fileName;
        end

        % Get the data size.
        imgSize = [vidReader.Height, vidReader.Width];
        nbFrames = vidReader.NumberOfFrames;
    
        if (any(frameIndexes<0))
            frameIndexes=[1:nbFrames];
        end
        frameIndexes=intersect([1:nbFrames],frameIndexes);
        colorEncodingFor16bits=vidReader.BitsPerPixel==24;
        
        % Create a 3D matrix from the video frames.
        dataCube=zeros([imgSize length(frameIndexes)],'single');
        for frameIndexIdx = 1:length(frameIndexes)
            frameIndex=frameIndexes(frameIndexIdx);
            vidFrame = read(vidReader,frameIndex);
            if (colorEncodingFor16bits)
                % 16 bit encoded in green and blue channel
                img = single(vidFrame(:,:,3))+single(vidFrame(:,:,2))*256;
            else
                img = mean(single(vidFrame),3);
            end
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
        if (ischar(fileName))
            delete(vidReader);
        end
    else
         warning off;
             mov=aviread(fileName,frameIdx);
         warning on;
         colorEncodingFor16bits=false;
         
         %Convert to matrix
         dataCubeCellArray=struct2cell(mov);
         dataCubeCellArray=dataCubeCellArray(1,1,:);
         dataCube=cell2mat(dataCubeCellArray);
         clear dataCubeCellArray;
	end
    
    if (colorEncodingFor16bits)
        maxValue=2^16-1;
    else
        maxValue=2^8-1;
    end
    if (normalize)
        dataCube=dataCube./maxValue;
    end
end