% dataCube=readDataCubeFromAviFile(fileName,projectionDimensionOrSubCube)
%
% Inputs:
%     fileName: string representing the avi file to read as a data cube
%     projectionDimensionOrSubCube: integer (default []). If not empty, the indicated dimension will be returned integrated.
%
% Returns a 3D matrix with the frames of an avi file in single precision
% and normalized to 1;
function dataCube=readDataCubeFromAviFile(fileName,projectionDimensionOrSubCube)
    if (nargin<2)
        projectionDimensionOrSubCube=[];
    end
    
    if (exist('VideoReader','class'))
       vidReader = VideoReader(fileName, 'tag', 'vidReader1');
 
       % Get the data size.
       imgHeight = get(vidReader, 'Height');
       imgWidth = get(vidReader, 'Width');
       nbFrames = get(vidReader, 'numberOfFrames');
       colorEncodingFor16bits=vidReader.BitsPerPixel==24;
        
       % Create a 3D matrix from the video frames.
       dataCube=zeros([imgHeight imgWidth nbFrames],'single');
       for frameIdx = 1 : nbFrames
           vidFrame = read(vidReader,frameIdx);
           if (colorEncodingFor16bits)
               % 16 bit encoded in green and blue channel
               img = single(vidFrame(:,:,3))+single(vidFrame(:,:,2))*256;
           else
               img = mean(single(vidFrame),3);
           end
           if (isempty(projectionDimensionOrSubCube))
               %Complete data cube
               dataCube(:,:,frameIdx)=img;
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
                       dataCube(:,:,frameIdx)=img;
                   else
                       dataCube(:,:,1)=dataCube(:,:,1)+img;
                   end
               else
                   %Crop to a subset of the data cube given by the matrix projectionDimensionOrSubCube.
                   if(frameIdx>=projectionDimensionOrSubCube(3,1) && frameIdx<=projectionDimensionOrSubCube(3,2))
                        img=img(projectionDimensionOrSubCube(1,1):projectionDimensionOrSubCube(1,2),projectionDimensionOrSubCube(2,1):projectionDimensionOrSubCube(2,2));
                        dataCube(:,:,frameIdx-projectionDimensionOrSubCube(3,1)+1)=img;
                   end
               end
           end
       end
    else
        warning off;
            mov=aviread(fileName);
        warning on;
        colorEncodingFor16bits=false;
        
        %Convert to matrix
        dataCubeCellArray=struct2cell(mov);
        dataCubeCellArray=dataCubeCellArray(1,1,:);
        dataCube=cell2mat(dataCubeCellArray);
        clear dataCubeCellArray;
    end
    
    if (colorEncodingFor16bits)
        dataCube=dataCube./(2^16-1);
    else
        dataCube=dataCube./(2^8-1);
    end
end