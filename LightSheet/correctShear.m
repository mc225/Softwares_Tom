% correctedDataCube=correctShear(dataCube,shearFactor)
%
% dataCube: a three-dimensional matrix
% shearFactor: the shear to be undone in units of pixels for the first and second dimension, respectively.
function dataCube=correctShear(dataCube,shearFactor)
    if (~isempty(shearFactor) && ~all(shearFactor==0))
        logMessage('Correcting a shear of (%0.3f,%0.3f) pixels...',shearFactor);
        xRange=[1:size(dataCube,1)];
        yRange=[1:size(dataCube,2)];
        for zIdx=1:size(dataCube,3)
            shearing=shearFactor*(zIdx-floor(1+size(dataCube,3)/2));
            dataCube(:,:,zIdx)=interp2(dataCube(:,:,zIdx),yRange+shearing(2),xRange.'+shearing(1),'*cubic',0);
        end
    end
end