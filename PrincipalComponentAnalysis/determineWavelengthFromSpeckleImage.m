% [wavelength testImagesInPrincipleComponentSpace] = determineWavelengthFromSpeckleImage(img,calibration)
%
% Determines the wavelength corresponding to a speckle image given the
% calibration data.
%
% Inputs:
%     img: a 2D image, or a 3D stack of speckle testImages
%     calibration: a struct with calibration data
%
% Outputs:
%     wavelength: the wavelength corresponding to the speckle image, or a
%                 vector with values for each image in a stack.
%     testImagesInPrincipleComponentSpace: matrix with in the columns the
%                 coordinates of the testImages in principal component basis
%                
%
function [wavelengths testImagesInPrincipleComponentSpace] = determineWavelengthFromSpeckleImage(testImages,calibration)
    testInputSize=size(testImages);
    if (length(testInputSize)<3)
        testInputSize(3)=1;
    end
    
    nbTrainingSamplesPerWavelength=numel(calibration.trainingWavelengths)/prod(testInputSize(4:end));
    
    % Preprocess
    testImages=testImages./repmat(sqrt(sum(sum(testImages.^2))),[testInputSize(1:2) 1]);
    testImages=testImages-repmat(calibration.meanNormalizedImage,[1 1 testInputSize(3:end)]);
    
    % Vectorize the testImages
    testImages=reshape(testImages,[prod(testInputSize(1:2)) prod(testInputSize(3:end))]);

    % Classify the input
    testImagesInPrincipleComponentSpace=calibration.principalComponentsInImageSpace*testImages;
    %Nearest neighbor
    wavelengths=[];
    for (imgIdx=1:prod(testInputSize(3:end)))
        ind=dsearchn(calibration.trainingImagesInPrincipalComponentSpace',testImagesInPrincipleComponentSpace(:,imgIdx)');
        wavelengths(end+1)=calibration.trainingWavelengths(ind);
        detectedWavelengthIdx=1+floor((ind-1)/nbTrainingSamplesPerWavelength);
    end
%     %Linear
%     wavelengths = classify(testImagesInPrincipleComponentSpace',calibration.trainingImagesInPrincipalComponentSpace',calibration.trainingWavelengths,'linear');
%     %Mahalanobis
%     wavelengths = classify(testImagesInPrincipleComponentSpace',calibration.trainingImagesInPrincipalComponentSpace',calibration.trainingWavelengths,'mahalanobis');
%       %Support Vector Machine
%     wavelengthIndexes = multisvm(calibration.trainingImagesInPrincipalComponentSpace',calibration.trainingWavelengths(:),testImagesInPrincipleComponentSpace');
%     wavelengths=calibration.trainingWavelengths(wavelengthIndexes);
    
    % Shape back to input dimensions
    if (length(testInputSize)>3)
        wavelengths=reshape(wavelengths,testInputSize(3:end));
    end
    nbPrincipalComponents=size(calibration.principalComponentsInImageSpace,1);
    testImagesInPrincipleComponentSpace=reshape(testImagesInPrincipleComponentSpace,[nbPrincipalComponents testInputSize(3:end)]);
end