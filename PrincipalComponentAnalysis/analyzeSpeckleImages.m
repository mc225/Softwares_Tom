% analyzeSpeckleImages(inputFolderName,outputFileName)
%
% inputFolderName: a string representing the input folder created by analyzeSpeckleImages.m
% outputFileName: the filename of the output .mat file.#
%                 Default: the same as the input folder.
%
function analyzeSpeckleImages(inputFolderName,outputFileName)
    close all;

    if (nargin<1 || isempty(inputFolderName))
 %       inputFolderName='measurement_2013-01-16_12_10_55_alumina1_50pm';
%          inputFolderName='measurement2013-01-15_13_20_53_alumina1_10pm';
%        inputFolderName='measurement_2013-01-16_12_40_07_beamprofile_50pm';
%        inputFolderName='D:\alumina2\measurement_2013-01-17_17_21_20_longAndRedone';
%        inputFolderName='D:\alumina2\measurement_2013-01-17_21_02_12_785.3032nm_repeated';
%         inputFolderName='D:\SpeckleWaveMeter\NarrowBand\FabryPerot\withoutDiffuserAnd2F\measurement_2013-01-27_12_53_23_2Fsetup';
%         inputFolderName='D:\SpeckleWaveMeter\BroadBand\Alumina\measurement_2013-01-27_21_28_41_gain0_100pm_695nm_705nm';
        inputFolderName='D:\SpeckleWaveMeter\BroadBand\FabryPerot\measurement_2013-01-27_18_47_13_gain0_100pm_695nm_700nm_and_1000pm_700nm_720nm';
    end
    if (nargin<2 || isempty(outputFileName))
        outputFileName=strcat(inputFolderName,'.mat');
    end
    
    %
    % Read measurement data from disk
    %
    % Detect which measurements have been done
    files=dir(strcat(inputFolderName,'/imagesForWavelength*nm.mat'));
    fileNames={files.name}; clear files;
    wavelengths=regexp(fileNames,'imagesForWavelength([\d]+(\.[\d]*)?)nm.mat','tokens','once');
    wavelengths=cellfun(@(c) str2double(c)*1e-9,wavelengths);
    
    % Use all wavelengths
    wavelengthSelection=logical(ones(size(wavelengths)));
    % Handle only selected wavelengths
    %wavelengthSelection=wavelengths<=785.170e-9 & wavelengths>=785.160e-9;
%     wavelengthSelection=mod(wavelengths,5e-12)==0 | abs(wavelengths-785.174e-9)<1e-15;
    
    % Some configuration
    nbOfTrainingImages=50;
    nbOfTestImages=10;
    calibration={};
    calibration.regionOfInterest=[100 0 300 300];
    maxNumberOfPrincipalComponents=9;
    signalToNoise=10.0;
    
    % Read all files and pick out a smaller section
    allImages=single([]);
    testWavelengths=[];
    for (fileIdx=find(wavelengthSelection))
        currentFileName=strcat(inputFolderName,'/',fileNames{fileIdx})
        load(currentFileName,'images');
        images=images(calibration.regionOfInterest(1)+[1:calibration.regionOfInterest(3)],calibration.regionOfInterest(2)+[1:calibration.regionOfInterest(4)],1:(nbOfTrainingImages+nbOfTestImages));
        newTrainingWavelengths=ones(1,nbOfTrainingImages)*wavelengths(fileIdx);
        newTestWavelengths=ones(1,nbOfTestImages)*wavelengths(fileIdx);
%         images=mean(images,3); % Average to safe memory
        if (~isempty(allImages))
            allImages(:,:,:,end+1)=images;
            calibration.trainingWavelengths(end+[1:nbOfTrainingImages])=newTrainingWavelengths;
            testWavelengths(end+[1:nbOfTestImages])=newTestWavelengths;
        else
            allImages=images;
            calibration.trainingWavelengths=newTrainingWavelengths;
            testWavelengths=newTestWavelengths;
        end
    end
    wavelengths=wavelengths(wavelengthSelection);
%     % Show the intensity spot for the 2F system FP
%     figure;
%     plot(wavelengths*1e9,1000*squeeze(mean(mean(allImages(125:165,200:232,1,:)))).'); xlim(wavelengths([1 end])*1e9);
    
    %
    % Pre-process the data and split it in a training and a test data set
    %
    %Split data set in training and test set
    trainingImages=allImages(:,:,1:nbOfTrainingImages,:);
    testImages=allImages(:,:,nbOfTrainingImages+[1:nbOfTestImages],:);
    clear allImages;
    testInputSize=size(testImages);
    trainingInputSize=size(trainingImages);
    %Normalize the training images to L2 norm 1 and subtract the mean
    trainingImages=trainingImages./repmat(sqrt(sum(sum(trainingImages.^2))),[trainingInputSize(1:2) 1 1]);
    calibration.meanNormalizedImage=mean(trainingImages(:,:,:),3);
    trainingImages=trainingImages-repmat(calibration.meanNormalizedImage,[1 1 trainingInputSize(3:end)]);
    
    logMessage('Using %d images for training and %d for testing.',[trainingInputSize(3) testInputSize(3)]);

    % Vectorize the images, so each images is a column in a 2D matrix
    vectorizedImages=reshape(trainingImages,[prod(trainingInputSize(1:2)) prod(trainingInputSize(3:end))]).';
    clear trainingImages;
    
    %
    % Calculate the principal components
    %
    [eigenVectors,eigenValues]=eig(vectorizedImages*vectorizedImages');
    
    %Determine how many principal components we want to use
    energyPerEigenValue=cumsum(eigenValues(end:-1:1));
    energyPerEigenValue=energyPerEigenValue/energyPerEigenValue(end);
    numberOfComponentsToUse=min(maxNumberOfPrincipalComponents,find(energyPerEigenValue>(1-1/signalToNoise),1,'first')); 
       
    % Show the principle components in 'image' space
    numberOfProbesToShow=9;
    imagesInNormalizedPCBasis=diag(1./sqrt(diag(eigenValues)))*eigenVectors'*vectorizedImages; 
    figure();
    for (probeIdx=1:numberOfProbesToShow)
        subplot(floor(sqrt(numberOfProbesToShow)),ceil(numberOfProbesToShow/floor(sqrt(numberOfProbesToShow))),probeIdx);
        imagesc(reshape(imagesInNormalizedPCBasis(end+1-probeIdx,:),trainingInputSize(1:2)));
        drawnow();
    end
    
    % Calculate the principal components and project the calibration images onto it
    calibration.principalComponentsInImageSpace=eigenVectors'*vectorizedImages;
    calibration.principalComponentsInImageSpace=calibration.principalComponentsInImageSpace(end-numberOfComponentsToUse+1:end,:);
    calibration.trainingImagesInPrincipalComponentSpace=calibration.principalComponentsInImageSpace*vectorizedImages';
    
    %
    % Save the principal components and other relevant data
    %
    logMessage('Saving processed data to file %s.',outputFileName);
    save(outputFileName,'-struct','calibration');

    %
    % Test
    % 
    relativeEigenValues=diag(eigenValues)/sum(diag(eigenValues));
    variabilityConsidered=sum(relativeEigenValues(end-numberOfComponentsToUse+1:end))
    
    %% Check if we have test inputs
    if (nbOfTestImages>0)
        [detectedWavelengths testImagesInPrincipleComponentSpace] = determineWavelengthFromSpeckleImage(testImages,calibration);
         
    	%
    	% Output results
    	%
	    figure();
        deltaLambda=diff(wavelengths(1:2));
        wavelengthIndexMismatch=(detectedWavelengths(:).'-testWavelengths)/deltaLambda;
        detectionEfficiency=sum(wavelengthIndexMismatch<1)/length(wavelengthIndexMismatch)
        hist(wavelengthIndexMismatch(:),50)
    
        figure();
        subplot(221);scatter(testImagesInPrincipleComponentSpace(end,:),testImagesInPrincipleComponentSpace(end-1,:),'+'); title('1-2');
        subplot(222);scatter(testImagesInPrincipleComponentSpace(end,:),testImagesInPrincipleComponentSpace(end-2,:),'+'); title('1-3');
        subplot(223);scatter(testImagesInPrincipleComponentSpace(end-1,:),testImagesInPrincipleComponentSpace(end-2,:),'+'); title('2-3');
        subplot(224);scatter3(testImagesInPrincipleComponentSpace(end,:),testImagesInPrincipleComponentSpace(end-1,:),testImagesInPrincipleComponentSpace(end-2,:),'+'); title('1-2-3');
    end
    %%
end