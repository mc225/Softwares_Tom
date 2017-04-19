% createResolutionPlots(sourceFolder)
%
%
function createResolutionPlots(sourceFolder)
    if nargin<1 || isempty(sourceFolder)
        sourceFolder=uigetdir('.','Select Data Folder...');
        if (isempty(sourceFolder))
            logMessage('User cancelled');
            return;
        end
    end
    
    sliceWidth=10e-6;
    perspectiveScaling=[0 500]; %[373.4900  722.6800];
    scanShear=[0 .16];
    
    folderDescs=dir(sourceFolder);
    for folderIdx=1:length(folderDescs)
        folderDesc=folderDescs(folderIdx);
        if folderDesc.isdir && ~ismember(folderDesc.name,{'.','..'})
            tokens=regexp(folderDesc.name,'(\d\d\d\d-\d\d-\d\d\ \d\d_\d\d_\d\d\.\d\d\d )?at(_?)(\d+)um$','tokens');
            if ~isempty(tokens)
                approximateOffset=str2double(tokens{1}{3})*1e-6;
                if ~strcmp(tokens{1}{2},'_')
                    approximateOffset=-approximateOffset;
                end
                inputFolder=fullfile(sourceFolder,folderDesc.name);
                targetFolder=fullfile([sourceFolder,'_sliced'],folderDesc.name);
                logMessage('Creating folder %s...',targetFolder);
                mkdir(targetFolder);
                fileDescs=dir(fullfile(inputFolder,'recording_*.mat'));
                inputFileNames=arrayfun(@(fn) fullfile(inputFolder,fn.name),fileDescs,'UniformOutput',false);
                outputFileNames=arrayfun(@(fn) fullfile(targetFolder,fn.name),fileDescs,'UniformOutput',false);
                for fileIdx=1:length(inputFileNames)
                    inputMatFile=matfile(inputFileNames{fileIdx},'Writable',false);
                    % Prepare slicing
                    limits=approximateOffset+[-0.5 0.5]*sliceWidth;
                    xRange=inputMatFile.xRange;
                    yRange=inputMatFile.yRange+63e-6;
                    zRange=inputMatFile.zRange;
                    selectedIndexes=find(yRange>=limits(1) & yRange<=limits(2));
                    %Slice
                    yRange=yRange(1,selectedIndexes);
                    recordedImageStack=inputMatFile.recordedImageStack(:,selectedIndexes,:);
                    restoredDataCube=inputMatFile.restoredDataCube(:,selectedIndexes,:);
                    lightSheetPsf=inputMatFile.lightSheetPsf(1,selectedIndexes,:);
                    ZOtf=inputMatFile.ZOtf;
                    setupConfig=inputMatFile.setupConfig;
                    delete(inputMatFile);

                    effectiveNA=setupConfig.excitation.fractionOfNumericalApertureUsed*setupConfig.excitation.objective.numericalAperture;
                    spFreqCutOff=(2*effectiveNA)/setupConfig.excitation.wavelength;

                    recordedImageStack=geometricCorrection(scanShear,perspectiveScaling,xRange,yRange,zRange,recordedImageStack);

                    recordedImageStack=recordedImageStack-median(recordedImageStack(:)); %Subtract background noise

                    % Find maximum intensity
                    axialProj=sum(recordedImageStack,3);
                    [~, maxI]=max(axialProj(:));
                    [maxRow maxCol]=ind2sub(size(axialProj),maxI);
                    %Select some pixels around it and integrate
                    %intensityTrace=recordedImageStack(maxRow+[-1:1],maxCol+[-1:1],:);
                    intensityTrace=recordedImageStack(:,:,:);
                    intensityTrace=squeeze(mean(mean(intensityTrace)));
                    % Calculate MTF
                    mtf=fft(intensityTrace([1:end end*ones(1,floor(end/2)) 1*ones(1,floor((1+end)/2))]));
                    mtf=mtf./mtf(1);
                    mtf=fftshift(abs(mtf));

                    fig=figure('Position',[50 50 1024 768]);
                    subplot(2,2,1);
                    plot(zRange*1e6,intensityTrace); xlabel('z [um]');
                    subplot(2,2,2);
                    plot(ZOtf*1e-3,mtf);
                    xlim([0 spFreqCutOff*1e-3]);
                    xlabel('\nu_z [cycles/mm]');
%                         showImage(squeeze(max(recordedImageStack,[],1)).',-1,yRange*1e6,zRange*1e6);
                    subplot(2,2,3);
                    imagesc(yRange*1e6,zRange*1e6,squeeze(max(recordedImageStack,[],1)).');
                    axis equal; xlabel('x [\mum]'); ylabel('z [\mum]');
                    subplot(2,2,4);
                    imagesc(xRange*1e6,zRange*1e6,squeeze(max(recordedImageStack,[],2)).');
                    axis equal; xlabel('y [\mum]'); ylabel('z [\mum]');
                    title(inputFileNames{fileIdx});
                    saveas(fig,[inputFileNames{fileIdx}(1:end-4),'.png'],'png');
                    close(fig);
                end
            end
        end
    end
end
    
function recordedImageStack=geometricCorrection(scanShear,scaling,xRange,yRange,zRange,recordedImageStack)
    cubicInterpolation=true;
    
    if (~all(scanShear==0) || ~all(scaling==0))
        % Geometrically correcting recorded data cube and light sheet
        logMessage('Geometrically shifting and deforming recorded date cube and light sheet by [%0.3f%%,%0.3f%%] and a magnification change of [%0.3f%%,%0.3f%%] per micrometer...',-[scanShear scaling*1e-6]*100);
        for (zIdx=1:size(recordedImageStack,3))           
            zPos=zRange(zIdx);
            sampleXRange = scanShear(1)*zPos + xRange*(1-scaling(1)*zPos);
            sampleYRange = scanShear(2)*zPos + yRange*(1-scaling(2)*zPos);
            if (cubicInterpolation)
                interpolatedSlice=interp2(yRange.',xRange,recordedImageStack(:,:,zIdx),sampleYRange.',sampleXRange,'*cubic',0);
            else
                interpolatedSlice=interp2(yRange.',xRange,recordedImageStack(:,:,zIdx),sampleYRange.',sampleXRange,'*linear',0);
           end
            recordedImageStack(:,:,zIdx)=interpolatedSlice;
        end
    end
end
    