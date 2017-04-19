%
% previewLightSheetRecordings(folders)
%
% Calls previewLightSheetRecording repeatedly
%
function previewLightSheetRecordings(folderNames,destinationEmails)
	if (nargin<1 || isempty(folderNames))
        folderNames={'H:\RESULTS\2013-11-27\Amphioxus'};
%         folderNames={'H:\RESULTS\2013-11-27\Zebra fish','H:\RESULTS\2013-11-27\Spheroid'}; % Jonny: Please update the e-mail addresses below
% 		folderNames=uigetdir(pwd(),'Select Recording(s)');
%         if (~iscell(folderNames) && ~ischar(folderNames) && ~isempty(folderNames))
%             logMessage('Folder selection canceled.');
%             return;
%         end
    end
    if (~iscell(folderNames))
        folderNames={folderNames};
    end
    if nargin<2
        destinationEmails={'tv2@st-andrews.ac.uk','ccl9@st-andrews.ac.uk','jn78@st-andrews.ac.uk'};
    end
    
    dataFileRegEx='recording_lambda([0-9]+)nm.+\.avi';
    
    % Recursively process all folders
    for folderName=folderNames,
        folderName=folderName{1};
        allFilesOrDirs=dir(folderName);
        allDirs={allFilesOrDirs(cell2mat({allFilesOrDirs.isdir}) & cellfun(@(nm) nm(1)~='.',{allFilesOrDirs.name})).name};
        allDataFiles={allFilesOrDirs(cellfun(@(nm) ~isempty(regexpi(nm,dataFileRegEx)),{allFilesOrDirs.name})).name};
        allDirs=cellfun(@(nm) fullfile(folderName,nm),allDirs,'UniformOutput',false);
        allDataFiles=cellfun(@(nm) fullfile(folderName,nm),allDataFiles,'UniformOutput',false);
        % Go through subfolders and add to file list if it has a time stamp, otherwise recurse
        for subDir=allDirs,
            subDir=subDir{1};
            if ~isempty(regexpi(subDir,'[0-9][0-9][0-9][0-9]+-[0-9][0-9]-[0-9][0-9] [0-9][0-9]_[0-9][0-9]_[0-9][0-9]\.[0-9][0-9]+$')) % If a time stamp
                subFolderContents=dir(subDir);
                subFolderContents={subFolderContents(cellfun(@(nm) ~isempty(regexpi(nm,dataFileRegEx)),{subFolderContents.name})).name};
                subFolderContents=cellfun(@(nm) fullfile(subDir,nm),subFolderContents,'UniformOutput',false);
                allDataFiles(end+[1:numel(subFolderContents)])=subFolderContents; % Append new found files
            else
                previewLightSheetRecordings(subDir,destinationEmails);
            end
        end
        % Combine files by wavelength
        allDataFileSets={};
        for dataFileIdx=1:numel(allDataFiles),
            dataFileSet={allDataFiles{dataFileIdx}}; % Still a single file at this point
            % Check if there is also another wavelength version
            wavelengthStr=regexpi(dataFileSet{1},dataFileRegEx,'tokens');
            wavelengthStr=wavelengthStr{1}{1};
            firstGreen=strcmp(wavelengthStr,'488');
            if firstGreen,
                otherWavelengthStr='532';
            else
                otherWavelengthStr='488';
            end
            fileToFind=strrep(dataFileSet{1},['lambda',wavelengthStr,'nm'],['lambda',otherWavelengthStr,'nm']);
            clear otherWavelengthStr;
            [~,fileToFind]=fileparts(fileToFind); %Use only the file name part
            otherFileIdx=find(cellfun(@(nm) ~isempty(regexpi(nm,fileToFind)),allDataFiles),1);
            if ~isempty(otherFileIdx)
                if otherFileIdx>dataFileIdx
                    dataFileSet{2}=allDataFiles{otherFileIdx};
                    if ~firstGreen && numel(dataFileSet)>=2, % Make sure that green is first if it exists
                        dataFileSet=dataFileSet([2 1]);
                    end
                    allDataFileSets{end+1}=dataFileSet;
                end
            else
                allDataFileSets{end+1}=dataFileSet;
            end
        end
        % Process all the collected data files sets
        for dataFileSetIdx=1:numel(allDataFileSets),
            dataFileSet=allDataFileSets{dataFileSetIdx};
            try
                outputFullFileName=previewLightSheetRecording(dataFileSet);
                logMessage('Just created %s.',outputFullFileName,destinationEmails,outputFullFileName);
            catch Exc
                Exc
                logMessage('Could not process %s.',dataFileSet{1});
            end
        end
    end
    
end