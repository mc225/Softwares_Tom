% copyMatFilesOnly(sourceFolder,targetFolder)
%
% Copies only the .mat files and removes the date-time folders, combining everything into the parent folder.
%
function status=copyMatFilesOnly(sourceFolder,targetFolder)
    if nargin<1 || isempty(sourceFolder)
        sourceFolder=uigetdir(pwd(),'Select Source Folder...');
        if (sourceFolder==0)
            logMessage('No source folder selected.');
            return;
        end
    end
    if nargin<2 || isempty(targetFolder)
        targetFolder=uigetdir('D:\Stored Files\RESULTS','Select Destination...');
        if (targetFolder==0)
            logMessage('No target folder selected.');
            return;
        end
    end
    
    status=true;
    
    % Copy all files 
    fileNameList=dir(fullfile(sourceFolder,'*.mat'));
    for (fileName={fileNameList.name})
        fileName=fileName{1};
        filePathAndName=fullfile(sourceFolder,fileName);
        logMessage(['Copying ',filePathAndName,' to ',targetFolder,'...']);
        if ~exist(fullfile(targetFolder,fileName),'file')
            status=copyfile(filePathAndName,targetFolder);
            if (~status)
                logMessage('Error copying file %s!',filePathAndName);
                return;
            end
        else
            logMessage('The file %s is already copied.',fullfile(targetFolder,fileName));
        end
    end
    % Recurse through folders
    fileDescList=dir(sourceFolder);
    for fileDescIdx=1:length(fileDescList)
        fileDesc=fileDescList(fileDescIdx);
        if fileDesc.isdir && ~strcmp(fileDesc.name,'.') && ~strcmp(fileDesc.name,'..')
            if regexp(fileDesc.name,'^\d\d\d\d-\d\d-\d\d\ \d\d_\d\d_\d\d.\d\d\d$')
                targetFolderRecursive=targetFolder;
            else
                targetFolderRecursive=fullfile(targetFolder,fileDesc.name);
                [status message]=mkdir(targetFolderRecursive);
                if ~status
                    logMessage(message);
                end
            end
            status=copyMatFilesOnly(fullfile(sourceFolder,fileDesc.name),targetFolderRecursive);
            if (~status)
                return;
            end
        end
    end
    
    % Results
    if nargout==0
        if status
            logMessage('Copied all files.');
        else
            logMessage('Could not copy (all) files.');
        end
        clear status;
    end
end