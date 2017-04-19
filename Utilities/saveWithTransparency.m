% saveWithTransparency(figHandle,fileName)
%
% Saves a figure as a png image with a transparency
%
function saveWithTransparency(figHandle,fileName)
    if (~strcmpi(fileName(end-3:end),'.png'))
        fileName=strcat(fileName,'.png');
    end
    % save the original settings
    oldBackGround = get(figHandle,'Color');
    invertHardCopyStatus=get(figHandle,'InvertHardCopy');
    
    set(figHandle,'InvertHardCopy','off');
    
    % specify an all-black background and record the image
    set(figHandle,'Color',[0 0 0]);
    print(figHandle,'-dpng','-r300',fileName);
    noColorImg = imread(fileName);
    % Specify an all-white background and record the image
    set(figHandle,'Color',[1 1 1]);
    print(figHandle,'-dpng','-r300',fileName);
    maxColorImg = imread(fileName);
    
    % Calculate the alpha value from the modulation
    alpha=maxColorImg-noColorImg;
    
    % Write the image with alpha channel
    imwrite(noColorImg,fileName, 'png', 'BitDepth', 16,'Alpha',1.0-double(alpha)./256);
    
    set(figHandle,'Color',oldBackGround);
    set(figHandle,'InvertHardCopy',invertHardCopyStatus);
end