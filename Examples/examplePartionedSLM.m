%
% examplePartionedSLM()
%
% This example shows how to create two SLM objects that use each a separate part of the same SLM.
% Both parts have different first order deflections and behave entirely independent.
% This method may be used to build a set-up with two spatial light modulators using only one physical device.
% Applications include:
% * horizontal and vertical polarization modulation
% * dispersion corrected modulation of short pulses
%
function examplePartionedSLM()
    close all;
    displayNumber=[]; %Set to [] to test with a Matlab window, or 2 to test with screen 2
    
    %Create a fake screen for testing if displayNumber==[]
    if (isempty(displayNumber))
        hFig=figure();
        displayNumber=axes('Parent',hFig);
        image(zeros(600,800),'Parent',displayNumber);
    end
    
    %Create three SLM objects that use the same output device
    leftSLM=PhaseSLM(displayNumber,[1/10 1/10]);
    leftSLM.regionOfInterest=leftSLM.regionOfInterest.*[1 1 1 0.5];
    leftSLM.stabilizationTime=0.0; %don't keep us waiting
    rightSLM=PhaseSLM(displayNumber,[-1/10 -1/10]);
    rightSLM.regionOfInterest=rightSLM.regionOfInterest.*[0 0 1 0.5]+[0 0.5*rightSLM.regionOfInterest(4) 0 0];
    rightSLM.stabilizationTime=0.0; %don't keep us waiting

    %Modulate both parts of the SLM separately
    for innerRadius=[0:10:190],
        leftSLM.modulate(parsePupilEquation('R<200');
        rightSLM.modulate(parsePupilEquation(sprintf('R<200 & R>%f',innerRadius));
    end

    % Close everything again
    pause(5);
    close all;
    DefaultScreenManager.instance.delete(); % Close all full-screens
end
