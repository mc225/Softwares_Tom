%
%
%
function twoParticleRecording()
    fixedParticlePosition=5e-6; % metric units
    movingParticlePositions=-[5:10, 9:6]*1e-6; % metric units
    
    maximumStepSize=.2e-6; % metric units
    
    % Don't go above unity for both combined
    fixedParticlePower=0.35;
    movingParticlePower=0.35;
    
    fixedParticleDeflection=positionToDeflection(fixedParticlePosition);

    % Configure SLM and load aberration correction
    slm=PhaseSLM(2);
    slm.referenceDeflectionFrequency=[1/10 1/10];
    correctionFunctionFileName='C:\Documents and Settings\sk495\My Documents\software\matlab\GUI\calibrateSetup_2013-09-24.mat';
    correctionFileContents=whos('-file',correctionFunctionFileName);
    load(correctionFunctionFileName,'measuredPupilFunction');
    amplificationLimit=1;
    if (any(strcmp(correctionFileContents,'initialCorrection')))
        load(correctionFunctionFileName,'initialCorrection');
    else
        initialCorrection=1;
    end
    slm.correctionFunction=calcCorrectionFromPupilFunction(measuredPupilFunction./initialCorrection,amplificationLimit);
    logMessage('Correction function updated using amplification limit %d.',amplificationLimit);
                
    movingParticleDeflection=positionToDeflection(movingParticlePositions(1));
    slm.modulate(@(X,Y) (sqrt(X.^2+Y.^2)<300).*(fixedParticlePower*(exp(2i*pi*fixedParticleDeflection*X)) + movingParticlePower*(exp(2i*pi*movingParticleDeflection*X))));
    previousMovingParticlePosition=movingParticlePositions(1);
    
    % Loop until CTRL-C
    stopMeasurement=false;
    while ~stopMeasurement
        % Loop through all particle positions
        for posIdx=[1:length(movingParticlePositions)]
            targetPosition=movingParticlePositions(posIdx);
            logMessage('Moving particle 2 from %0.3f to %0.3f um.',[previousMovingParticlePosition targetPosition]*1e6);
            % Move slowly
            for movingParticlePosition=previousMovingParticlePosition:sign(targetPosition-previousMovingParticlePosition)*maximumStepSize:targetPosition
                movingParticleDeflection=positionToDeflection(movingParticlePosition);
                slm.modulate(@(X,Y) (sqrt(X.^2+Y.^2)<300).*(fixedParticlePower*(exp(2i*pi*fixedParticleDeflection*X)) + movingParticlePower*(exp(2i*pi*movingParticleDeflection*X))));
            end
            % Go to target position
            movingParticleDeflection=positionToDeflection(targetPosition);
            slm.modulate(@(X,Y) (sqrt(X.^2+Y.^2)<300).*(fixedParticlePower*(exp(2i*pi*fixedParticleDeflection*X)) + movingParticlePower*(exp(2i*pi*movingParticleDeflection*X))));
            previousMovingParticlePosition=targetPosition;
            % Record what needs to be done
            result=recordData(movingParticlePositions(posIdx));
            if ~result
                stopMeasurement=true;
            end
        end
        % Wait a bit for next measurement
        logMessage('Waiting for next movement.');
        pause(5);
    end

end

function result=recordData(movingParticlePosition)
    logMessage('Recording data...');
    pause(2);
    result=true;
end

function  deflection=positionToDeflection(position)
    slmPixelPitch=20e-6;
    telescopeMagnifications=[100/175 250/250];
    objectiveMagnification=100;
    objectiveTubeLensFocalLength=200e-3;
    wavelength=1070e-9;
    
    objectiveFocalLength=objectiveTubeLensFocalLength/objectiveMagnification;
    slmPixelPitchAtBackAperture=prod(telescopeMagnifications)*slmPixelPitch;
    tanOfUnityDeflection=wavelength/slmPixelPitchAtBackAperture;
    focalShiftOfUnityDeflection=objectiveFocalLength*tanOfUnityDeflection;
   
    deflection=position/focalShiftOfUnityDeflection;
end