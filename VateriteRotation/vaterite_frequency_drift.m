% Signal analysis for vaterite experiments: measuring the time-dependent rotation frequency and optical coupling
% MRDSV, Oct 2012

clear all; close all; clc;
cd 'D:\Matthieu\...';      % Directory
% addpath 'D:\Matthieu\...';

% Data are saved in 3 column text files, where col 1 is time, col 2 is channel 1, col 3 is channel 2.
% Files names follow the root 'yXXz.txt', where e.g. y = prefix, XX the laser intensity (/100 A), and z = suffix.

Intensity = 200;                            % Type the laser pump current (/100 A)
Root = Intensity;
Prefix = 'B';
Suffix = 'E';
HWP = 15;                                   % HWP angle (default: 55°)
Power = 0.55*(53.415 - 44.254*sind(4*HWP+47.51))*(918.13*Intensity/100 - 374.92)/95.6;  % Laser power in sample (mW)

% Parameters

NP = 1;                     % Evaluating both p and np files? If 0, uses the imposed value of Sigma_np
Sig_NP = -0.786;            % Imposed value of Sigma_np
TEST = 0;                   % Displaying all found extrema? (for test purposes only)
run = 9;                    % Number of points on which the running mean is calculated (must be odd por practical reasons)
Max = 1;                    % Choose 1 (0) for looking for maxima (minima) of Signal 1 (and opposite for Signal 2)
Noisy = 0;                  % If a signal is too noisy type its index (1 or 2), the extrema will be looked for on the other one
Saw = 0;                    % If the signal is sawtooth-shaped type 1 (all other parameters will be forced as below)
if Saw == 1
    run = 1;
    Max = 1;
    Noisy = 1;
end

% Extracting the data

R = num2str(Root);
eval(['load ' Prefix R Suffix '.txt' ';']);
Data = eval([Prefix R Suffix]);

Time = Data(:,1);
Voltage1 = Data(:,2);
Voltage2 = Data(:,3);
Signal1 = Voltage1(Voltage1 > 0);          % Removes the negative data points (meaningless)
Signal2 = Voltage2(Voltage2 > 0);
Time = Time(Voltage1 > 0);

% Smoothing the signal and computing its variation

RunningMean1 = tsmovavg(Signal1','s',run);                              % Running mean over 'run' data points
RunningMean1 = RunningMean1(find(RunningMean1 >=0 | RunningMean1 < 0)); % Removes the 'NaN' points
RunningMean2 = tsmovavg(Signal2','s',run);
RunningMean2 = RunningMean2(find(RunningMean2 >=0 | RunningMean2 < 0));

dRM1 = [1 diff(RunningMean1)];              % Differentiation of the signal running mean
dRM2 = [1 diff(RunningMean2)];

% Finding the extrema (i.e. where the differential = 0)

% For convenience I chose to look for the maxima (and minima) of signal 1 (and 2). The conditions are:
% (i) the derivative changes sign once (and only) before and after the investigated data point;
% (ii) the signal value is bigger (smaller) than a threshold defined from the local values of the signal ('25' points before and after);
% (iii) each extremum must be followed by at least one opposite.

Z1 = [];
Z2 = [];
Zeros1 = [];
Zeros2 = [];
for n = 25:length(dRM1)-25
    Zmax1 = n( find( dRM1(n-1) >= 0 && dRM1(n) >= 0 && dRM1(n+1) < 0 ... 
        && RunningMean1(n) > (max(RunningMean1(n-24:n+25))+min(RunningMean1(n-24:n+25)))/2 )) + floor(run/2);
    Zmin1 = n( find( dRM1(n-1) <= 0 && dRM1(n) <= 0 && dRM1(n+1) > 0 ... 
        && RunningMean1(n) < mean(RunningMean1(n-24:n+25)) )) + floor(run/2);
    if Saw ~= 1
        Zmin2 = n( find( dRM2(n-1) <= 0 && dRM2(n) <= 0 && dRM2(n+1) > 0 ...
            && RunningMean2(n) < mean(RunningMean2(n-24:n+25))/1.0 )) + floor(run/2);
        Zmax2 = n( find( dRM2(n-1) >= 0 && dRM2(n) >= 0 && dRM2(n+1) < 0 ... 
            && RunningMean2(n) > (max(RunningMean2(n-24:n+25))+min(RunningMean2(n-24:n+25)))/2 )) + floor(run/2);
    else
%         Zmin2 = n( find( dRM2(n) < 0 && dRM2(n+1) >= 0 )) + floor(run/2);
        Zmin2 = n( find( dRM2(n) < 0 && dRM2(n+1) >= 0 ...
            && RunningMean2(n) < mean(RunningMean1(n-24:n+25))/4 )) + floor(run/2);
        Zmax2 = n( find( dRM2(n) >= 0 && dRM2(n+1) <= 0 )) + floor(run/2);
    end
    Z1 = [Z1; (-1)^(Max+1)*Zmax1; (-1)^Max*Zmin1];
    Z2 = [Z2; (-1)^(Max+1)*Zmin2; (-1)^Max*Zmax2];
end
if Noisy > 0
    eval(['Z' num2str(Noisy) '= Z' num2str(real(2 + i^(Noisy))) ';']);
end

for k = 1:length(Z1)-1
    z = k(find(Z1(k) > 0 && Z1(k+1) < 0));
    Zeros1 = [Zeros1;Z1(z)];
end
for j = 1:length(Z2)-1
    zz = j(find(Z2(j) > 0 && Z2(j+1) < 0));
    Zeros2 = [Zeros2;Z2(zz)];
end

% Converting the extrema to the particle rotation time, and determining the period

if TEST == 1
    Rotation1 = Zeros1;                 % For practical reasons I keep the possibility of displaying all found extrema
    Rotation2 = Zeros2;
else
    Rotation1 = Zeros1(mod(1:length(Zeros1),2) > 0);    % I arbitrarily choose the odd indexes of the 'Zeros' vectors
    Rotation2 = Zeros2(mod(1:length(Zeros2),2) > 0);
end

% Small correction for the (rare) case of an extremum found on only one signal (close to temporal cutoff)
if length(Rotation1) ~= length(Rotation2)
    Rotation1 = Rotation1(1:min(length(Rotation1),length(Rotation2)));
    Rotation2 = Rotation2(1:min(length(Rotation1),length(Rotation2)));
end
RotationTime1 = Time(Rotation1);
RotationTime2 = Time(Rotation2);
Period1 = diff(RotationTime1);
Period2 = diff(RotationTime2);
Frequency1 = 1./Period1;
Frequency2 = 1./Period2;

% Depicting results

figure
hold on
plot(Time,Signal1, 'r')
plot(Time(round(run/2):length(RunningMean1)+run-round(run/2)),RunningMean1', 'g')
plot(RotationTime1,Signal1(Rotation1), '+k')
plot(Time,Signal2, 'b')
plot(Time(round(run/2):length(RunningMean2)+run-round(run/2)),RunningMean2', 'y')
plot(RotationTime2,Signal2(Rotation2), 'ok')
plot(Time,Signal1+Signal2, '-.c')
title('Detection of the rotating periods')

% Evaluation of the temporal variation of DeltaSigma (from period to period)

if TEST ~= 1
    for t = 1:length(Period1)
        RunningSignal1 = Signal1(Rotation1(t):Rotation1(t+1));
        RunningSignal2 = Signal2(Rotation2(t):Rotation2(t+1));
        RunningInt1(t) = sum(RunningSignal1);
        RunningInt2(t) = sum(RunningSignal2);
        RunningSigma(t) = (RunningInt2(t) - RunningInt1(t))/(RunningInt2(t) + RunningInt1(t));

        if NP == 1
            eval(['load np' num2str(Intensity) '.txt' ';']);
            Data_np = eval(['np' num2str(Intensity)]);
            Time_np = Data_np(:,1);
            Voltage1_np = Data_np(:,2);
            Voltage2_np = Data_np(:,3);
            Signal1_np = Voltage1_np(Voltage1_np > 0);
            Signal2_np = Voltage2_np(Voltage2_np > 0);
            Time_np = Time_np(Voltage1_np > 0);
            Integral1_np = sum(Signal1_np);
            Integral2_np = sum(Signal2_np);
            sigma_np = (Integral2_np - Integral1_np)/(Integral2_np + Integral1_np);
        else
            sigma_np = Sig_NP;
        end
        RunningDeltaSigma(t) = sigma_np - RunningSigma(t);
        OpticalTorque(t) = abs(Power*RunningDeltaSigma(t)/1.762);
    end
    Integral1 = sum(Signal1);
    Integral2 = sum(Signal2);
    sigma = (Integral2 - Integral1)/(Integral2 + Integral1);
    DeltaSigma = sigma_np - sigma;

    % Displaying results

    figure
    subplot(1,2,1)
    hold on
    plot(RotationTime1,[NaN;Frequency1], '+-k')
    plot(RotationTime2,[NaN;Frequency2], 'o-b')
    title('Evolution of the rotation frequency')
    legend('Ch 1', 'Ch 2', 'location','SouthEast')
    subplot(1,2,2)
    plot(RotationTime1,[NaN RunningDeltaSigma], [min(RotationTime1) max(RotationTime1)],[DeltaSigma DeltaSigma])
    title('Evolution of the optical coupling')

    figure
    hold on
    plot(OpticalTorque,Frequency1, '+-k')
    plot(OpticalTorque,Frequency2, 'o-b')
    title('Rotation frequency vs optical torque')
    legend('Ch 1', 'Ch 2', 'location','SouthEast')

    Results = [RotationTime1(2:length(RotationTime1)) RunningDeltaSigma' OpticalTorque' Frequency1];
end