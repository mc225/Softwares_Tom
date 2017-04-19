% Signal analysis for vaterite experiments: calculating the signal integral of the 2 polarization components and extracting the particle
% rotation frequency.
% MRDSV, July 2012

clear all; close all; clc;
cd 'D:\Matthieu\...';      % Directory
% addpath 'D:\Matthieu\...';

% Data are saved in 3 column text files, where col 1 is time, col 2 is channel 1, col 3 is channel 2.
% Files names follow the root 'yXXz.txt', where e.g. y = prefix, XX the laser intensity (/100 A), and z = suffix.
% 'np' files correspond to the cases without a particle (for signal calibration).

Intensity = 80;                             % Type the laser pump current (/100 A)
Root = Intensity;
Prefix = 'c';
Suffix = '';
HWP = 29;                                   % HWP angle (default: 55°)
Power = 0.55*(53.415 - 44.254*sind(4*HWP+47.51))*(918.13*Intensity/100 - 374.92)/95.6;  % Laser power in sample (mW)
NP = 1;                                     % Evaluating both p and np files? If 0, uses the imposed value of Sigma_np
Sig_NP = -0.781;                            % Imposed value of Sigma_np
Fmin = 0.05;                                % Limits of the spectrum to be displayed
Fmax =04.0;

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
dt = diff(Time);
SumSignal = Signal1 + Signal2;

% Plotting the signal
figure
plot(Time,Signal1, Time,Signal2, Time,SumSignal,'--y');

% Calculating the integral - approximated as a discrete sum as the signal is discrete.
% Note: the time increment is actually not required as it cancels in sigma (ratio).
Integral1 = sum(Signal1)*dt(1);
Integral2 = sum(Signal2)*dt(1);
sigma = (Integral2 - Integral1)/(Integral2 + Integral1);

if NP == 1
    eval(['load np' num2str(Intensity) '.txt' ';']);
    Data_np = eval(['np' num2str(Intensity)]);
    Time_np = Data_np(:,1);
    Voltage1_np = Data_np(:,2);
    Voltage2_np = Data_np(:,3);
    Signal1_np = Voltage1_np(Voltage1_np > 0);
    Signal2_np = Voltage2_np(Voltage2_np > 0);
    Time_np = Time_np(Voltage1_np > 0);
    dt_np = diff(Time_np);
    Integral1_np = sum(Signal1_np)*dt_np(1);
    Integral2_np = sum(Signal2_np)*dt_np(1);
    sigma_np = (Integral2_np - Integral1_np)/(Integral2_np + Integral1_np);
    figure
    plot(Time_np,Signal1_np, Time_np,Signal2_np);
else
    sigma_np = Sig_NP;
end
Delta_sigma = sigma_np - sigma;
OpticalTorque = abs(Power*Delta_sigma/1.762);           % given in pN µm
OpticalData = [Power mean(SumSignal) Delta_sigma OpticalTorque];
OtherData = [sigma sigma_np mean(Signal1_np) mean(Signal2_np)];

% Fourier spectrum
Fsample = 1/dt(1);
NFFT = 2^nextpow2(length(Time));
f = Fsample/2*linspace(0,1,NFFT/2+1);
F1 = fft(Signal1,NFFT)/length(Time);
F2 = fft(Signal2,NFFT)/length(Time);

% Plotting the spectrum
figure
plot(f,2*abs(F1(1:NFFT/2+1)), f,2*abs(F2(1:NFFT/2+1)))
xlim([Fmin Fmax])