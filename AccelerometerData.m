clc
clear
close all

% Load in walking data
load('/Users/yang/MATLAB-Drive/MobileSensorData/walking1010.mat')
time = Acceleration.Timestamp;

% Conversion of date-time format to seconds starting from 0
initial_min = double(minute(time(1)));
initial_sec = double(second(time(1)));

for i=1:length(time)
    if double(minute(time(i)))>initial_min
        t(i) = double(second(time(i))) + 60 - initial_sec;
    else
        t(i) = double(second(time(i)))-initial_sec;
    end
end

% Removing constant bias and gravitational acceleration
accelX = Acceleration.X + 0.1;
accelY = Acceleration.Y + 0.03;
accelZ = Acceleration.Z - 9.85;


plot(t,accelX,'r',t,accelY,'g',t,accelZ,'b')
legend("Acceleration X", "Acceleration Y", "Acceleration Z")
xlabel("Time (s)")
ylabel("Acceleration (m/s^2)")
title("3 Axis Acceleration Data")

accelNorm = sqrt(accelX.^2 + accelY.^2 + accelZ.^2);

fs = 50; %Sampling Frequency
figure();
plot(t, accelNorm)
xlabel("Time (s)")
ylabel("Acceleration (m/s^2)")
title("Magnitude of Acceleration")


% Take FFT of the discrete time signal and shift frequency spectral density
% according to the sampling rate
y = fft(accelNorm, length(accelNorm));
y = fftshift(y);
f = (-length(y)/2:(length(y)-1)/2)*fs/length(y);

m = abs(y);
realfreq=f(f>0);
realmag=m(f>0);
figure();
plot(realfreq,realmag)
title("Power Spectral Density")
xlabel("Frequency (Hz)")
ylabel("Amplitude")

% Find all the local peaks of the Frequency Spectrum
[peaks, freqPeaks] = findpeaks(realmag);
maxPeak = 0;
maxFreqIndex = 0;

% Set minimum and maximum frequency ranges to avoid constant gravity noise
% and high frequency noises
minFreq = 0.5;
maxFreq = 5;

% Iterate through Frequency peaks to find the largest amplitude (ignoring
% the constant gravity around 0Hz) to find dominant walking cadence
for i = 1:length(freqPeaks)-1
    if peaks(i) > maxPeak && realfreq(freqPeaks(i)) > minFreq && realfreq(freqPeaks(i)) < maxFreq
        maxPeak = peaks(i);
        maxFreqIndex = freqPeaks(i);
    end
end
hold on
plot (realfreq(maxFreqIndex), maxPeak, 'o')

% Set a range of frequencies to take according to the dominant cadence
freqWidth = 0.25;
freqRange = [realfreq(maxFreqIndex) - realfreq(maxFreqIndex)*freqWidth, realfreq(maxFreqIndex) + realfreq(maxFreqIndex)*freqWidth];
if freqRange(1)<minFreq
    freqRange(1) = minFreq;
end

if freqRange(2)>maxFreq
    freqRange(2) = maxFreq;
end

% Plot the range of frequency allowed in the bandpass filter
hold on;
ylimits = double(ylim);
a = area([freqRange(1), freqRange(1), freqRange(2), freqRange(2)], [0, ylimits(2), ylimits(2), 0]);
a.FaceAlpha = 0.2;

% 6th order bandpass filter for the determined frequency range
[b,a] = butter(6, freqRange/(fs/2), 'bandpass');
y=filter(b,a,accelNorm);
figure();
plot(t, y)
title("Filtered Acceleration Magnitude")
xlabel("Time (s)")
ylabel("Acceleration (m/s^2)")

% Calculate the minimum peaks to determine a step based on the quarter of 
% the measured maximum acceleration magnitude
[potSteps, potTimeStep] = findpeaks(y);
sortedAccelerations = sort(potSteps, 'descend');
accelTopQuarter = sortedAccelerations(1:round(length(sortedAccelerations)/4));
minPeak = 1;

% Set a lower minimum amplitude to account for lower acceleration 
% magnitudes when taking a step with the leg not containing the phone
minPeakMeasured = 0.5*sum(accelTopQuarter)/length(accelTopQuarter);

% Ensure the calculated minimum peak is not lower than lower bound to
% account for time instances where no steps were taken
if minPeakMeasured < minPeak
    minPeakMeasured = minPeak;
end

% Iterate through the filtered signal's local peaks to find steps
j = 1;
for i = 1:length(potSteps)-1
    if potSteps(i) > minPeakMeasured
        steps(j) = potSteps(i);
        timeStep(j) = potTimeStep(i);
        j = j+1;
    end
end

% Plot the detected steps on the filtered time signal
numSteps = length(steps);
fprintf("The pedometer has detected %d steps\n", numSteps);
hold on
plot(t(timeStep), steps, 'o');
