
close all

% Load in walking data
time = Acceleration.Timestamp;

% Conversion of date-time format to seconds starting from 0
initial_min = double(minute(time(1)));
initial_sec = double(second(time(1)));
t = zeros(length(time));

for i=1:length(time)
    if double(minute(time(i)))>initial_min
        t(i) = double(second(time(i))) + 60 - initial_sec;
    else
        t(i) = double(second(time(i)))-initial_sec;
    end
end

% Removing constant bias and gravitational acceleration
accelX = Acceleration.X + 0.05;
accelY = Acceleration.Y + 10;
accelZ = Acceleration.Z + 0.07;

%plot(t,accelZ)
p1 = plot(t,accelX,'r');
hold on 
p2 = plot(t,accelY,'g');
p3 = plot(t,accelZ,'b');
hold off
legend([p1(1) p2(1) p3(1)], {'Acceleration X', 'Acceleration Y', 'Acceleration Z'})
xlabel("Time (s)")
ylabel("Acceleration (m/s^2)")
title("3 Axis Acceleration Data")
%title("Z-Axis Acceleration Data")
accelNorm = sqrt(accelX.^2 + accelY.^2 + accelZ.^2);

fs = 50; %Sampling Frequency (higher than the nyquist frequency)
figure();
plot(t, accelNorm)
xlabel("Time (s)")
ylabel("Acceleration (m/s^2)")
title("Magnitude of Acceleration")

% Set minimum and maximum frequency ranges to avoid constant gravity noise
% and high frequency noises
minFreq = 1;
maxFreq = 4;

% Take FFT of the discrete time signal and shift frequency spectral density
% according to the sampling rate
Y = fft(accelNorm, length(accelNorm));
Y = fftshift(Y);
f = (-length(Y)/2:(length(Y)-1)/2)*fs/length(Y); % zero-centered frequency range
m = abs(Y).^2/length(Y); % zero-centered power
figure();
plot(f,m)
xlim([minFreq maxFreq])
title("Relative Amplitude vs Frequency")
xlabel("Frequency (Hz)")
ylabel("Amplitude")

% Find Power Spectral Density 
r = xcorr(accelNorm);
PSD = fft(r);
PSD = fftshift(PSD);
m = real(PSD);
f = (-length(PSD)/2:(length(PSD)-1)/2)*fs/length(PSD);
figure();
plot(f,m);
xlim([minFreq maxFreq])
title("Power Spectral Density")
xlabel("Frequency (Hz)")
ylabel("Power Density (dB/Hz)")

% Find all the local peaks of the Frequency Spectrum
[peaks, freqPeaks] = findpeaks(m);
maxPeak = 0;
maxFreqIndex = 0;



% Iterate through Frequency peaks to find the largest amplitude (ignoring
% the constant gravity around 0Hz) to find dominant walking cadence
for i = 1:length(freqPeaks)-1
    if peaks(i) > maxPeak && f(freqPeaks(i)) > minFreq && f(freqPeaks(i)) < maxFreq
        maxPeak = peaks(i);
        maxFreqIndex = freqPeaks(i);
    end
end
hold on
plot (f(maxFreqIndex), maxPeak, 'o')

% Set a range of frequencies to take according to the dominant cadence
freqWidth = 0.25;
freqRange = [f(maxFreqIndex) - f(maxFreqIndex)*freqWidth, f(maxFreqIndex) + f(maxFreqIndex)*freqWidth];
if freqRange(1)<minFreq
    freqRange(1) = minFreq;
end

if freqRange(2) > maxFreq
    freqRange(2) = maxFreq;
end

% Plot the range of frequency allowed in the bandpass filter
hold on;
ylimits = double(ylim);
a = area([freqRange(1), freqRange(1), freqRange(2), freqRange(2)], [0, ylimits(2), ylimits(2), 0]);
a.FaceAlpha = 0.2;

% 6th order bandpass filter for the determined frequency range
[b,a] = butter(6, freqRange/(fs/2), 'bandpass');
figure();
G=tf(b,a,1/fs); 
freqz(b,a,[],fs);
title("6th Order Bandpass Filter Bode Plot")
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
minPeakMeasured = 0.45*sum(accelTopQuarter)/length(accelTopQuarter);

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

% Plot the detected steps on the filtered time signal if detected
if exist('steps','var')
    fprintf("The pedometer has detected %d steps\n", length(steps));
    hold on
    plot(t(timeStep), steps, 'o');
    yline(minPeakMeasured,'-','Threshold');
    xL=xlim;
    yL=ylim;
    str = "The pedometer has detected " +num2str(length(steps)) +" steps";
    text(0.99*xL(2),0.99*yL(2),str,'HorizontalAlignment','right','VerticalAlignment','top')
else
    fprintf("The pedometer did not detect any steps\n")
end