% Load the LTE waveform with desired configuration
load('lte.mat');

% Print the configuration to the window for viewing
rmcconfig

% Create a time axis based on the LTE sampling rate
% Note that the sampling rate is in Hz, however the signal
% is not necessarilly one second of data so the number
% of evenly spaced points must be truncated
t = linspace(0,.01,(rmcconfig.SamplingRate)/(rmcconfig.SamplingRate/length(rmcwaveform)));

% Plot the real version of the waveform to see the
% signal in the time domain
plot(t,real(rmcwaveform));
title('RMC Waveform in Time');
xlabel('Time (s)');
ylabel('Signal Strength');

% Create new plot for Resource Block Grid
figure;

% From LTE_Downlink.m, we can visualize the resource grid
% of the signal
% Populated resource grid, returned as a numeric 3-D array of resource elements
% for several subframes across all configured antenna ports

% This is the total grid of all 100 RB's for 1 frame
% 12 subcarriers * 100 RB's = 1200 subcarriers
% 7 symbols/slot (normal CP) * 20 slots/frame = 140 symbols
mesh(abs(rmcgrid));
view(2);
title('Resource Element Grid');
ylabel('Subcarriers');
xlabel('Symbols');

% New figure to show a subset of the resource grid in 3 dimensions
figure;
% Grab 12 sub carriers and 2 slots (14 symbols)
rmcgrid_subset = rmcgrid(1188:1200,1:14);
% Mesh as before, using new colors and 3d view
mesh(abs(rmcgrid_subset));
mycolors = [0 0 0; 1 1 0; 1 0 0];
colormap(mycolors);
title('Resource Element Grid Subset');
ylabel('Subcarriers');
xlabel('Symbols');
zlabel('Magnitude');

% New figure for 2d view of resource grid subset
figure;
mesh(abs(rmcgrid_subset));
mycolors = [0 0 0; 1 1 0; 1 0 0];
colormap(mycolors);
% Change view to 2d (viewing from above)
view(2);
title('Resource Element Grid Subset');
ylabel('Subcarriers');
xlabel('Symbols');

% Calculations to plot frequency spectrum of LTE signal

%%Time specifications:
Fs = rmcconfig.SamplingRate;                      % samples per second
dt = 1/Fs;                     % seconds per sample
StopTime = length(rmcwaveform)/rmcconfig.SamplingRate;                  % seconds
t = (0:dt:StopTime-dt)';
N = size(t,1);

% Fourier Transform to get Frequency Components
X = fftshift(fft(rmcwaveform));
   
%%Frequency specifications:
dF = Fs/N;                      % hertz
f = -Fs/2:dF:Fs/2-dF;           % hertz
%%Plot the spectrum:
figure;
plot(f,abs(X)/N);
xlabel('Frequency (in hertz)');
title('Magnitude Response');


