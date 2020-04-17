%% Cannot run this without LTE Toolbox
%  The following example shows how to create a 20MHz, QPSK, 3/4 rate
%  waveform corresponding to transmission mode 8 ('Port7-8' transmission
%  scheme) with full allocation 
dataStream = [];     % Define the input user data stream (was [1 0 0 1] )
params = struct();          % Initialize the parameter structure
params.NDLRB = 100;         % 20 MHz bandwidth
params.CellRefP = 1;        
params.PDSCH.PRBSet = (0:params.NDLRB-1)'; % Full allocation
params.PDSCH.TargetCodeRate = 3/4; % The target code rate
params.PDSCH.TxScheme = 'Port0'; % Transmission mode 0
params.PDSCH.NLayers = 1;          % 1 layer transmission
params.PDSCH.Modulation = 'QPSK';  % Modulation scheme
params.PDSCH.NSCID = 0;            % Scrambling identity
params.PDSCH.NTxAnts = 1;          % 1 transmit antennas
params.PDSCH.W = lteCSICodebook(params.PDSCH.NLayers,...
                        params.PDSCH.NTxAnts,0).'; % Precoding matrix

% Now use lteRMCDL to populate other parameter fields
fullParams = lteRMCDL(params);
% Generate the waveform using the full parameter set 'fullParams'
[rmcwaveform, rmcgrid, rmcconfig] = lteRMCDLTool(fullParams,dataStream);
% dlWaveform is the time domain waveform, dlGrid is the resource grid and
% dlParams is the full set of parameters used in the waveform generation.

% Populated resource grid, returned as a numeric 3-D array of resource elements
% for several subframes across all configured antenna ports
mesh(abs(rmcgrid))
view(2)
title('Resource Element Grid')
ylabel('Subcarriers')
xlabel('Symbols')

disp('LTE Downlink Configuration')
rmcconfig
%%
% 100 RBs.  20 MHz Signal has 100 RBs 
% 12 Subcarriers per RB
% So, 100 * 12 = 1200 Subcarriers

% 10 subframes
% 2 Slots per Subframe
% 7 Symbols per Slot
% So, 10 * 2 * 7 = 140 Symbols

% 1 ms per Subframe
% So we have 10 * 1 = 10 ms of time domain signal

% fs = 30,720,000
% number of samples = 307,200 in time domain signal
% So, yes we have 10 ms of time domain signal
% 307,200 / 140 ~= 2194.28 samples per symbol

symbol_size1 = 2208;
symbol_size2 = 2192;

r=zeros(1,length(rmcwaveform));
for n=1:length(rmcwaveform) - (15360+symbol_size2-1) % 15360 = 2208 + 6 * 2192
    if(mod(n,1000) == 0)
        n  % to see progress
    end
    xl = rmcwaveform(n:n+symbol_size2-1);
    xm = rmcwaveform(n+15360 : n+15360+symbol_size2-1);
%    r(n) = (1/symbol_size2) * max(abs(xcorr(xl,xm)));
%    r(n) = abs(xl' * conj(xm));
%    r(n) = abs(transpose(xl) * conj(xm));
%    r(n) = sum(abs(xl .* conj(xm)))/symbol_size2;
     r(n) = (1/symbol_size2) * sum(xl .* conj(xm)); % might be this
    
end
% This marks the odd slots
idx1 = linspace(1,length(rmcwaveform),20+1);
idx1=idx1(1:20);
idx1=idx1+15360; % 2208 + 6 * 2192
% This marks the even slots
idx2 = linspace(1,length(rmcwaveform),20+1);
idx2=idx2(1:20);
idx2=idx2+8784;  % 2208 + 3 * 2192

hold off
plot(abs(r))
hold
stem(idx1, max(abs(r))*ones(1,length(idx1)))
stem(idx2, max(abs(r))*ones(1,length(idx2)))
hold off
%saveas(gcf,'corr_with_simple_multiply.jpg')
%saveas(gcf,'max_corr.jpg')
%save('lte.mat', 'rmcwaveform', 'rmcgrid', 'rmcconfig');

% Combine idx1 and idx2 and sort them to generate a comb to grab the best N
% correlations that fit the known pilot pattern.  Round since some are not
% integers
idx = sort([round(idx1) round(idx2)]);
% For now, assume that these are the best spots
r(idx);
C = sum(r(idx)) / 37 % where for is the number of r terms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%compute decision threshold for a given SNR_dB

SNR_dB = -20;
symbol_size2 = 2192;

%first compute noise variance
L=length(rmcwaveform);
SNR = 10^(SNR_dB/10); %SNR to linear scale
Esym=sum(abs(rmcwaveform).^2)/(L); %Calculate actual symbol energy
N0=Esym/SNR; %Find the noise spectral density
if(isreal(rmcwaveform))
    noiseSigma = sqrt(N0);%Standard deviation for AWGN Noise when x is real
    n = noiseSigma*randn(1,L);%computed noise
else
    noiseSigma=sqrt(N0/2);%Standard deviation for AWGN Noise when x is complex
    n = noiseSigma*(randn(1,L) + 1i*randn(1,L));%computed noise
end

noise_variance = var(n);

%next compute c varaince

C_variance = (noise_variance^2) / (37 * symbol_size2 * (1200^2)); % where 40 is the number of r terms, symbol_size 2 is the length of an OFDM symbol, 1200 is the number of OFDM subcarriers

%compute decision threshold

prob_false_alarm = .5;

decision_threshold = sqrt(-C_variance*log(prob_false_alarm))

%% First run this section to compute the variance
%Compute the variance of the test statistic for the desired SNR
symbol_size1 = 2208;
symbol_size2 = 2192;

SNR_dB = -23; % desired SNR_dB
%prob_false_alarm = .7; % desired probability of a false alarm

c = zeros(1,50);
tic
for i = 1:50
    i
    [noisey_rmcwaveform, noise, noise_variance] = add_awgn_noise(transpose(rmcwaveform), SNR_dB, i+100);

    noisey_rmcwaveform = transpose(noisey_rmcwaveform);
    noise = transpose(noise);

    r=zeros(1,length(noisey_rmcwaveform));

    for n=1:length(noisey_rmcwaveform) - (15360+symbol_size2-1) % 15360 = 2208 + 6 * 2192
        if(mod(n,1000) == 0)
            n;  % to see progress
        end
        xl = noisey_rmcwaveform(n:n+symbol_size2-1);
        xm = noisey_rmcwaveform(n+15360 : n+15360+symbol_size2-1);

    %    r(n) = max(abs(xcorr(xl,xm)));
    %    r(n) = abs(xl' * conj(xm));
    %    r(n) = abs(transpose(xl) * conj(xm));
        r(n) = (1/symbol_size2) * sum(xl .* conj(xm));  % i belive this is the right way to do it

    end
    % This marks the odd slots
    idx1 = linspace(1,length(noisey_rmcwaveform),20+1);
    idx1=idx1(1:20);
    idx1=idx1+15360; % 2208 + 6 * 2192
    % This marks the even slots
    idx2 = linspace(1,length(noisey_rmcwaveform),20+1);
    idx2=idx2(1:20);
    idx2=idx2+8784;  % 2208 + 3 * 2192

    %hold off
    %plot(r)
    %hold
    %stem(idx1,0.1*ones(1,length(idx1)))
    %stem(idx2,0.1*ones(1,length(idx2)))
    %hold off
    %saveas(gcf,'corr_with_simple_multiply.jpg')
    %saveas(gcf,'max_corr.jpg')
    %save('lte.mat', 'rmcwaveform', 'rmcgrid', 'rmcconfig');

    % Combine idx1 and idx2 and sort them to generate a comb to grab the best N
    % correlations that fit the known pilot pattern.  Round since some are not
    % integers
    idx = sort([round(idx1) round(idx2)]);
    % For now, assume that these are the best spots
    r(idx);
    C = sum(r(idx)) / 37; % where for is the number of r terms
    c(i) = abs(C);
    %next compute c varaince
    
end

C_variance = var(c);

%decision_threshold = sqrt(-C_variance*log(prob_false_alarm));
toc
%C_variance = (noise_variance^2) / (37 * symbol_size2 * (1200^2)); % where 40 is the nuber of r terms, symbol_size 2 is the length of an OFDM symbol, 1200 is the number of OFDM subcarriers

%compute decision threshold

%decision_threshold = sqrt(-C_variance*log(prob_false_alarm))

%%
% Add noise to LTE signal
%https://www.gaussianwaves.com/2015/06/how-to-generate-awgn-noise-in-matlaboctave-without-using-in-built-awgn-function/
%Authored by Mathuranathan Viswanathan
%How to generate AWGN noise in Matlab/Octave by Mathuranathan Viswanathan
%is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
%You must credit the author in your work if you remix, tweak, and build upon the work below

%SNR_dB = -20;

%rmcwaveform_noise = add_awgn_noise(rmcwaveform', SNR_dB);

%rmcwaveform_noise = rmcwaveform_noise';

function [noisy_signal, noise, noise_variance]  = add_awgn_noise(x,SNR_dB,seed)
     %y=awgn_noise(x,SNR) adds AWGN noise vector to signal 'x' to generate a
     %resulting signal vector y of specified SNR in dB
     rng(seed);%set the random generator seed to default (for comparison only)
     L=length(x);
     SNR = 10^(SNR_dB/10); %SNR to linear scale
     Esym=sum(abs(x).^2)/(L); %Calculate actual symbol energy
     N0=Esym/SNR; %Find the noise spectral density
     if(isreal(x))
         noiseSigma = sqrt(N0);%Standard deviation for AWGN Noise when x is real
         n = noiseSigma*randn(1,L);%computed noise
     else
         noiseSigma=sqrt(N0/2);%Standard deviation for AWGN Noise when x is complex
         n = noiseSigma*(randn(1,L) + 1i*randn(1,L));%computed noise
     end
     noise_variance = var(n);
     noise = n;
     noisy_signal = x + n; %received signal
end



