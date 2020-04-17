% This code calculates the test statistic variance for various SNRs

% Load the generated LTE communications signal from LTE_Downlink.m
% In order to regenerate this signal, the LTE Toolbox is needed
% This loads the RMC (Reference Measurement Channel) Waveform (along with the grid and configuration)
load('lte.mat')

symbol_size1 = 2208;
symbol_size2 = 2192;

% 8 different SNR's will be generated as can be seen below
variances = zeros(1,8);
% Loop counter for saving the calculated variance
k = 1

% This for loop iterates through various SNR values
% going from -18 to -4 by 2's gives 8 values
for SNR_dB = -18 : 2 : -4

    tic % Used for timing the loop  
    c = zeros(1,50);

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

            idx = sort([round(idx1) round(idx2)]);
            % For now, assume that these are the best spots
            r(idx);
            C = sum(r(idx)) / 37; % where for is the number of r terms
            c(i) = abs(C);
            %next compute c varaince

        end

    variances(k) = var(c);
    
    k = k + 1
    toc

end

% Add noise to LTE signal
%https://www.gaussianwaves.com/2015/06/how-to-generate-awgn-noise-in-matlaboctave-without-using-in-built-awgn-function/
%Authored by Mathuranathan Viswanathan
%How to generate AWGN noise in Matlab/Octave by Mathuranathan Viswanathan
%is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
%You must credit the author in your work if you remix, tweak, and build upon the work below

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