load('lte.mat') %load LTE signal
%Calculate the test statistic variance for various SNRs
%Compute the variance of the test statistic for the desired SNR

symbol_size1 = 2208;
symbol_size2 = 2192;

variances = zeros(1,8);
k = 1 % the kth SNR

for SNR_dB = -18 : 2 : -4

    tic   
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
                xl = noisey_rmcwaveform(n:n+symbol_size2-1); % a vector of symbol_size2 length
                xm = noisey_rmcwaveform(n+15360 : n+15360+symbol_size2-1); % a vector of symbol_size2 length that is 7 symbols away
                
                r(n) = (1/symbol_size2) * sum(xl .* conj(xm)); % calculate cross correlation

            end
            % This marks the odd slots
            idx1 = linspace(1,length(noisey_rmcwaveform),20+1);
            idx1=idx1(1:20);
            idx1=idx1+15360;  % 2208 + 6 * 2192 =  Distance between two symbols that have identical pilot tone information
            % This marks the even slots
            idx2 = linspace(1,length(noisey_rmcwaveform),20+1);
            idx2=idx2(1:20);
            idx2=idx2+8784;  % 2208 + 3 * 2192 = Distance between two symbols that both have pilot tone information

            idx = sort([round(idx1) round(idx2)]);
            % For now, assume that these are the best spots
            r(idx);
            C = sum(r(idx)) / 37; % where 37 is the number of r terms
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