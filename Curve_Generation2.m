load('lte.mat') %load LTE signal
load('variances_minus_18_to_4.mat') % Load the variances computed with the variance calculator

P_FA = 0.01;
tests_per_point = 100;
number_of_points = 8; % according to the number of SNR's
num_misdetect = zeros(1, number_of_points); % each index will correspond to 100 tests for a given SNR

point_number = 0;

for SNR_dB = -18 : 2 : -4 % loop through SNRS
    tic    
    point_number = point_number + 1
    decision_threshold = sqrt(-1*variances(point_number)*log(P_FA));
    
    for k = 1 : tests_per_point % compute the test statistic 100 times per point
       k 
       %Build signal with noise
       [noisey_rmcwaveform, noise, noise_variance] = add_awgn_noise(transpose(rmcwaveform), SNR_dB, (k+1000)*-1*SNR_dB); % change the seed
       noisey_rmcwaveform = transpose(noisey_rmcwaveform);
       
       %Generate test statistic

       C = generate_test_statistic(noisey_rmcwaveform);
       
       %Compare to threshold
       
       if C < decision_threshold
           num_misdetect(point_number) = num_misdetect(point_number) + 1; %Accumulate count of misdetect
       end
        
    end
    toc
end

SNR_dB = -18 : 2 : -4;
plot(SNR_dB, num_misdetect/100)
ylim([10^-2 1])
set(gca, 'YScale', 'log')
xlabel 'SNR dB'
ylabel 'Probability of Misdetection'
title 'For P_F_A = 0.01: Probability of Misdetection vs. SNR'
grid on

%%
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

function test_statistic = generate_test_statistic(signal)
    symbol_size1 = 2208;
    symbol_size2 = 2192;
       r=zeros(1,length(signal));

        for n=1:length(signal) - (15360+symbol_size2-1) % 15360 = 2208 + 6 * 2192
            if(mod(n,1000) == 0)
                n;  % to see progress
            end
            
            xl = signal(n:n+symbol_size2-1); % a vector of symbol_size2 length
            xm = signal(n+15360 : n+15360+symbol_size2-1); % a vector of symbol_size2 length that is 7 symbols away
            
            r(n) = (1/symbol_size2) * sum(xl .* conj(xm)); % calculate cross correlation

        end
        % This marks the odd slots
        idx1 = linspace(1,length(signal),20+1); % from plots, we expect 20 pairs of alike pilot tone symbols
        idx1=idx1(1:20);
        idx1=idx1+15360; % 2208 + 6 * 2192 =  Distance between two symbols that have identical pilot tone information
        % This marks the even slots
        idx2 = linspace(1,length(signal),20+1); % from plots, we expect 20 pairs of alike pilot tone symbols
        idx2=idx2(1:20);
        idx2=idx2+8784;  % 2208 + 3 * 2192 = Distance between two symbols that both have pilot tone information
        
        idx = sort([round(idx1) round(idx2)]);
        % For now, assume that these are the best spots
        r(idx);
        C = sum(r(idx)) / 37; % where 37 is the number of r terms
        C = abs(C);
        test_statistic = C;
end
