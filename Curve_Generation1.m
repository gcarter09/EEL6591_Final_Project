%Then run this section to get results for various probabilities of false
%alarm
%C_variance_minus_30dB = 3.8470e-10;
C_variance_minus_25dB = 6.6890e-11; %This is a function of the SNR (original 3.8943e-11)
%C_variance_minus_24dB = 4.1696e-11;
C_variance_minus_23dB = 2.3425e-11;
C_variance_minus_20dB = 4.2770e-12; 
%C_variance_minus_15dB = 6.4222e-13; 

SNR_dB = -23;
number_of_points = 18; % according to the number of P_FA's
tests_per_point = 100;


num_detect = zeros(1, number_of_points);
num_FA = zeros(1, number_of_points);


point_number = 0;
tic
for P_FA = 0.05 : 0.05: 0.9
   point_number = point_number + 1
   
   %set threshold
   decision_threshold = sqrt(-1*C_variance_minus_23dB*log(P_FA));
   
   for k = 1 : tests_per_point
       
       %Build signal with noise
       [noisey_rmcwaveform, noise, noise_variance] = add_awgn_noise(transpose(rmcwaveform), SNR_dB, k+100000); % change the seed
       noisey_rmcwaveform = transpose(noisey_rmcwaveform);
       
       %Generate test statistic

       C = generate_test_statistic(noisey_rmcwaveform);
       
       %Compare to threshold
       
       if C > decision_threshold
           num_detect(point_number) = num_detect(point_number) + 1; %Accumulate count of detect
       end
       
       
       %Build just noise
       
       %[noisey_rmcwaveform, noise, noise_variance] = add_awgn_noise(transpose(rmcwaveform), SNR_dB, k+100); % k is the seed.  Dont want same random seed
       %noise = transpose(noise);
       
       %Generate test statistic
       
       %C = generate_test_statistic(noise);
       
       %Accumulate count of false alarm
       %if C > decision_threshold
       %    num_FA(point_number) = num_FA(point_number) + 1; %Accumulate count of false alarm
       %end
   
   end
    
end
toc
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
            xl = signal(n:n+symbol_size2-1);
            xm = signal(n+15360 : n+15360+symbol_size2-1);
            
            r(n) = (1/symbol_size2) * sum(xl .* conj(xm));  % i belive this is the right way to do it

        end
        % This marks the odd slots
        idx1 = linspace(1,length(signal),20+1);
        idx1=idx1(1:20);
        idx1=idx1+15360; % 2208 + 6 * 2192
        % This marks the even slots
        idx2 = linspace(1,length(signal),20+1);
        idx2=idx2(1:20);
        idx2=idx2+8784;  % 2208 + 3 * 2192
        
        idx = sort([round(idx1) round(idx2)]);
        % For now, assume that these are the best spots
        r(idx);
        C = sum(r(idx)) / 37; % where for is the number of r terms
        C = abs(C);
        test_statistic = C;
end