The below files are the results of running simulations:

1) variances_minus_18_to_4.mat

This is a computed variance of test statistics generated from -18 to -4 dB SNRs

2) num_misdetect_minus_18_to_4.mat

This is the number of times (out of 100) that the algorithm
failed to detect an LTE signal when there was one present

3) Probability Misdetection Curve

Shows probability of misdetection vs SNR.


How to obtain results from scratch:

1) Run Variance_Calculator.m as is

Code takes ~50 mins

This will give you result variances_minus_18_to_4.mat

2) Run Curve_Generation2.m as is

Code takes ~1.5 hours

This will give you num_misdetect_minus_18_to_4.mat and Probability_Misdetection_Curve.fig