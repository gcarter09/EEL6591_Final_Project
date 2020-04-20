For the Prob{detect} vs. Prob{false alarm} graph:

Results:

1) SNR_minus_20_db_results.mat
2) SNR_minus_23_db_results.mat
3) SNR_minus_25_db_results.mat

The number of times (out of 100) that we correctly identified an LTE signal when there was indeed one present

4) Probability Detection Curve

Shows Prob{detect} vs. Prob{false alarm}


How to obtain results from scratch:

1) Run Variance_Calculator.m
	
Run this for the desired SNRs of -20, -23, and -25

2) Run Curve_Generation1.m

Specify which SNR and its variance.  Detection data will be generated for that SNR

___________________________________________________________________________________________________________
For the Prob{misdetect} vs. SNR graph:

Results:

1) variances_minus_18_to_4.mat

This is a computed variance of test statistics generated from -18 to -4 dB SNRs

2) num_misdetect_minus_18_to_4.mat

This is the number of times (out of 100) that we failed to detect an LTE signal when there was indeed one present

3) Probability Misdetection Curve

Shows probability of misdetection vs SNR.


How to obtain results from scratch:

1) Run Variance_Calculator.m as is

Code takes ~50 mins

This will give you result variances_minus_18_to_4.mat

2) Run Curve_Generation2.m as is

Code takes ~1.5 hours

This will give you num_misdetect_minus_18_to_4.mat and Probability_Misdetection_Curve.fig