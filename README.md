# Convolutional PCA
A PCA for multiple time series. It looks for a small set of series, whose filtered versions can explain most of the variances of observations. This kind of PCA is known as PCA in the frequency domain, spectral PCA, or dynamic PCA. Existing work solves this problem exclusively in the frequency domain. I reformulated it in the time or z-domain, thus  completely avoiding issues like aliasing inherent to frequency domain processing, and making such PCA more explainable with reasonable decomposition constraints like sparsity.    

Please check the enclosed pdf file for technical details. You need Matlab or Octave to run these demos.
### Introduction to demos
##### demo_verification_of_theory
It verifies two most import properties of our PCA: 

1) Conservation of variances: (total variances) = (variances of principal components) + (variances of residual errors), regardless of the order L of filter W(z). 
2) W(z) is approximately (due to finite L) paraunitary.      

##### demo_ConvPCA_adaptive
A toy online signal detection demo showing convergence curve and evolution of filter taps. An LMS like stochastic gradient descent algorithm is used to update the filter coefficients. As shown here, it is possible to recover the original mixing coefficients if they are sparse enough. 

![alt text](https://github.com/lixilinx/ConvPCA/blob/master/toy.png)

##### demo_delay_estimation
A Time Difference of Arrival (TDOA) estimation demo. The data are two microphone recordings of a single speech source in a highly reverberant conference room. Convolutional PCA can be used to estimate the two source-to-mic relative responses. By imposing strong sparsity constraints on W(z), delay gap between the two peaks in these two responses can show the TDOA. Comparison results with generalized cross-correlation (GCC), old but still widely used today, is shown as below (the true TDOA is about 12 sampling periods at 16 KHz sampling rate). The four estimated delays in the order of increasing sample sizes are (7, 11, 12, 11) , (5, 0, 12, 13), and (8, 12, 12, 12) for cross correlation, GCC-PHAT, and convolutional PCA, respectively.

![alt text](https://github.com/lixilinx/ConvPCA/blob/master/tde.png)

##### demo_ecg
A biosignal analysis demo. The electrocardiogram (ECG) data of a pregnant woman is downloaded from [here](http://homes.esat.kuleuven.be/~smc/daisy/). Please check the original source for further details. Analysis results of different PCAs with different number of principal components are shown as below. Convolutional PCA tends to preserve more variances, and thus information, with the same number of components. It still keeps fetus ECG well with only two principal components. Aliasing is inherent to frequency domain processing. It introduces hardly explainable noises, and could destroy the conservation of variances for PCA done in the frequency domain.      

![alt text](https://github.com/lixilinx/ConvPCA/blob/master/ecg.png)
