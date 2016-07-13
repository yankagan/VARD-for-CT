# Large-scale image reconstruction for transmission tomography with automatic regularization


This is a MATLAB demo of the Variational Automatic Relevance Determination (VARD) algorithm for image reconstruction in transmission tomography under Poisson noise that was published in [1]. VARD is an extension of Automatic Relevance Determination (ARD) [2,3] to Poisson noise models with Beer's law used in computed tomography. It also scales significantly better than prior ARD algorithms. 
 
The key advantage of ARD or VARD is the lack of any tuning parameters. Typically, in iterative image reconstruction algorithms, one has to choose values for tuning parameters that control the balance between the fit of the model to the data and the penalty that promotes solutions with desired properties (common choices are the Lasso and Total Variation penalties). One is then forced to perform many trials for different values of tuning parameters in order to get a good performance. 
 
VARD uses a penalty that promotes sparsity in the pixel/voxel differences domain and its weight relative to the data-fit is determined automatically from the data. VARD saves the agonizing pain of choosing values for the tuning parameters; it can also save a lot of time spend on this task.
 
In this demo we compare VARD to the prior image reconstruction algorithms in [4] and [3].  All the algorithms in the demo use separable surrogates which reduce the algorithms to parallel line searches, making them feasible for large scale problems such as x-ray CT (they scale linearly with the number of pixels instead of cubically as in [3]). 

The VARD algorithm can be sped up considerably using various techniques (e.g., using a better initialization, ordered subsets) but these are not considered in the current version.

References
----------
1. Y. Kaganovsky, S. Han, S. Degirmenci, D. G. Politte, D. J. Brady, J. A. O'Sullivan and L. Carin, Alternating Minimization Algorithm with Automatic Relevance Determination for Transmission Tomography under Poisson Noise,  SIAM Journal on Imaging Sciences, 8(3), 2015.  
2. R. M. Neal, Bayesian Learning for Neural Networks, PhD Thesis, 1995. 
3. M. E. Tipping, Sparse Bayesian Learning and the Relevance Vector Machine, Journal of Machine Learning Research, 2001. 
4. J. A. O'Sullivan and J. Benac, Alternating Minimization Algorithm for Transmission Tomography, IEEE Trans. on Medical Imaging, 2007.

How to Use the Code
--------------------
The code is provided for EDUCATIONAL purposes only!
(1) For a demo that includes only the VARD algorithm, please run the file: Demo_VARD.m 
(2) For a demo that includes a comparison between VARD and other algorithms, please run the file: Demo_all.m

