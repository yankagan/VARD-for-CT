# Large-scale image reconstruction for transmission tomography with automatic regularization


This is a demo of the Variational Automatic Relevance Determination algorithm for image reconstruction in transmission tomography under Poisson noise, published in \emph{Y. Kaganovsky, S. Han, S. Degirmenci, D. G. Politte, D. J. Brady, J. A. O'Sullivan and L. Carin, Alternating Minimization Algorithm with Automatic Relevance Determination for Transmission Tomography under Poisson Noise,  SIAM Journal on Imaging Sciences, 8(3), 2015.} We consider x-ray computed tomography (CT), but the algorithm applies more generally, e.g., it can be used for electron tomography.  VARD stands for Variational Automatic Relevance Determination and it is an extension of Automatic Relevance Determination (ARD) [2,3] to Poisson noise models with Beer's law. 
 
The key advantage of ARD or VARD is the lack of any tuning parameters. Typically, in penalized likelihood algorithms, one has to choose values for tuning parameters that control the balance between data-fit and the penalty due to prior knowledge (common choices are the l1 and TV penalties). One is then forced to perform many trials for different values of parameters and to choose the values by trial-and-error.
 
VARD uses a prior that promotes sparsity in the pixel/voxel differences domain and the weight of the prior relative to the data-fit is determined automatically from the data. VARD saves the agonizing pain of choosing values for the tuning parameters; it can also save a lot of time spend on this task.
 
In this demo we compare VARD to the maximum likelihood estimator in [4] and also to a penalized likelihood estimator [1].  All the algorithms in the demo use separable surrogates which reduce the algorithms to parallel line searches, making them feasible for large scale problems such as x-ray CT.
The algorithms can be sped up considerably using various techniques (e.g., using a better initialization, ordered subsets) but these are not considered here.
