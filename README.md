## Semi-Blind Denoising of EEG signals

This repository contains code for my remote project "Semi-Blind Denoising of EEG Signals", and I developed two innovative algorithms to clean extremely noisy EEG signals recorded during Galvanic Vestibular Stimulation in an MRI machine. In the first method, I employed denoising source separation, an expectation-maximization approach, using adaptive noise cancellation as the expectation step to extract noise components. These components were fed into another adaptive noise cancellation process to further denoise each EEG channel. The second method involved semi-blind Independent Component Analysis (ICA), where I regularized the loss function of the Infomax algorithm to ensure that the target noise component was as correlated as possible with the stimulation signal. The methods are applied after removing first principal component and after applying EEGLAB's fMRIB to remove the gradient artifact since the signals are extremely noisy. The results for the first algorithm are provided in Plots, the second method performs better on non-stationary stimulations and the results will be provided later.

## Sample ANC-DSS denoising

![alt text](https://github.com/arminnv/Semi-blind-denoising-of-EEG-signals/blob/a2a32592d3f6dc735d367e926940298512b9295b/ANC_DSS/Plots/Exp1-flag11.png?raw=true)
