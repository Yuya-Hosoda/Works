# Speech Bandwidth Extension using Data Hiding based on Discrete Hartley Transform Domain

# Description:
This work performs a speech bandwidth extension method using data hiding based on discrete Hartley transform domain.

# Abstract
The public switching telephone network restricts a speech signal to the narrow bandwidth (NB) of 0.3–3.4 kHz, which results in the quality decline due to the missing upper bandwidth (UB) spectrum of 3.4–7 kHz. This paper proposes a speech bandwidth extension method that reconstructs the missing UB spectrum using side information. The sender side obtains a UB spectral envelope and relative gains between NB and UB excitation signals as side information. Side information is then converted into a binary signal using two codebooks. Using speech steganography based on the discrete Hartley transform (DHT) domain, the proposed method robustly embeds the binary signal into an amplitude spectrum of the NB speech signal in the high-frequency bandwidth of 3.4–4.6 kHz to produce a composite narrow bandwidth (CNB) speech signal. On the receiver side, the missing UB spectrum is reconstructed using side information extracted from the CNB speech signal. Theoretical and simulation analysis shows that side information is retrieved from the CNB speech signal accurately. Subjective listening tests and objective measures also show that the proposed method enhances the quality of the NB speech signal by reconstructing the missing UB spectrum.

# Contents
- Program : This folder contains scripts to perform speech bandwidth extension methods.


____________________________________________________________________________
All the scripts are successfully tested on MATLAB 2020a
____________________________________________________________________________
