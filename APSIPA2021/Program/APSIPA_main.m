%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Yuya HOSODA, Arata KAWAMURA, Youji IIGUNI,                                                   
% "Pitch Estimation Algorithm for Narrowband Speech Signal using Phase Differences between Harmonics,"  
% in Proc. APSIPA2021,2021 (to be accepted).                                           
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Outline 
% This paper proposes a pitch estimation algorithm for a narrowband speech signal using phase difference 
% between harmonics. A narrowband speech signal has an incomplete harmonic structure due to 
% bandwidth limitation, which degrades the pitch estimation accuracy. In this paper, we focus on 
% the fact that phase differences between harmonics are constant by approximating a speech signal 
% based on a sinusoidal model. Since the narrowband speech signal has a partial harmonic structure, 
% the proposed method selects the pitch such that the phase differences between harmonics are constant 
% on the narrow bandwidth. Also, we take a frame additive average method for the phase differences 
% between harmonics to improve the robustness against noise. Experimental results show that 
% the proposed method estimates the pitch from the narrowband speech signal under noisy environments 
% more accurately than the traditional methods.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Program branch
% EUSIPCO_Main.m
% |-Proposed_pitch_estimation.m (Proposed method)
% |-YIN.m (YIN algorithm)
% |-SRH.m (SRH algorithm)
% |-Evaluation.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input / Output
% Input  : a speech signal at a sampling rate of 16 kHz
%          a noise signal at a sampling rate of 16 kHz
% Output : Pitch [Hz]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Relatied works
%  [1] S. Gonzalez and M. Brookes,
%  ``PEFAC-A pitch estimation algorithm robust to high levels of noise,"
%  IEEE/ACM Trans. Audio, Speech, Language Process., vol.22, no.2, pp.518-530, 2014.
%  [2] A. e Cheveign\'{e} and H. Kawahara,
%  ``YIN, a fundamental frequency estimator for speech and music,"
%  J. Acoust. Soc. Am., vol.111, no.4, pp.1917-1930, 2002.
%  [3] T. Drugman and A. Alwan,
%  ``Joint robust voicing detection and pitch estimation based on residual harmonics,"
%  in Proc. INTERSPEECH, 2011, pp.1973-1976.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
%   Parameter   %
%%%%%%%%%%%%%%%%%

% Input Signal-to-Noise ratio (SNR) level [dB]
Input_SNR = 0;


%%%%%%%%%%%%%%%%%%
%   Load files   %
%%%%%%%%%%%%%%%%%%

% Load input speech signal at a sampling rate of 16 kHz
[Input, ~] = audioread('Input_speech.wav');

% Original pitch of the input speech signal (analyzed frame 32 ms)
Original_Pitch = load('Original_pitch.txt');

% In this program, we assume that Voice Active Detection (VAD) has been known
% We set a voice active frame (VAD=1) such as the original pitch (> 50 Hz) exists.
Original_VAD = Original_Pitch>50;


%%%%%%%%%%%%%%%%%%%%%%%
%   Noise addiction   %
%%%%%%%%%%%%%%%%%%%%%%%

% Load a noise signal at a sampling rate of 16 kHz
[Noise, ~] = audioread('Noise.wav');

% Calculate a noise power with input SNR level
Noise_Power = norm(Input)/norm(Noise)*10^(-Input_SNR/20);

% Add the noise signal to the input speech signal with input SNR level
Input = Input+Noise_Power.*Noise;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Public switching telephone network   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filtering processes
[High_pass_1st_a,High_pass_1st_b] = butter(2,  300/8000,'high');
High_pass_2nd = fir1(800, 0.0325, 'high', chebwin(801, 80));
[Low_pass_a,Low_pass_b] = butter(10,  3300/8000);

Input_NB_UB_1st  = filter(High_pass_1st_a, High_pass_1st_b, Input);
Input_NB_UB     = filter(High_pass_2nd, 1, Input_NB_UB_1st);
Input_NB        = filter(Low_pass_a, Low_pass_b, Input_NB_UB);

% Down-sampling at 8 kHz
Input_NB = downsample(Input_NB, 2);
Fs = 8000;

% Normalization
Input_NB = Input_NB./max(abs(Input_NB));

% G.711 encoding & decoding
Input_NB = mu2lin(lin2mu(Input_NB)); 


%%%%%%%%%%%%%%%%%%%%%%%%
%   Pitch estimation   %
%%%%%%%%%%%%%%%%%%%%%%%%

% PEFAC algorithm [1]
% Frame length 96ms / Flame shift 10ms
% Defalt pitch estimation algorithm (installed in Matlab R2020a)
PEF_Pitch = pitch(Input_NB, Fs, 'Method', 'PEF', 'Range', [50, 400], 'WindowLength', round(0.096*Fs),'OverlapLength', round(0.086*Fs), "MedianFilterLength", 1);
PEF_Pitch = [zeros(2, 1); PEF_Pitch; zeros(length(Original_VAD)-length(PEF_Pitch)-2, 1)];
PEF_Pitch = PEF_Pitch.*Original_VAD;

% YIN algorithm [2]
% Frame length 40ms / Flame shift 10ms
YIN_Pitch = YIN(Input_NB, Original_VAD, Fs, 0.040*Fs, 0.010*Fs);

% SRH algorithm [3]
% Frame length 96ms / Flame shift 10ms 
SRH_Pitch = SRH(Input_NB, Original_VAD(3:end-2), Fs, Fs*0.096, Fs*0.010);
SRH_Pitch = [zeros(2, 1); SRH_Pitch; zeros(length(Original_VAD)-length(SRH_Pitch)-2, 1)];

% Proposed pitch estimation algorithm
Prop_Pitch = Proposed_pitch_estimation(Input_NB, Original_VAD, Fs, 0.040*Fs, 0.010*Fs);


%%%%%%%%%%%%%%%%%%
%   Evaluation   %
%%%%%%%%%%%%%%%%%%

Evaluation_PEF  = Evaluation(Original_Pitch, PEF_Pitch);
Evaluation_SRH  = Evaluation(Original_Pitch, SRH_Pitch);
Evaluation_YIN  = Evaluation(Original_Pitch, YIN_Pitch);
Evaluation_Prop = Evaluation(Original_Pitch, Prop_Pitch);



