%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Yuya HOSODA, Arata KAWAMURA, Youji IIGUNI,                                                   
% "Speech Bandwidth Extension using Data Hiding based on Discrete Hartley Transform-Domain,"  
% Circuits, Syst., Signal Process., (under review).                                           
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Outline 
% The public switching telephone network restricts a speech signal to the narrow bandwidth (NB) 
% of 0.3-3.4 kHz, which results in the quality decline due to missing the upper bandwidth (UB) 
% spectrum of 3.4-7 kHz. This paper proposes a speech bandwidth extension method that reconstructs 
% the missing UB spectrum using side information. The sender side obtains a UB spectral envelope 
% and relative gains between NB and UB excitation signals as side information. Side information 
% is then converted into a binary signal using two codebooks. Using speech steganography based 
% on the discrete Hartley transform (DHT) domain, the proposed method robustly embeds the binary 
% signal into an amplitude spectrum of the NB speech signal in the high-frequency bandwidth of 
% 3.4--4.6 kHz to produce a composite narrow bandwidth (CNB) speech signal. On the receiver side, 
% the missing UB spectrum is reconstructed using side information extracted from the CNB speech signal. 
% Theoretical and simulation analysis shows that side information is retrieved from the CNB speech 
% signal accurately. Subjective listening tests and objective measures also show that the proposed 
% method enhances the quality of the NB speech signal by reconstructing the missing UB spectrum.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Program branch
% CSSP_Main.m
%  |-SetUp.m
%  |-Extract_BinarySignal_LSF.m
%  | -LSF_to_BinarySignal.m
%  |-Extract_LSF.m
%  |-Extract_BinarySignal_LSFandGain.m
%  | -LSFandGain_to_BinarySignal.m
%  |-Extract_BinarySignal_Gain.m
%  | -Gain_to_BinarySignal.m
%  |-Data_Hiding_DFT.m
%  |-Data_Hiding_DHT.m
%  |-Extract_BinarySignal_DFT.m
%  |-Extract_BinarySignal_DHT.m
%  |-Generate_UBspeech_SingleGain.m
%  |-Generate_UBspeech_Gains.m
%   -Generate_WBspeech.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input / Output
% Input  : Original WB speech signal at sampling rate of 16 kHz
% Output : Generated WB speech signal at sampling rate of 16 kHz
%           WB_Sp_DFT_ConvBit 
%           -WB speech using TDDH based on DFT domain with LSF and a relative gain (Conv.) 
%           WB_Sp_DHT_ConvBit
%           -WB speech using TDDH based on DHT domain with LSF and a relative gain
%           WB_Sp_DHT_PropBit
%           -WB speech using TDDH based on DHT domain with LSF and relative gains
%           WB_Sp_DFT_PropBit
%           -WB speech using TDDH based on DFT domain with LSF and relative gains (Prop.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
%%%%%%%%%%%%%%%%%%%%%%%
%   Global variable   %
%%%%%%%%%%%%%%%%%%%%%%%

% Setting global variable
SetUp();
    

%%%%%%%%%%%%%%%%%
%   Parameter   %
%%%%%%%%%%%%%%%%%

% Smapling rate for NB speech signal
Fs_NBspeech = 8000;

% The number of frame samples (20 ms)
Frame_length = Fs_NBspeech*0.020;

% Strength for Data Hiding    
Strength_DataHiding = 0.01;
    
% Length of binary signal
BinarySignal_length = 12;

% Length of binary signal for LSF(Prop.)
BinarySignal_length_LSF = 4;

% Length of binary signal for relative gain(Prop.)
BinarySignal_length_Gain = 8;

% Length of hidden vector to be embedded into magnitude in bandwidth of 3.4-4kHz(Conv)
HiddenVector_length_Conv = 12;

% Length of hidden vector to be embedded into amplitude in bandwidth of 3.4-4.6kHz(Prop)
HiddenVector_length_Prop = 24;

%%%%%%%%%%%%%%%%%
%   Load file   %
%%%%%%%%%%%%%%%%%

% Load input speech file (WB speech signal in bandwidth of 0.3-7 kHz)
input_WBspeech_name = 'input_speech.wav';
[Input_WBspeech, Fs_WBspeech] = audioread(input_WBspeech_name);


%%%%%%%%%%%%%%%%%%%
%   Filter Pass   %
%%%%%%%%%%%%%%%%%%%

% Generate NB speech signal in bandwidth of 0.3-3.4 kHz
Input_NBspeech = lowpass(Input_WBspeech, 3400, Fs_WBspeech);

% Generate UB speech signal in bandwidth of 3.4-7 kHz
Input_UBspeech = Input_WBspeech-Input_NBspeech;

% Frequency shift (3.4-7 kHz -> 0-3.4 kHz)
time = 1:length(Input_UBspeech);
Input_UBspeech_shift = 2*Input_UBspeech.*cos(2*pi*3400*time'/Fs_WBspeech);
Input_UBspeech_shift = lowpass(Input_UBspeech_shift, Fs_NBspeech/2, Fs_WBspeech);

% Down-sampling(16 kHz -> 8 kHz)
Input_NBspeech_8kHz = downsample(Input_NBspeech, 2);
Input_UBspeech_8kHz = downsample(Input_UBspeech_shift, 2);

Input_NBspeech_8kHz = Input_NBspeech_8kHz(1:floor(length(Input_NBspeech_8kHz)/Frame_length)*Frame_length);
Input_UBspeech_8kHz = Input_UBspeech_8kHz(1:floor(length(Input_UBspeech_8kHz)/Frame_length)*Frame_length);


%%%%%%%%%%%%%%%%%%%%%%%
%   Feature Extract   %
%%%%%%%%%%%%%%%%%%%%%%%

% Obtain a binary signal for LSF(Prop.)
BinarySignal_LSF_UBspeech = Extract_BinarySignal_LSF(Input_UBspeech_8kHz, Frame_length, BinarySignal_length_LSF);

% Extract LSF from UB speech signal(Conv.)
LSF_UBspeech = Extract_LSF(Input_UBspeech_8kHz, Frame_length);


%%%%%%%%%%%%%%%%%%%%%
%   Gain Encoding   %
%%%%%%%%%%%%%%%%%%%%%

% Obtain a binary signal for LSF and relative gain(Conv.)
BinarySignal_LSFGain_Conv = Extract_BinarySignal_LSFandGain(Input_UBspeech_8kHz, Input_NBspeech_8kHz, LSF_UBspeech, Frame_length, BinarySignal_length);

% Obtain a binary signal for relative gain(Prop.)
BinarySignal_Gain = Extract_BinarySignal_Gain(Input_UBspeech_8kHz, Input_NBspeech_8kHz, BinarySignal_LSF_UBspeech, Frame_length,  BinarySignal_length_Gain);

% Incorpolate binary signals(Prop.)
BinarySignal_LSFGain_Prop = [BinarySignal_LSF_UBspeech, BinarySignal_Gain];


%%%%%%%%%%%%%%%%%%%
%   Data Hiding   %
%%%%%%%%%%%%%%%%%%%

% CNB speech using TDDH based on DHT domain with relative gains(Prop.)
CNBspeech_DHT_Gains      = Data_Hiding_DHT(Input_NBspeech_8kHz, BinarySignal_LSFGain_Prop, Frame_length, Strength_DataHiding, HiddenVector_length_Prop);

% CNB speech using TDDH based on DHT domain with a relative gain 
CNBspeech_DHT_SingleGain = Data_Hiding_DHT(Input_NBspeech_8kHz, BinarySignal_LSFGain_Conv, Frame_length, Strength_DataHiding, HiddenVector_length_Prop);

% CNB speech using TDDH based on DFT domain with relative gains 
CNBspeech_DFT_Gains      = Data_Hiding_DFT(Input_NBspeech_8kHz, BinarySignal_LSFGain_Prop, Frame_length, Strength_DataHiding, HiddenVector_length_Conv);

% CNB speech using TDDH based on DFT domain with a relative gain(Conv.)
CNBspeech_DFT_SingleGain = Data_Hiding_DFT(Input_NBspeech_8kHz, BinarySignal_LSFGain_Conv, Frame_length, Strength_DataHiding, HiddenVector_length_Conv);


%%%%%%%%%%%%%%%%%%%%%%%
%   Coding/Encoding   %
%%%%%%%%%%%%%%%%%%%%%%%

% Default : Clean environment without speech codecs

CNBspeech_DHT_Gains_Receive      = CNBspeech_DHT_Gains;
CNBspeech_DHT_SingleGain_Receive = CNBspeech_DHT_SingleGain;
CNBspeech_DFT_SingleGain_Receive = CNBspeech_DFT_SingleGain;
CNBspeech_DFT_Gains_Receive      = CNBspeech_DFT_Gains;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Binary signal Extraction  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract binary signal from CNB speech
[BinarySignal_DHT_Gains,      NBspeech_DHT_Gains]      = Extract_BinarySignal_DHT(CNBspeech_DHT_Gains_Receive,  Frame_length, HiddenVector_length_Prop, BinarySignal_length);
[BinarySignal_DHT_SingleGain, NBspeech_DHT_SingleGain] = Extract_BinarySignal_DHT(CNBspeech_DHT_SingleGain_Receive,  Frame_length, HiddenVector_length_Prop, BinarySignal_length);
[BinarySignal_DFT_SingleGain, NBspeech_DFT_SingleGain] = Extract_BinarySignal_DFT(CNBspeech_DFT_SingleGain_Receive, Frame_length, HiddenVector_length_Conv, BinarySignal_length, Strength_DataHiding);
[BinarySignal_DFT_Gains,      NBspeech_DFT_Gains]      = Extract_BinarySignal_DFT(CNBspeech_DFT_Gains_Receive,  Frame_length, HiddenVector_length_Conv, BinarySignal_length, Strength_DataHiding);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   UB speech Generation   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Conventional UB speech generation using side information
UBspeech_DFT_SingleGain = Generate_UBspeech_SingleGain(BinarySignal_DFT_SingleGain, NBspeech_DFT_SingleGain', Frame_length);
UBspeech_DHT_SingleGain = Generate_UBspeech_SingleGain(BinarySignal_DHT_SingleGain, NBspeech_DHT_SingleGain', Frame_length);

% Proposed UB speech generation using side information
UBspeech_DHT_Gains = Generate_UBspeech_Gains(BinarySignal_DHT_Gains, NBspeech_DHT_Gains', Frame_length, BinarySignal_length_LSF);
UBspeech_DFT_Gains = Generate_UBspeech_Gains(BinarySignal_DFT_Gains, NBspeech_DFT_Gains', Frame_length, BinarySignal_length_LSF);


%%%%%%%%%%%%%%
%   Output   %
%%%%%%%%%%%%%%

% WB speech generation
WBspeech_DFT_SingleGain = Generate_WBspeech(NBspeech_DFT_SingleGain', UBspeech_DFT_SingleGain);
WBspeech_DHT_SingleGain = Generate_WBspeech(NBspeech_DHT_SingleGain', UBspeech_DHT_SingleGain);
WBspeech_DHT_Gains = Generate_WBspeech(NBspeech_DHT_Gains', UBspeech_DHT_Gains);
WBspeech_DFT_Gains = Generate_WBspeech(NBspeech_DFT_Gains', UBspeech_DFT_Gains);

audiowrite('WBspeech_DFT_SingleGain.wav', WBspeech_DFT_SingleGain,  16000);
audiowrite('WBspeech_DHT_SingleGain.wav', WBspeech_DHT_SingleGain,  16000);
audiowrite('WBspeech_DHT_Gains.wav', WBspeech_DHT_Gains,  16000);
audiowrite('WBspeech_DFT_Gains.wav', WBspeech_DFT_Gains,  16000);






