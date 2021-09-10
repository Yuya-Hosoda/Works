function [BinarySignal, NBspeech] = Extract_BinarySignal_DFT(CNBspeech, Frame_length, HiddenVector_length, BinarySignal_length, Strength_DataHiding)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Function outline
% This function extracts a binary signal from CNB speech, which is embedded
% into the amplitude spectrum in the bandwidth of 3.4-4.6 kHz based on DHT domain.
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%   Variable Setting   %
%%%%%%%%%%%%%%%%%%%%%%%%

% The number of frames
Number_of_Frame = length(CNBspeech)/(Frame_length);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          %
%   Extract binary signal from CNB speech  %
%                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PN matrix
PN_Sequence = hadamard(HiddenVector_length);
PN_Sequence = PN_Sequence(:, 1:BinarySignal_length);

% Sample index
Sample_index = 1;

% Initialize output
NBspeech = zeros(1, Frame_length*Number_of_Frame);
BinarySignal = char(zeros(Number_of_Frame, BinarySignal_length));


for Frame_index = 1:Number_of_Frame
    
    % CNB Speech signal to be analyzed
    CNBspeech_analyzed = CNBspeech(Sample_index:Sample_index+Frame_length-1);
    
    % DFT 
    DFT_CNBspeech_analyzed = fft(CNBspeech_analyzed);
    Magnitude_CNBspeech_analyzed = abs(DFT_CNBspeech_analyzed);
    Phase_CNBspeech_analyzed = angle(DFT_CNBspeech_analyzed);
   
    % Extract Spread Spectrum Code from Descrete Hartley Transform (3.4-4.6 kHz)
    Half_Magnitude_NBspeech_analyzed = Magnitude_CNBspeech_analyzed(1:Frame_length/2);
    HiddenVector = Half_Magnitude_NBspeech_analyzed(end+1-HiddenVector_length:end)-BinarySignal_length*Strength_DataHiding;
    Half_Magnitude_NBspeech_analyzed(end+1-HiddenVector_length:end) = zeros(1, HiddenVector_length);
    Magnitude_NBspeech_analyzed = [Half_Magnitude_NBspeech_analyzed, 0, Half_Magnitude_NBspeech_analyzed(end:-1:2)];
         
    % Generate NB speech signal
    NBspeech(Sample_index:Sample_index+Frame_length-1) = real(ifft(Magnitude_NBspeech_analyzed .* exp(1i * Phase_CNBspeech_analyzed)));

    
    % Convert hidden vector into binary signal
    SignSignal = sign(PN_Sequence'*HiddenVector');
    BinarySignal(Frame_index, :) = SignSignal_convert_BinarySignal(SignSignal, BinarySignal_length);
        
    Sample_index = Sample_index+Frame_length;
end

% Avoid discontinuity between frames by over-lapping
NBspeech = OverLapping(NBspeech, Frame_length*4);
