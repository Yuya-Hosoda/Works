function CNBspeech = Data_Hiding_DFT(Input_NBspeech_8kHz, BinarySignal_LSFGain, Frame_length, Strength_DataHiding, HiddenVector_length)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Function outline
% This function embeds a binary signal in to a magnitude spectrum of the NB speech signal
% in the bandwidth of 3.4-4 kHz using speech steganography based on the discrete Fourier transform domain.
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Global variable Setting   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The number of frames
Number_of_Frame = length(Input_NBspeech_8kHz)/(Frame_length);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                              %
%   Embed binary signal using data hiding based on DFT domain  %
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize output
CNBspeech = zeros(1, Frame_length*Number_of_Frame);

% Sample index
Sample_index = 1;

% Length of binary signal
BinarySignal_length = size(BinarySignal_LSFGain, 2);

% PN matrix
PN_Sequence =  hadamard(HiddenVector_length);
PN_Sequence = PN_Sequence(:, 1:BinarySignal_length);

for Frame_index = 1:Number_of_Frame
    
    % NB Speech signal to be analyzed
    Input_NBspeech_8kHz_analyzed = Input_NBspeech_8kHz(Sample_index:Sample_index+Frame_length-1);
    
    % Descrete Hartley Transform
    DFT_Input_NBspeech_8kHz_analyzed = fft(Input_NBspeech_8kHz_analyzed);
    Magnitude_Input_NBspeech_8kHz_analyzed = abs(DFT_Input_NBspeech_8kHz_analyzed);
    Phase_Input_NBspeech_8kHz_analyzed = angle(DFT_Input_NBspeech_8kHz_analyzed);
        
    % Convert binary signal into sign signal
    BinarySignal = zeros(1, BinarySignal_length);    
    
    for Bit_index = 1:BinarySignal_length
        if BinarySignal_LSFGain(Frame_index,Bit_index)=='0'
            BinarySignal(Bit_index) = -1;
        else
            BinarySignal(Bit_index) = 1;
        end
    end
    
    % Generate hidden vector
    HiddenVector = PN_Sequence*BinarySignal'+BinarySignal_length;
        
    % Data hiding based on DFT domain
    Half_Magnitude_HiddenVector_Input_NBspeech_8kHz_analyzed = Magnitude_Input_NBspeech_8kHz_analyzed(1:Frame_length/2);
    Half_Magnitude_HiddenVector_Input_NBspeech_8kHz_analyzed(end+1-BinarySignal_length:end) = Strength_DataHiding.*HiddenVector;
    Magnitude_HiddenVector_Input_NBspeech_8kHz_analyzed = [Half_Magnitude_HiddenVector_Input_NBspeech_8kHz_analyzed; 0; Half_Magnitude_HiddenVector_Input_NBspeech_8kHz_analyzed(end:-1:2)];
    
    % Generate CNB speech signal
    CNBspeech(Sample_index:Sample_index+Frame_length-1) = real(ifft(Magnitude_HiddenVector_Input_NBspeech_8kHz_analyzed .* exp(1i * Phase_Input_NBspeech_8kHz_analyzed)));
            
    Sample_index = Sample_index+Frame_length;
end