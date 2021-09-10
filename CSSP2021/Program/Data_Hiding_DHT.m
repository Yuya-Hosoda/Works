function CNBspeech = Data_Hiding_DHT(Input_NBspeech_8kHz, BinarySignal_LSFGain, Frame_length, Strength_DataHiding, HiddenVector_length)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Function outline
% This function embeds a binary signal in to an amplitude spectrum of the NB speech signal
% in the bandwidth of 3.4-4.6 kHz using speech steganography based on the discrete Hartley transform domain.
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Global variable Setting   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% The number of frames
Number_of_Frame = length(Input_NBspeech_8kHz)/(Frame_length);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                              %
%   Embed binary signal using data hiding based on DHT domain  %
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize output
CNBspeech = zeros(1, Frame_length*Number_of_Frame);

% Sample index
Sample_index = 1;

% Length of binary signal
BinarySignal_length = size(BinarySignal_LSFGain, 2);

% Generate PN sequense
PN_Sequence =  hadamard(HiddenVector_length);
PN_Sequence = PN_Sequence(:, 1:BinarySignal_length);

for Frame_index = 1:Number_of_Frame
    
    % NB Speech signal to be analyzed
    Input_NBspeech_8kHz_analyzed = Input_NBspeech_8kHz(Sample_index:Sample_index+Frame_length-1);
    
    % Descrete Hartley Transform
    DHT_Input_NBspeech_8kHz_analyzed = real(fft((1+1i).*Input_NBspeech_8kHz_analyzed));
    
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
    HiddenVector = Strength_DataHiding*PN_Sequence*BinarySignal';
    
    % Data hiding based on DHT domain
    DHT_HiddenVector_Input_NBspeech_8kHz_analyzed = DHT_Input_NBspeech_8kHz_analyzed;
    DHT_HiddenVector_Input_NBspeech_8kHz_analyzed(Frame_length/2+1-HiddenVector_length/2:Frame_length/2+HiddenVector_length/2) = HiddenVector;
    
    % Generate CNB speech signal using Inverse Descrete Hartley Transform    
    CNBspeech(Sample_index:Sample_index+Frame_length-1) = real(ifft((1-1i).*DHT_HiddenVector_Input_NBspeech_8kHz_analyzed));
        
    Sample_index = Sample_index+Frame_length;
end