function BinarySignal_Gain = Extract_BinarySignal_Gain(Input_UBspeech_8kHz, Input_NBspeech_8kHz, BinarySignal_LSF_UBspeech, Frame_length, BinarySignal_length_Gain)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Function outline
% This function calculates a relative gain in each sub-frame and converts it into a binary
% signal using a codebook.
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Global variable Setting   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% import codebook for LSF
global codebook_LSF
codebook = codebook_LSF;

% The number of sub-frames
Number_of_SubFrame = 4;

% Sub-frame
SubFrame_length = floor(Frame_length/Number_of_SubFrame);

% The number of frames
Number_of_Frame = floor(length(Input_UBspeech_8kHz)/(Frame_length));

% Margin process for window function
Input_UBspeech_8kHz = [zeros(Frame_length/4, 1); Input_UBspeech_8kHz; zeros(Frame_length+Frame_length/2, 1)];
Input_NBspeech_8kHz = [zeros(Frame_length/4, 1); Input_NBspeech_8kHz; zeros(Frame_length+Frame_length/2, 1)];

% Window function
Window_function = hann(Frame_length+Frame_length/2);

% Pole for autoregressive model
AR_Pole = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                          %
%   Extract relative gain and convert into binary signal   %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Sample index
Sample_index = 1;

% Initialize output
BinarySignal_Gain = char(zeros(Number_of_Frame, BinarySignal_length_Gain));


for Frame_index = 1:Number_of_Frame
    
    % Initialize matrix for relative gains with sub-frames
    Gain_SubFrames = zeros(Number_of_SubFrame, 1);
    
    for SubFrame_index = 1:Number_of_SubFrame
        
        % UB Speech signal to be analyzed
        Input_UBspeech_8kHz_analyzed = Input_UBspeech_8kHz(Sample_index+(SubFrame_index-1)*SubFrame_length:Sample_index+Frame_length+Frame_length/2-1+(SubFrame_index-1)*SubFrame_length).*Window_function;
        Input_UBspeech_8kHz_analyzed = Input_UBspeech_8kHz_analyzed(Frame_length/4+1:Frame_length/4+Frame_length);

        % Avoid system error for Linear Predictive Coding             
        if sum(abs(Input_UBspeech_8kHz_analyzed))==0
            Input_UBspeech_8kHz_analyzed = Input_UBspeech_8kHz_analyzed + 0.000001*randn(Frame_length, 1);
        end

        % Obtain AR coefficients for UB speech signal from a codebook
        AR_UBspeech = lsf2poly(codebook(bin2dec(BinarySignal_LSF_UBspeech(Frame_index, :))+1, :)');

        % UB excitation signal
        Excitation_UBspeech = Input_UBspeech_8kHz_analyzed-filter([0 -AR_UBspeech(2:end)], 1, Input_UBspeech_8kHz_analyzed);

        % NB Speech signal to be analyzed
        Input_NBspeech_8kHz_analyzed = Input_NBspeech_8kHz(Sample_index+(SubFrame_index-1)*SubFrame_length:Sample_index+Frame_length+Frame_length/2-1+(SubFrame_index-1)*SubFrame_length).*Window_function;
        Input_NBspeech_8kHz_analyzed = Input_NBspeech_8kHz_analyzed(Frame_length/4+1:Frame_length/4+Frame_length);

        % Avoid system error for Linear Predictive Coding             
        if sum(abs(Input_NBspeech_8kHz_analyzed))==0
            Input_NBspeech_8kHz_analyzed = Input_NBspeech_8kHz_analyzed + 0.000001*randn(Frame_length, 1);
        end

        % Autoregressive coefficients for NB speech signal   
        AR_NBspeech = lpc(Input_NBspeech_8kHz_analyzed, AR_Pole);

        % NB excitation signal
        Excitation_NBspeech = Input_NBspeech_8kHz_analyzed-filter([0 -AR_NBspeech(2:end)],1,Input_NBspeech_8kHz_analyzed);

        % Calculate relative gain
        Gain_SubFrames(SubFrame_index) = 20*log10(sum(Excitation_UBspeech.^2))-20*log10(sum(Excitation_NBspeech.^2));

    end
        
    % Convert relative gain into Binary signal 
    BinarySignal_Gain(Frame_index, :) = Gain_to_BinarySignal(Gain_SubFrames, BinarySignal_length_Gain);
    
    Sample_index = Sample_index+Frame_length;
    
end