function BinarySignal_LSF_UBspeech = Extract_BinarySignal_LSF(Input_UBspeech_8kHz, Frame_length, BinarySignal_length_LSF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Function outline
% This function extracts LSF from UB speech signal with frame-shift and
% down-sampling at a sampling rate of 8 kHz. LSF is converted into a binary
% signal using a coodbook.
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             %
%   Global variable Setting   %
%                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The number of frames
Number_of_Frame = floor(length(Input_UBspeech_8kHz)/(Frame_length));

% Window function
Window_function = hann(Frame_length+Frame_length/2);

% Margin process for window function
Input_UBspeech_8kHz = [zeros(Frame_length/4, 1); Input_UBspeech_8kHz; zeros(Frame_length/4, 1)];

% Pole for autoregressive model
AR_Pole = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                           %
%   Extract a binary signal for LSF from UB speech signal   %
%                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sample index
Sample_index = 1;

% Initialize output
BinarySignal_LSF_UBspeech = char(zeros(Number_of_Frame, BinarySignal_length_LSF));

for Frame_index = 1:Number_of_Frame
    
    % UB Speech signal to be analyzed
    Input_UBspeech_8kHz_analyzed = Input_UBspeech_8kHz(Sample_index:Sample_index+Frame_length+Frame_length/2-1).*Window_function;
    Input_UBspeech_8kHz_analyzed = Input_UBspeech_8kHz_analyzed(Frame_length/4+1:Frame_length/4+Frame_length);
    
    % Avoid system error for Linear Predictive Coding
    if sum(abs(Input_UBspeech_8kHz_analyzed))==0
        Input_UBspeech_8kHz_analyzed = Input_UBspeech_8kHz_analyzed + 0.0001*randn(Frame_length, 1);
    end
    
    % Autoregressive coefficients for UB speech signal    
    AR_UBspeech = lpc(Input_UBspeech_8kHz_analyzed, AR_Pole);
         
    % Convert AR coefficients into LSF
    LSF_UBspeech = poly2lsf(AR_UBspeech);
    
    % Convert LSF into a binary signal
    BinarySignal_LSF_UBspeech(Frame_index, :) = LSF_to_BinarySignal(LSF_UBspeech, BinarySignal_length_LSF);
    
    Sample_index = Sample_index+Frame_length;
end
