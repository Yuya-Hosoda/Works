function BinarySignal_LSFGain = Extract_BinarySignal_LSFandGain(Input_UBspeech_8kHz, Input_NBspeech_8kHz, LSF_UBspeech, Frame_length, BinarySignal_length)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Function outline
% This function extracts a relative gain between NB and UB excitation
% signals. LSF and relative gains are grouped and convered into a binary
% signal using a coodbook.
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             %
%   Global variable Setting   %
%                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The number of frames
Number_of_Frame = length(Input_UBspeech_8kHz)/Frame_length;

% Window function
Window_function = hann(Frame_length+Frame_length/2);

% Margin process for window function
Input_UBspeech_8kHz = [zeros(Frame_length/4, 1); Input_UBspeech_8kHz; zeros(Frame_length/4, 1)];
Input_NBspeech_8kHz = [zeros(Frame_length/4, 1); Input_NBspeech_8kHz; zeros(Frame_length/4, 1)];

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
BinarySignal_LSFGain = char(zeros(Number_of_Frame, BinarySignal_length));

for Frame_index = 1:Number_of_Frame
    
    % UB Speech signal to be analyzed
    Input_UBspeech_8kHz_analyzed = Input_UBspeech_8kHz(Sample_index:Sample_index+Frame_length+Frame_length/2-1).*Window_function;
    Input_UBspeech_8kHz_analyzed = Input_UBspeech_8kHz_analyzed(Frame_length/4+1:Frame_length/4+Frame_length);
    
    % Avoid system error for Linear Predictive Coding 
    if sum(abs(Input_UBspeech_8kHz_analyzed))==0
        Input_UBspeech_8kHz_analyzed = Input_UBspeech_8kHz_analyzed + 0.000001*randn(Frame_length, 1);
    end
                    
    % Autoregressive coefficients for UB speech signal   
    AR_UBspeech = lsf2poly(LSF_UBspeech(Frame_index, :)');
        
    % UB excitation signal
    Excitation_UBspeech = Input_UBspeech_8kHz_analyzed-filter([0 -AR_UBspeech(2:end)], 1, Input_UBspeech_8kHz_analyzed);
    
    % UB Speech signal to be analyzed
    Input_NBspeech_8kHz_analyzed = Input_NBspeech_8kHz(Sample_index:Sample_index+Frame_length+Frame_length/2-1).*Window_function;
    Input_NBspeech_8kHz_analyzed = Input_NBspeech_8kHz_analyzed(Frame_length/4+1:Frame_length/4+Frame_length);
    
    % Avoid system error for Linear Predictive Coding     
    if sum(abs(Input_NBspeech_8kHz_analyzed))==0
        Input_NBspeech_8kHz_analyzed = Input_NBspeech_8kHz_analyzed + 0.000001*randn(Frame_length, 1);
    end

    % Autoregressive coefficients for NB speech signal   
    AR_NBspeech = lpc(Input_NBspeech_8kHz_analyzed, AR_Pole);
        
    % NB excitation signal
    Excitation_NBspeech = Input_NBspeech_8kHz_analyzed-filter([0 -AR_NBspeech(2:end)],1,Input_NBspeech_8kHz_analyzed);
    
    % Relative gain
    Gain = 20*log10(sum(Excitation_UBspeech.^2))-20*log10(sum(Excitation_NBspeech.^2));
    
    % Group feature vectors
    LSFandGain = [LSF_UBspeech(Frame_index, :), Gain];
    
    % Convert LSF and Gain into Binary signal 
    BinarySignal_LSFGain(Frame_index, :) = LSFandGain_to_BinarySignal(LSFandGain, BinarySignal_length);
    
    Sample_index = Sample_index+Frame_length;
end