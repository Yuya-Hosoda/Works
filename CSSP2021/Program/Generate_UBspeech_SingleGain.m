function UBspeech = Generate_UBspeech_SingleGain(BinarySignal, NBspeech, Frame_length)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Function outline
% This function generates UB speech signal using side information extracted
% from the CNB speech signal. The Ub speech signal is generated using
% a relative gain without overlapping.
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%   Variable Setting   %
%%%%%%%%%%%%%%%%%%%%%%%%

% Import codebook
global codebook_LSFandGain
codebook = codebook_LSFandGain;

% The number of frames
Number_of_frame = length(NBspeech)/Frame_length;

% Window function
Window_function = hann(Frame_length+Frame_length/2);

% Pole of autoregressive model
AR_Pole = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              %
%   Generate UB speech signal  %
%                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize output
UBspeech   = zeros(length(NBspeech), 1);

% Margin process for window function
NBspeech = [zeros(Frame_length/4, 1); NBspeech; zeros(Frame_length/4, 1)];

% Speech sample index
Sample_index = 1;

for Frame_index = 1:Number_of_frame
    
    % NB Speech signal to be analyzed
    NBspeech_analyzed = NBspeech(Sample_index:Sample_index+Frame_length+Frame_length/2-1).*Window_function;
    NBspeech_analyzed = NBspeech_analyzed(Frame_length/4+1:Frame_length/4+Frame_length);

    % Avoid LPC calculation error
    if sum(abs(NBspeech_analyzed))==0
        NBspeech_analyzed = NBspeech_analyzed + 0.000001*randn(Frame_length, 1);
    end    
    
    % Calculate AR coefficients of NB speech signal
    AR_NBspeech = lpc(NBspeech_analyzed, AR_Pole);
    
    % NB excitation signal
    Excitation_NBspeech = NBspeech_analyzed-filter([0 -AR_NBspeech(2:end)],1,NBspeech_analyzed);
    
    % Convert binary signal into side information
    Side_information = codebook(bin2dec(BinarySignal(Frame_index, :))+1, :);
    UBspeech_LSF  = Side_information(1:AR_Pole);
    SingleGain = Side_information(end);
           
    % Generate UB excitation signal
    Excitation_UBspeech = Excitation_NBspeech'.*sqrt(10.^(SingleGain/20));
          
    % Convert LSF into AR coefficients
    UBspeech_AR = lsf2poly(UBspeech_LSF);
        
    % Generate UB speech signal
    UBspeech(Sample_index:Sample_index+Frame_length-1) = filter(1, [1 UBspeech_AR(2:end)], Excitation_UBspeech);
    
    Sample_index = Sample_index+Frame_length;
    
end