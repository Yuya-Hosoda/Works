function UBspeech = Generate_UBspeech_Gains(BinarySignal, NBspeech, Frame_length, BinarySignal_length_LSF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Function outline
% This function generates UB speech signal using side information extracted
% from the CNB speech signal. The Ub speech signal is generated using
% relative gains with overlapping.
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%   Variable Setting   %
%%%%%%%%%%%%%%%%%%%%%%%%

% import codebooks
global codebook_Gain
global codebook_LSF

% Overlapping rate(3/4 Overlapping)
OverLapping_length = 4;

% The number of frame samples
Number_of_frame = length(NBspeech)/Frame_length;

% Frame_shift
FrameShift = floor(Frame_length/OverLapping_length);

% Margin process for window function
NBspeech = [zeros(Frame_length/4, 1); NBspeech; zeros(Frame_length+Frame_length/2, 1)];

% Pole of autoregressive model
AR_Pole = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              %
%   Generate UB speech signal  %
%                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Window function
Window_function = hann(Frame_length+Frame_length/2);

% Window function weight
win_weight = 0;
for i = 1:OverLapping_length 
    win_weight = win_weight + Window_function(Frame_length/4+1+FrameShift*(i-1));
end

% Speech Sample index
Sample_index = 1;

% Initialize output
UBspeech   = zeros(length(NBspeech), 1);
UBspeech = [zeros(Frame_length/4, 1); UBspeech; zeros(Frame_length+Frame_length/2, 1)];

for Frame_index = 1:Number_of_frame

    BinarySignal_LSF  = BinarySignal(Frame_index, 1:BinarySignal_length_LSF);
    BinarySignal_Gain = BinarySignal(Frame_index, BinarySignal_length_LSF+1:end);
    
    UBspeech_LSF  = codebook_LSF(bin2dec(BinarySignal_LSF)+1, :);
    Gains = codebook_Gain(bin2dec(BinarySignal_Gain)+1, :);
    
    for SubFrame_index = 1:OverLapping_length
        
        % NB Speech to be analyzed
        NBspeech_analyzed = NBspeech(Sample_index+(SubFrame_index-1)*FrameShift:Sample_index+Frame_length+Frame_length/2-1+(SubFrame_index-1)*FrameShift).*Window_function;
        NBspeech_analyzed = NBspeech_analyzed(Frame_length/4+1:Frame_length/4+Frame_length);

        % Avoid LPC calculation error
        if sum(abs(NBspeech_analyzed))==0
            NBspeech_analyzed = NBspeech_analyzed + 0.000001*randn(Frame_length, 1);
        end    

        % Calculate AR coefficients for NB speech signal
        AR_NBspeech = lpc(NBspeech_analyzed, AR_Pole);

        % NB excitation signal
        Excitation_NB = NBspeech_analyzed-filter([0 -AR_NBspeech(2:end)],1,NBspeech_analyzed);

        % Generate UB excitation signal
        Excitation_UB = Excitation_NB'.*sqrt(10.^(Gains(SubFrame_index)/20));

        % Convert LSF into AR coefficients
        AR_UBspeech = lsf2poly(UBspeech_LSF);

        % Generate UB speech signal with overlapping
        UBspeech(Sample_index+(SubFrame_index-1)*FrameShift:Sample_index+Frame_length+(SubFrame_index-1)*FrameShift-1) = UBspeech(Sample_index+(SubFrame_index-1)*FrameShift:Sample_index+Frame_length+(SubFrame_index-1)*FrameShift-1) + filter(1, [1 AR_UBspeech(2:end)], Excitation_UB)'./win_weight;
    end
    
    Sample_index = Sample_index+Frame_length;
    
end

% Remove margine for window function
UBspeech = UBspeech(Frame_length/4+1:Frame_length/4+Number_of_frame*Frame_length);