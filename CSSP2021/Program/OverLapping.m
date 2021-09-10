function NBspeech_OL = OverLapping(NBspeech, Frame_length)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Function outline
% This function acts the overlapping process for the NB speech signal to
% avoid the discontinunity between frames.
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%   Variable Setting   %
%%%%%%%%%%%%%%%%%%%%%%%%

% The length of the bandwidth where the hidden vector was embedded
ZeroBandwidth_length = round(600/(8000/Frame_length));

% The number of frames
Number_of_frame = length(NBspeech)/Frame_length;

% Overlapping rate (Defalt:7/8 Overlapping)
OverLapping_length = 8;

% Window function
Window_function = hann(Frame_length);

% Window function weight
Window_function_weight = 0;

for OverLapping_index=1:OverLapping_length
    Window_function_weight = Window_function_weight + Window_function(1+(OverLapping_index-1)*Frame_length/OverLapping_length);
end

% Speech index
Sample_index = 1;

% Initialize output
NBspeech_OL = zeros(size(NBspeech, 1), size(NBspeech, 2));

for Frame_index=1:Number_of_frame*OverLapping_length-OverLapping_length+1
    
    % NB speech signal to be analyzed
    NBspeech_analyzed = NBspeech(Sample_index:Sample_index+Frame_length-1)'.*Window_function;
    
    % DFT
    DFT_NBspeech_analyzed = fft(NBspeech_analyzed);
    Magnitude_NBspeech_analyzed = abs(DFT_NBspeech_analyzed);
    Phase_NBspeech_analyzed = angle(DFT_NBspeech_analyzed);
    
    % Zero padding
    Harf_Magnitude_NBspeech_analyzed = Magnitude_NBspeech_analyzed(1:Frame_length/2);
    Harf_Magnitude_NBspeech_analyzed(end+1-ZeroBandwidth_length:end) = zeros(ZeroBandwidth_length, 1);
    Magnitude_NBspeech_analyzed = [Harf_Magnitude_NBspeech_analyzed; 0; Harf_Magnitude_NBspeech_analyzed(end:-1:2)];
            
    % Overlapping
    NBspeech_OL(Sample_index:Sample_index+Frame_length-1) = NBspeech(Sample_index:Sample_index+Frame_length-1) + real(ifft(Magnitude_NBspeech_analyzed .* exp(1i * Phase_NBspeech_analyzed))')./Window_function_weight;
    
    Sample_index = Sample_index + Frame_length/OverLapping_length;
end