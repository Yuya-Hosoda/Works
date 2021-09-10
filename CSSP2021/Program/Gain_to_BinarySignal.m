function BinarySignal_Gain = Gain_to_BinarySignal(Gain_SubFrames, BinarySignal_length_Gain)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Function outline
% This function convert relative gains into a binary signal using a codebook.
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Global variable Setting   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read grobal variable for codebook
global codebook_Gain
codebook = codebook_Gain;

% Select the most approcoate codebook index by Mean Square Error
[~, Codebook_index] = min(sum((codebook - ones(size(codebook, 1), 1)*Gain_SubFrames').^2, 2));

% Convert codebook index(Decimal) into Binary signal
BinarySignal_Gain = dec2bin(Codebook_index-1, BinarySignal_length_Gain);


