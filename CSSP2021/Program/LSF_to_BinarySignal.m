function BinarySignal_LSF_UBspeech = LSF_to_BinarySignal(LSF_UBspeech, BinarySignal_length_LSF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Function outline
% This function convert LSF into a binary signal using a codebook.
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Global variable Setting   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read grobal variable for codebook
global codebook_LSF
codebook = codebook_LSF;

% Select the most approcoate codebook index by Mean Square Error
[~, Codebook_index] = min(sum((codebook - ones(size(codebook, 1), 1)*LSF_UBspeech').^2, 2));

% Convert codebook index(Decimal) into Binary signal
BinarySignal_LSF_UBspeech = dec2bin(Codebook_index-1, BinarySignal_length_LSF);


