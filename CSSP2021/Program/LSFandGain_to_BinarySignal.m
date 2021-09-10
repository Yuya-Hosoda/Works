function BinarySignal_LSFandGain = LSFandGain_to_BinarySignal(LSFandGain, BinarySignal_length)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Function outline
% This function convert LSF and relative gain into a binary signal using a codebook.
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Global variable Setting   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read grobal variable for codebook
global codebook_LSFandGain
codebook = codebook_LSFandGain;

% Select the most approcoate codebook index by Mean Square Error
[~, Codebook_index] = min(sum((codebook - ones(size(codebook, 1), 1)*LSFandGain).^2, 2));

% Convert codebook index(Decimal) into Binary signal
BinarySignal_LSFandGain = dec2bin(Codebook_index-1, BinarySignal_length);
