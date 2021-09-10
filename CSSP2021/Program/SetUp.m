function SetUp()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Function outline
% This function sets a global variables for codebooks. These codebooks have been generated using
% "Make_Codebook.m".
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Global variable Setting   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global codebook_LSFandGain
global codebook_Gain
global codebook_LSF
 
% Codebook for LSF and a relative gain
codebook_LSFandGain  = importdata('codebook\codebook_LSFandGain.txt');

% Codebook for LSF
codebook_LSF  = importdata('codebook\codebook_LSF.txt');

% Codebook for relative gains
codebook_Gain = importdata('codebook\codebook_Gain.txt');
 