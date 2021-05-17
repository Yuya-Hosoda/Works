function Output_Phase_Flip = Flip(NB_Phase, frame)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes phase pectrum using duplication of flipped NB phase spectrum
% Input:  NB_Phase: Existing NB phase spectrum
%         frame: frame length [sample]
% Output: Output_Phase_Flip: Rreconstructed WB phase spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Initialization
Output_Phase_Flip = NB_Phase;

% Duplication of flipped NB phase spectrum
Output_Phase_Flip(frame/4+1:frame/2, :) = -1*NB_Phase(frame/4:-1:1, :);
