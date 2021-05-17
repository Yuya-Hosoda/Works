function Output_GL = GL(NB_Phase, Iteration_time, frame, frame_shift, Fs, WB_Amp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes phase pectrum using Griffin-Lim
% Input:  NB_Phase: Existing NB phase spectrum
%         Iteration_time: Iteration time for Griffin-Lim
%         frame: Frame length [sample]
%         frame_shift: Frame shift [sample]
%         Fs: Sampling Frequency [Hz]
%         WB_Amp: Reconstructed WB amplitude
% Output: Output_GL: Rreconstructed WB audio signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   D. Griffin, J. Lim,                                             %
%	"Signal estimation from modified short-time Fourier transform   %
%	IEEE Trans. Acoust., Speech, Signal Process.,                   %
%   vol.ASSP-32, no.2, pp.236-242, 1984.                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialization(Random Phase spectrum)
Output_Phase_half_Random = 2*pi*(rand(size(NB_Phase, 1), size(NB_Phase, 2)));
Output_Phase_half_Random(1:frame/4, :) = NB_Phase(1:frame/4, :);
Output_Phase_Random = [-1*Output_Phase_half_Random(end-1:-1:2, :); Output_Phase_half_Random];
Output_Phase_Random(end, :) = pi*ones(1, size(NB_Phase, 2));
Output_S_Random = WB_Amp.*exp(1i*Output_Phase_Random);
Input_GL = real(istft(Output_S_Random, Fs, 'Window', hann(frame), 'OverlapLength',frame-frame_shift));

% Griffin-Lim
for now=1:Iteration_time
    Input_GL_S = stft(Input_GL, Fs, 'Window', hann(frame), 'OverlapLength', frame-frame_shift, 'FFTLength', frame);

    % Estimate Phase spectrum
    Input_GL_Phase = angle(Input_GL_S);

    % iSTFT
    Output_GL_S = WB_Amp.*exp(1i*Input_GL_Phase);
    Output_GL = real(istft(Output_GL_S, Fs, 'Window', hann(frame), 'OverlapLength',frame-frame_shift));
    
    % Iteration
    Input_GL = Output_GL;
end    
   