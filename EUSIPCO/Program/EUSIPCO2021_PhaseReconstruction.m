%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase reconstruction method for artificial bandwidth extension approach
% This program assumes that missing UB amplitude spectrum has been known
% and has been reconstructed perfectly. 
%
% Program branch
%   Flip.m:                    Duplication of Flipped NB phase spectrum
%   GL.m:                      Griffin-Lim
%   STFTPI.m:                  STFTPI approach
%   STFTPI_prop.m(Our method): STFTPI approach with F0 estimation error correction
%                              & harmonic descimination
%
% Input & Output
%   Input: Multiple sinusoidal wave signals(F0:C5-C6)
%   Output:Bandwidth extended WB audio signal
%       Output_Flip:        Duplication of Flipped NB phase spectrum
%       Output_GL:          Griffin-Lim
%       Output_STFTPI:      STFTPI approach
%       Output_STFT_prop:   Proposed method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Yuya HOSODA, Arata KAWAMURA, Youji IIGUNI,                                         %
%	"Phase Reconstruction for Artificial Bandwidth Extension toward Audio Signal,"     %
%	EUSIPCO, pp.-, 2021.(to be accepted)                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%
%   Setting   %
%%%%%%%%%%%%%%%

% Sampling frequency
Fs = 16000;

% Time index
time = 0:1/Fs:1-1/Fs;
time = time';

% Frame length(32ms)
frame = Fs*32/1000;

% Frame shift(4ms)
frame_shift = frame/8;

% Fundamental Frequency
F0_all = [523.251, 587.33, 659.255, 698.456, 783.991, 880, 987.767, 1046.502];

% Input signal
Input = [];

for j=1:8
    
    Input_now = zeros(length(time), 1);
    f0_now = F0_all(j);

    % Number of harmonics
    H_max = floor(8000/f0_now);
    
    % Phase[-2π,2π]
    Phi = (rand(H_max, 1)-0.5)*2*pi;

    for h_num = 1:H_max
        % Attenuation coefficient
        H_bias = 0.9*h_num/(1-H_max)+(0.1-H_max)/(1-H_max);
        Input_now = Input_now + H_bias*sin(2*pi*h_num*f0_now*time+Phi(h_num));
    end

    Input = [Input; Input_now; zeros(frame_shift, 1)];
end

% Normalization
Input = Input./max(abs(Input));
Input = [zeros(frame, 1); Input; zeros(frame, 1)];

% Fumdamental frequency estimation ()
F0_frame_all = pitch(Input, Fs, 'Method', 'LHS', 'Range', [50, 1200], 'WindowLength', round(0.032*Fs),'OverlapLength', round(0.028*Fs), "MedianFilterLength", 10);
F0_frame_all = F0_frame_all';

% Voice active detection(1:Voiced / 0:Unvoiced)
VAD_frame_all = zeros(1, length(F0_frame_all));
VAD_frame_p_all = harmonicRatio(Input,Fs, 'Window',hann(round(Fs.*0.032),'periodic'), 'OverlapLength', round(0.028*Fs));
VAD_frame_p_all = smoothdata(VAD_frame_p_all);

for i=1:length(VAD_frame_p_all)
    if VAD_frame_p_all(i)>0.7
        VAD_frame_all(i)=1;
    end
end


%%%%%%%%%%%
%  STFT   %
%%%%%%%%%%%

Input_S = stft(Input, Fs, 'Window', hann(frame), 'OverlapLength', frame-frame_shift, 'FFTLength', frame);

% Amplitude spectrum
Input_Amp = abs(Input_S);
Input_Amp_half = Input_Amp(frame/2:end, :);

% Phase spectrum
Input_Phase = angle(Input_S);
Input_Phase_half = Input_Phase(frame/2:end, :);

% NB Phase spectrum
Input_Phase_half(frame/4+1:frame/2, :) = zeros(frame/4, size(Input_Phase_half, 2));


%%%%%%%%%%%%%%%%%%
%   Phase Flip   %
%%%%%%%%%%%%%%%%%%

% Duplication of flipped NB phase spectrum
Output_Phase_half_Flip = Flip(Input_Phase_half, frame);

% Estimated WB phase spectrum
Output_Phase_Flip = [-1*Output_Phase_half_Flip(end-1:-1:2, :); Output_Phase_half_Flip];
Output_Phase_Flip(end, :) = pi*ones(1, size(Input_Phase_half, 2));

% iSTFT
Output_S_Flip = Input_Amp.*exp(1i*Output_Phase_Flip);
Output_Flip = real(istft(Output_S_Flip, Fs, 'Window', hann(frame), 'OverlapLength',frame-frame_shift));


%%%%%%%%%%%%%%%%%%%
%   Griffin-Lim   %
%%%%%%%%%%%%%%%%%%%

% Iteration times
Iteration_time = 100;

% Griffin-Lim
Output_GL = GL(Input_Amp, Input_Phase_half, Iteration_time, frame, frame_shift, Fs);


%%%%%%%%%%%%%%%%%%%%%%%
%   Proposed method   %
%%%%%%%%%%%%%%%%%%%%%%%

% STFTPI
Output_Phase_STFTPI_half = STFTPI(Input_Phase_half, F0_frame_all, VAD_frame_all, frame, frame_shift, Fs);

% Estimate Phase spectrum
Output_Phase_STFTPI_half(1:frame/4, :) = Input_Phase_half(1:frame/4, :);
Output_Phase_STFTPI = [Output_Phase_STFTPI_half(end-1:-1:2, :)*-1; Output_Phase_STFTPI_half];
Output_Phase_STFTPI(end, :) = Input_Phase(end, :);

% iSTFT
Output_S_STFTPI = Input_Amp.*exp(1i*Output_Phase_STFTPI);
Output_STFTPI = real(istft(Output_S_STFTPI, Fs, 'Window', hann(frame), 'OverlapLength',frame-frame_shift));


%%%%%%%%%%%%%%%%%%%%%%%
%   Proposed method   %
%%%%%%%%%%%%%%%%%%%%%%%

% STFTPI
Output_Phase_STFTPI_prop_harf = STFTPI_prop(Input_Phase_half, F0_frame_all, VAD_frame_all, frame, frame_shift, Fs, Input_Amp_half, Output_Phase_half_Flip);

% Estimate Phase spectrum
Output_Phase_STFTPI_prop_harf(1:frame/4, :) = Input_Phase_half(1:frame/4, :);
Output_Phase_STFTPI_prop = [Output_Phase_STFTPI_prop_harf(end-1:-1:2, :)*-1; Output_Phase_STFTPI_prop_harf];
Output_Phase_STFTPI_prop(end, :) = Input_Phase(end, :);

% iSTFT
Output_S_STFTPI_prop = Input_Amp.*exp(1i*Output_Phase_STFTPI_prop);
Output_STFTPI_prop = real(istft(Output_S_STFTPI_prop, Fs, 'Window', hann(frame), 'OverlapLength',frame-frame_shift));



