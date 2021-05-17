function Output_Phase_STFTPI = STFTPI(NB_Phase, F0, VAD, frame, frame_shift, Fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes phase pectrum using STFTPI approach
% Input:  NB_Phase: Existing NB phase spectrum
%         F0: Estimated fundamenta frequencies [Hz]
%         VAD: Voice Active Detection
%         frame: Frame length [sample]
%         frame_shift: Frame shift [sample]
%         Fs: Sampling Frequency [Hz]
% Output: Output_Phase_STFTPI: Rreconstructed WB phase spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Martin Krawczyk, Timo Gerkmann,                     %
%	"STFT Phase Reconstruction in Voiced Speech for     %
%	an Improved Single-Channel Speech Enhancement",     %
%	IEEE/ACM Trans. Audio, Speech, Lang. Process.,	    %
%	vol.22, no.12, pp.1931-1940, 2014.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Y. Wakabayashi, T. Fukumori, M. Nakayama,           %
%   T. Nishiura and Y. Yamashita                        %
%	"Single-channel speech enhancement with phase       %
%   reconstruction based on phase distortion averaging,"%
%   IEEE/ACM Trans. Audio, Speech, Lang. Process.,      %
%   vol.26, no.9, pp.1559--1569, 2018.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


%%%%%%%%%%%%%%%
%   Setting   %
%%%%%%%%%%%%%%%

% FFT half size
FFT_harf_Bin = size(NB_Phase, 1);

% FFT index distance [Hz]
FFT_Freq_distance = Fs/frame;

% Frequency index [Hz]
Freq_Hz = (0:frame/2)*(Fs/frame);
Freq_Hz = Freq_Hz';

% Variable for window function (hann)
a = 0.5;

% Number of analysis frames
Number_frames = size(NB_Phase, 2);

% Initialization
Output_Phase_STFTPI = zeros(FFT_harf_Bin, Number_frames);


%%%%%%%%%%%%%%
%   STFTPI   %
%%%%%%%%%%%%%%

for Frame_now = 1:Number_frames
        
    % VAD switching (V->UV:1 / otherwise:0)
    VAD_Switch =  VAD(Frame_now) > VAD(max(1,Frame_now-1));
    
    % Harmonic half distance [Hz]
    Freq_harf_Hz = F0(Frame_now)/2;
    
    % Harmonic half distance [sample]
    Freq_harf_Bin  = floor(Freq_harf_Hz/FFT_Freq_distance) + 1;

    % Voiced frame
    if VAD(Frame_now) 
        
        % Number of harmonics
        Number_Harmonics  = floor((Fs/2) / F0(Frame_now));  
        
        % Harmonics [Hz]
        Harmonic_Hz = (F0(Frame_now):F0(Frame_now):Number_Harmonics*F0(Frame_now))';
        
        % Harmonics [sample]
        Harmonic_Bin = 1 + round(Harmonic_Hz/FFT_Freq_distance);           
                
        
        if ~VAD_Switch 
            % Harmonics exist in (l-1)-th frame (Voiced frame in (l-1)-th frame)
            % Eq.(7)
            Output_Phase_STFTPI(Harmonic_Bin, Frame_now) = Output_Phase_STFTPI(Harmonic_Bin, Frame_now-1) + 2.*pi.*Harmonic_Hz.*frame_shift/Fs;
        else 
            % Harmonics do not exist in (l-1)-th frame (Unvoiced frame in (l-1)-th frame)
            % Initialization by random phase [-2π, 2π]
            Output_Phase_STFTPI(Harmonic_Bin, Frame_now) = 2*pi*(rand(length(Harmonic_Bin), 1)-0.5);
        end
        
        % DFT for window function
        Phi_W0 = Window_Phase(Freq_Hz(Harmonic_Bin), Harmonic_Hz, frame, Fs, a);
     
        % fundamental frequency > FFT index distance 
        if F0(Frame_now)>FFT_Freq_distance
            % Frequency adjacent to harmonic
            for Bandwidth_Bin = -Freq_harf_Bin:Freq_harf_Bin
                % Maximum frequency adjacent to harmonic
                Bandwidth_Bin_max = max(2, min(Harmonic_Bin+Bandwidth_Bin, FFT_harf_Bin-1));
                % Frequency adjacent to harmonic exists
                if (Bandwidth_Bin ~= 0)&&(~any(Bandwidth_Bin_max == Harmonic_Bin))   
                   % Eq.(8)
                   Phi_W = Window_Phase(Freq_Hz(Bandwidth_Bin_max), Harmonic_Hz, frame, Fs, a);
                   Output_Phase_STFTPI(Bandwidth_Bin_max, Frame_now) = Output_Phase_STFTPI(Harmonic_Bin, Frame_now) - Phi_W0 + Phi_W;
                end
            end
            % Frequency adjacent to harmonic lower than FFT index distance 
            if ((Harmonic_Bin(1)-Freq_harf_Bin) >2) 
                for Bandwidth_Bin = Harmonic_Bin(1)-Freq_harf_Bin-1:-1:2
                    % Eq.(8)
                    Phi_W(1) = Window_Phase(Freq_Hz(Bandwidth_Bin), Harmonic_Hz(1), frame, Fs, a);
                    Output_Phase_STFTPI(Bandwidth_Bin, Frame_now) = Output_Phase_STFTPI(Harmonic_Bin(1), Frame_now) - Phi_W0(1) + Phi_W(1);         
               end
            end
            % Frequency adjacent to harmonic higher than maximun harmonic 
            if ((Harmonic_Bin(end)+Freq_harf_Bin) < FFT_harf_Bin-1)
                for Bandwidth_Bin = Harmonic_Bin(end)+Freq_harf_Bin+1:FFT_harf_Bin-1
                    % Eq.(8)
                    Phi_W(end) = Window_Phase(Freq_Hz(Bandwidth_Bin), Harmonic_Hz(end), frame, Fs, a);
                    Output_Phase_STFTPI(Bandwidth_Bin,Frame_now) = Output_Phase_STFTPI(Harmonic_Bin(end),Frame_now) - Phi_W0(end) + Phi_W(end);
               end
            end
        end   
    end
end

% Phase wrapping
Output_Phase_STFTPI = angle(exp(1i*Output_Phase_STFTPI));

end

% DFT for window function
function Phi_W = Window_Phase(Freq_Bin, Freq_Hz, frame, Fs, a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes DFT for window function
% Input:  Freq_Bin: Frequency [sample]
%         Freq_Hz: Frequency [Hz]
%         frame: length of the window in samples
%         Fs: Sampling frequency [Hz]
%         a: Variable for window function (hann window in a=0.5)
% Output: Phi_W: Phase spetrum for window function using DFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalized angular frequency
Omega = 2*pi*(Freq_Bin-Freq_Hz)/Fs;

% Initialization
Phi_W = zeros(length(Omega), 1);

% Phase spetrum for window function using DFT
for Window_m=1:2
    Phi_W = Phi_W + 0.25*exp(-1i*(frame-1)*Omega/2).*sin(frame*Omega/2).*(exp(-1i*pi*Window_m/frame)./sin(0.5*(Omega-2*pi*Window_m/frame)) + exp(1i*pi*Window_m/frame)./sin(0.5*(Omega+2*pi*Window_m/frame)));
end

% Phase wrapping
Phi_W = angle(Phi_W);

% exclude NaNs such as "0/0"
    if any(isnan(Phi_W))
        Phi_W(sin(0.5*Omega)==0) = 0;
        if a~=1
            Phi_W(sin(0.5*(Omega - 2*pi/frame)) == 0) = -pi;
            Phi_W(sin(0.5*(Omega + 2*pi/frame)) == 0) = +pi;
        end
    end    
end


