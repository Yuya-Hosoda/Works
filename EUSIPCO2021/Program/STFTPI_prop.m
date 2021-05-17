function [Output_Phase_Prop] = STFTPI_prop(NB_Phase, F0, VAD, frame, frame_shift, Fs, WB_Amp, Output_Phase_Flip)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes phase pectrum using proposed STFTPI approach
% Input:  NB_Phase: Existing NB phase spectrum
%         F0: Estimated fundamenta frequencies [Hz]
%         VAD: Voice Active Detection
%         frame: Frame length [sample]
%         frame_shift: Frame shift [sample]
%         Fs: Sampling Frequency [Hz]
%         WB_Amp: Reconstructed WB amplitude
%         Output_Phase_Flip: Duplication of Flipped NB phase spectral
% Output: Output_Phase_STFTPI: Rreconstructed WB phase spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

% Threshold for fundamental frequency correction
Th_F0_correction = FFT_Freq_distance/2;

% Threshold for harmonic discrimination
Th_Harmonic_Discrimination = 1;

% Initialization
Output_Phase_Prop = Output_Phase_Flip;
Harmonic_Phase_pre = [];
Number_Emphasized_Harmonic = 0;


%%%%%%%%%%%%%%%%%%%%%%%
%   Proposed STFTPI   %
%%%%%%%%%%%%%%%%%%%%%%%

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
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Fundamental freqency estimation error correction   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Phase spectrum for NB harmonics in l-th frame
        Harmonic_Phase_now = NB_Phase(Harmonic_Bin(Harmonic_Hz<4000), Frame_now);
        
        % Number of harmonics in l-th frame coinsides with Number of harmonics in (l-1)-th frame
        if length(Harmonic_Phase_now) == length(Harmonic_Phase_pre)
          
            % Phase difference between harmonics
            Harmonic_Phase_difference_now = Harmonic_Phase_now(2:end)-Harmonic_Phase_now(1:end-1);
            Harmonic_Phase_difference_pre = Harmonic_Phase_pre(2:end)-Harmonic_Phase_pre(1:end-1);
            
            % Fundamental frequency error correction
            Corrected_F0 = (mod(Harmonic_Phase_difference_now-Harmonic_Phase_difference_pre, 2*pi)*Fs/(2*pi*frame_shift));
            
            % Fundamental frequency candidates
            F0_candidate = 0:Fs/frame_shift:4000;
            
            % Corrected fundamental frequency candidates
            Corrected_F0_candidate = F0_candidate+median(Corrected_F0);
            
            % Selection for the corrected fundamental frequency with the smallest distance for estimated fundamental frequency 
            [~, F0_correct_num] = min(abs(Corrected_F0_candidate-F0(Frame_now)));
            F0_correct_prov = Corrected_F0_candidate(F0_correct_num);
            
            if abs(F0_correct_prov-F0(Frame_now))<Th_F0_correction
                % Amount for fundamental frequency correction is within threshold
                % Correct estimated fundamental frequency
                Corrected_F0 = F0_correct_prov;
            else
                % Amount for fundamental frequency correction is outside threshold
                % Keep estimated fundamental frequency
                Corrected_F0 = F0(Frame_now);
            end
        else
            % Number of harmonics in l-th frame does not coinside with Number of harmonics in (l-1)-th frame
            % Keep estimated fundamental frequency
            Corrected_F0 = F0(Frame_now);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Harmonic Discrimination   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Amplitude spectrum for harmonics
        Harmonic_Amp = WB_Amp(Harmonic_Bin, Frame_now);
        
        for Harmonic_num=1:Number_Harmonics
            % Minimum/Maximum adjacent frequency Bin to harmonic
            Delat_min = max(Harmonic_Bin(Harmonic_num)-Freq_harf_Bin, 1);
            Delat_max = min(Harmonic_Bin(Harmonic_num)+Freq_harf_Bin, frame/2);
            
            % Averaged Amplitude
            Averaged_Hamonic_Amp = mean(WB_Amp(Delat_min:Delat_max, Frame_now));
            
            % Amplitude spectrum for harmonic is bigger than Averaged Amplitude
            if Harmonic_Amp(Harmonic_num)/max(eps, Averaged_Hamonic_Amp) > Th_Harmonic_Discrimination
                % Discriminate that the harmonic should be emphasized 
                Number_Emphasized_Harmonic = Number_Emphasized_Harmonic+1;
            end
        end
        
        % Restriction muximum harmonic lower than Fs/2 [Hz]
        if Number_Emphasized_Harmonic*Corrected_F0>Fs/2
            Number_Emphasized_Harmonic = floor(8000/Corrected_F0);
        end
        
        % Corrected harmonics [Hz]
        Corrected_Harmonic_Hz = (Corrected_F0:Corrected_F0:Number_Emphasized_Harmonic*Corrected_F0)';
        
        % Corrected harmonics [sample]
        Corrected_Harmonic_Bin = 1 + round(Corrected_Harmonic_Hz/FFT_Freq_distance); 

        % Memory phase spectrum of harmonics
        Harmonic_Phase_pre = Harmonic_Phase_now;
        
        if ~VAD_Switch             
            % Harmonics exist in (l-1)-th frame (Voiced frame in (l-1)-th frame)
            % Eq.(16)
            Output_Phase_Prop(Corrected_Harmonic_Bin, Frame_now) = Output_Phase_Prop(Corrected_Harmonic_Bin, Frame_now-1) + 2.*pi.*Corrected_Harmonic_Hz.*frame_shift/Fs;
        else 
            % Harmonics do not exist in (l-1)-th frame (Unvoiced frame in (l-1)-th frame)
            % Initialization by random phase [-2π, 2π]
            Output_Phase_Prop(Corrected_Harmonic_Bin, Frame_now) = 2*pi*(rand(length(Corrected_Harmonic_Bin), 1)-0.5);
        end
        
        % DFT for window function
        Phi_W0 = Window_Phase(Freq_Hz(Corrected_Harmonic_Bin), Corrected_Harmonic_Hz, frame, Fs, a);
  
        % fundamental frequency > FFT index distance 
        if F0(Frame_now)>FFT_Freq_distance
            % Frequency adjacent to harmonic
            for Bandwidth_Bin = -Freq_harf_Bin:Freq_harf_Bin
                % Maximum frequency adjacent to harmonic
                Bandwidth_Bin_max = max(2, min(Corrected_Harmonic_Bin+Bandwidth_Bin, FFT_harf_Bin-1));
                % Frequency adjacent to harmonic exists
                if (Bandwidth_Bin~=0)&&(~any(Bandwidth_Bin_max==Corrected_Harmonic_Bin))   
                    % Eq.(17)
                    Phi_W = Window_Phase(Freq_Hz(Bandwidth_Bin_max), Corrected_Harmonic_Hz, frame, Fs, a);
                    Output_Phase_Prop(Bandwidth_Bin_max,Frame_now) = Output_Phase_Prop(Corrected_Harmonic_Bin,Frame_now) - Phi_W0 + Phi_W;
                end
            end
            % Frequency adjacent to harmonic lower than FFT index distance 
            if ((Corrected_Harmonic_Bin(1)-Freq_harf_Bin) >2) 
                for Bandwidth_Bin = Corrected_Harmonic_Bin(1)-Freq_harf_Bin-1:-1:2
                    % Eq.(17)
                    Phi_W(1) = Window_Phase(Freq_Hz(Bandwidth_Bin), Corrected_Harmonic_Hz(1), frame, Fs, a);
                    Output_Phase_Prop(Bandwidth_Bin,Frame_now) = Output_Phase_Prop(Corrected_Harmonic_Bin(1),Frame_now) - Phi_W0(1) + Phi_W(1);                
                end
            end
            % Frequency adjacent to harmonic higher than maximun harmonic 
            if ((Corrected_Harmonic_Bin(end)+Freq_harf_Bin) < FFT_harf_Bin-1)
                for Bandwidth_Bin = Corrected_Harmonic_Bin(end)+Freq_harf_Bin+1:FFT_harf_Bin-1
                    % Eq.(17)
                    Phi_W(end) = Window_Phase(Freq_Hz(Bandwidth_Bin), Corrected_Harmonic_Hz(end), frame, Fs, a);
                    Output_Phase_Prop(Bandwidth_Bin,Frame_now) = Output_Phase_Prop(Corrected_Harmonic_Bin(end),Frame_now) - Phi_W0(end) + Phi_W(end);
                end
            end
        end   
    end
end

% Phase wrapping
Output_Phase_Prop = angle(exp(1i*Output_Phase_Prop));

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


