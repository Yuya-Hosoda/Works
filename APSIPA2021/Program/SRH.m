function SRH_Pitch = SRH(Input, Original_VAD, Fs, Frame_length, Frame_shift)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Function outline
% This function estimates a pitch using SRH algorithm.
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Relatied works
%  T. Drugman and A. Alwan,
%  ``Joint robust voicing detection and pitch estimation based on residual harmonics,"
%  in Proc. INTERSPEECH, 2011, pp.1973-1976.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
%   Parameter   %
%%%%%%%%%%%%%%%%%

% Pole pf linear predictive coding (LPC)
LPC_pole = 10;

% Number of the analyzed frames
Number_frame = floor(length(Input)/Frame_shift - (Frame_length/Frame_shift)) + 1;

% Initialize output
SRH_Pitch = zeros(length(Original_VAD), 1);

% Window function
Window_function = hann(Frame_length);

% Number of the harmonics to be calculated
Harmonic_max = 5;


%%%%%%%%%%%%%%%%%%%%%
%   SRH algorithm   %
%%%%%%%%%%%%%%%%%%%%%

for Frame_index = 1:Number_frame
    
    % Pitch rage (50-400 Hz)
    Pitch_min = 50;
    Pitch_max = 400;

    % Input to be analyzed
    Input_analyzed = Input(1+(Frame_index-1)*Frame_shift:Frame_length+(Frame_index-1)*Frame_shift).*Window_function;
   
    % Avoid LPC culculation error
    if sum(abs(Input_analyzed))<eps 
       Input_analyzed = Input_analyzed + eps*randn(Frame_length, 1);
    end
   
    % Culculate LPC
    LPC = lpc(Input_analyzed, LPC_pole);
    
    % Transform LPC to line spectral frequencies (LSF)
    LSF = poly2lsf(LPC);
    
    % Transform LPC to autoregressive (AR) coefficients
    AR_coefficient = lsf2poly(LSF);
    
    % Calculate an excitation signal
    Excitation = filter([0 -AR_coefficient(2:end)], 1, Input_analyzed);

    % Calculate a magnitude of the excitation signal
    Magnitude_Excitation = abs(fft(Excitation));
    
    % Normalization
    Magnitude_Excitation  = Magnitude_Excitation ./sum(Magnitude_Excitation);

    % Calculate sum of the residual harmonics (SRH)
    SRH_Pitch_candidate = zeros(Pitch_max, 1);

    for Pitch_candidate = Pitch_min:Pitch_max
       SRH_Pitch_candidate(Pitch_candidate) = Magnitude_Excitation(round(Pitch_candidate/(Fs/Frame_length))+1);

       % Avoid sub-harmonics
       for Harmonic_now = 2:Harmonic_max
           SRH_Pitch_candidate(Pitch_candidate) = SRH_Pitch_candidate(Pitch_candidate)+ Magnitude_Excitation(round((Pitch_candidate*Harmonic_now)/(Fs/Frame_length))+1) -  Magnitude_Excitation(round((Pitch_candidate*(Harmonic_now-0.5))/(Fs/Frame_length))+1);
       end
    end

    % Select pitch candidate with the highest SRH
    if Original_VAD(Frame_index)>0
       [~, SRH_index] = max(SRH_Pitch_candidate);
       SRH_Pitch(Frame_index) = SRH_index;
    end

    % Averaged pitch from past frames
    Pitch_Mean = mean(SRH_Pitch(SRH_Pitch>0));

    % Update pitch range
    if isnan(Pitch_Mean)
        Pitch_min = 50;
        Pitch_max = 400;
    else
        Pitch_min = round(max(Pitch_Mean/2, 50));
        Pitch_max = round(min(Pitch_Mean*2, 400));
    end

    % Calculate SRH within the updated pitch range
    SRH_Pitch_candidate = zeros(Pitch_max, 1);

    for Pitch_candidate = Pitch_min:Pitch_max
       SRH_Pitch_candidate(Pitch_candidate) = Magnitude_Excitation(round(Pitch_candidate/(Fs/Frame_length))+1);

       for Harmonic_now = 2:Harmonic_max
           SRH_Pitch_candidate(Pitch_candidate) = SRH_Pitch_candidate(Pitch_candidate)+ Magnitude_Excitation(round((Pitch_candidate*Harmonic_now)/(Fs/Frame_length))+1) -  Magnitude_Excitation(round((Pitch_candidate*(Harmonic_now-0.5))/(Fs/Frame_length))+1);
       end
    end

    % Select pitch candidate with the highest SRH
    if Original_VAD(Frame_index)>0
       [~, SRH_index] = max(SRH_Pitch_candidate);
       SRH_Pitch(Frame_index) = SRH_index;
    end      
end

