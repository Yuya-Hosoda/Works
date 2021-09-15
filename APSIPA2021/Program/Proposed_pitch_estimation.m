function Prop_Pitch = Proposed_pitch_estimation(Input, Original_VAD, Fs, Frame_length, Frame_shift)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Function outline
% This function estimates a pitch using the proposed pitch estimation algorithm.
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Yuya HOSODA, Arata KAWAMURA, Youji IIGUNI,                                                   
% "Pitch Estimation Algorithm for Narrowband Speech Signal using Phase Differences between Harmonics,"  
% in Proc. APSIPA2021,2021 (under review).           
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of analysis frames
Number_frame = floor(length(Input)/Frame_shift - (Frame_length/Frame_shift)) + 1;

% F0 (Output)
Prop_Pitch = zeros(length(Original_VAD), 1);

% Window function (Hanning)
Window_function = hann(Frame_length);

% Threshold of the local minimum possibility for pitch candidate selection
Local_Min_Probability_Th = 0.5;

% Number of frame additive average
Number_frame_additive_average = 3;

% STFT
Input_STFT = stft(Input, Fs, 'Window', Window_function, 'OverlapLength', Frame_length-Frame_shift, 'FFTLength', Frame_length);

% Phase spectrum
Input_Phase = angle(Input_STFT(Frame_length/2+1:end, :));

% Pitch candidates (50-400 Hz)
Harmonic_Freq = 20:160;
Harmonic_Freq = Fs./Harmonic_Freq;

% Lowest harmonics for each pitch candidates on narrow bandwidth (NB) [Hz]
NB_Harmonic_Freq = ceil(300./Harmonic_Freq).*Harmonic_Freq;

% Number of additive average between harmonics
Number_additive_average_Harmonic = 5;

% Frequency index for Harmonics on NB
NB_harmonic = zeros(Number_additive_average_Harmonic, length(Harmonic_Freq));

for Harmonic_num=1:Number_additive_average_Harmonic
    NB_harmonic(Harmonic_num, :) = NB_Harmonic_Freq + (Harmonic_num-1)*Harmonic_Freq;
end

NB_harmonic = round(NB_harmonic./(Fs/Frame_length));

% Similarity
Similarity_additive_average_between_harmonics = zeros(Number_frame-1, length(Harmonic_Freq));

% Viterbi Flag
Vitabi_frag = 0;

% Pitch range
Pitch_max = 400;
Pitch_min = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Proposed algorithm   %
%%%%%%%%%%%%%%%%%%%%%%%%%%

for Frame_index = 1:Number_frame
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Pitch candidate selection using YIN algorithm   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Input
    Input_analyzed = Input((Frame_index-1)*Frame_shift+1:(Frame_index-1)*Frame_shift+Frame_length);
            
    % diference function
    Difference_function = zeros(1, Frame_length/2);
    
    % 差分累積和
    for Lag_index=0:Frame_length/2-1
        for Sample_index=1:Frame_length/2
            Difference_function(Lag_index+1) = Difference_function(Lag_index+1) + (Input_analyzed(Sample_index)-Input_analyzed(Sample_index+Lag_index)).^2;
        end
    end

    % Cumulative mean normalized difference function
    Normalized_Difference_function = ones(1, Frame_length/2);
    
    for Lag_index=1:Frame_length/2-1
        Normalized_Difference_function(Lag_index+1) = Lag_index*Difference_function(Lag_index+1)/(sum(Difference_function(2:Lag_index+1)));
    end
    
    % Local minimum of the cumulative mean normalized difference function
    Lag_index_all = 0:Frame_length/2-1;
    
    % Local minimum probability
    [~, Local_Min_Probability] = islocalmin(Normalized_Difference_function);
    
    % Local minimum possibility > Threshold -> Pitch candidate 
    Pitch_candidate_Lag_index = Lag_index_all(Local_Min_Probability>Local_Min_Probability_Th);
        
    % Pitch candidate between 50-400 Hz
    Pitch_candidate_Lag_index = Pitch_candidate_Lag_index(Pitch_candidate_Lag_index>(Fs/Pitch_max));
    Pitch_candidate_Lag_index = Pitch_candidate_Lag_index(Pitch_candidate_Lag_index<(Fs/Pitch_min));
    
    if isempty(Pitch_candidate_Lag_index)
        [~, Pitch_candidate_Lag_index] = max(Local_Min_Probability(round(Fs/Pitch_max):round(Fs/Pitch_min)));
        Pitch_candidate_Lag_index = Pitch_candidate_Lag_index+round(Fs/Pitch_max)-1;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %   Phase difference   %
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    % Phase spectrum at l frame
    Input_Phase_analyzed = Input_Phase(:, Frame_index);
    
    % Phase spectrum for NB harmonics
    NB_harmonic_Phase = Input_Phase_analyzed(NB_harmonic);

    % Phase difference between harmonic
    NB_Phase_Difference = NB_harmonic_Phase(2:end, :)-NB_harmonic_Phase(1:end-1, :);
    
    if Frame_index>1
        % Phase difference between harmonics on NB
        NB_harmonic_Phase_Difference = NB_Phase_Difference-NB_Phase_Difference_Past;
        
        % Theorotical Phase difference
        Theorotical_Phase_Difference = 2*pi*Harmonic_Freq*Frame_shift/Fs;

        % Similarity between phase difference on NB and thorotical phase difference
        Similarity = zeros(size(NB_harmonic_Phase_Difference, 1), size(NB_harmonic_Phase_Difference, 2));
                
        % Addictive average between harmonics for similarity
        for Harmonic_num = 1:size(NB_harmonic_Phase_Difference, 1)
            Similarity(Harmonic_num, :) = cos(NB_harmonic_Phase_Difference(Harmonic_num, :)-Theorotical_Phase_Difference);
        end
        Similarity_additive_average_between_harmonics(Frame_index-1, :) = mean(Similarity);
    end    
    
    % The picth candidates exsisted in the previous frame
    if Vitabi_frag>0
        if Frame_index>2
            % After third frame
            % Taking frame additive average for similarities
            Similarity_frame_additive_average = mean(Similarity_additive_average_between_harmonics(Frame_index-min(Number_frame_additive_average, Frame_index-1):Frame_index-1, :));
        else
            % At second frame
            Similarity_frame_additive_average = Similarity_additive_average_between_harmonics(Frame_index-1, :);
        end
        
        % Scale convert for similarity : -1-1 -> 0-2
        Similarity_Scale_convert = Similarity_frame_additive_average(Pitch_candidate_Lag_index-19)+1;
        
        % Pitch candidates [Hz]
        Pitch_candidate = Fs./Pitch_candidate_Lag_index;
        
        % Transition Propability
        Transition_probability = zeros(length(Pitch_candidate), length(Pitch_candidate_past));

        for Pitch_candidate_index = 1:length(Pitch_candidate)
            for Pitch_candidate_past_index = 1:length(Pitch_candidate_past)
                Transition_probability(Pitch_candidate_index, Pitch_candidate_past_index) = max(1-(abs((Pitch_candidate(Pitch_candidate_index)-Pitch_candidate_past(Pitch_candidate_past_index))/285))^(0.3679), eps)*Viterbi_score_past(Pitch_candidate_past_index);
            end
        end
        
        % Viterbi algorithm
        Viterbi_score = max(Transition_probability, [], 2).*max(Similarity_Scale_convert', 0);
        
        % Pitch selection
        if Original_VAD(Frame_index)>0
            [~, Pitch_index] =  max(Viterbi_score);
            Prop_Pitch(Frame_index) = Pitch_candidate(Pitch_index);
            
            % Update pitch candidates and Viterbi score
            Pitch_candidate_past = Pitch_candidate;
            Viterbi_score_past = Viterbi_score./sum(Viterbi_score);
        end
        
    else
        % first frame
        % Select lowest pitch candidate using YIN algorithm
        if Original_VAD(Frame_index)>0         
            d_t_dash_harmonic = Normalized_Difference_function(Pitch_candidate_Lag_index+1);
            [~, min_d_num] = min(d_t_dash_harmonic);
            Prop_Pitch(Frame_index) = Fs./Pitch_candidate_Lag_index(min_d_num);
            
            % Update pitch candidates and Viterbi score
            Pitch_candidate_past = Prop_Pitch(Frame_index);
            Viterbi_score_past = 1;
            
            Vitabi_frag = 1;
        end
    end
    
    % Update Phase difference
    NB_Phase_Difference_Past = NB_Phase_Difference;    
end
