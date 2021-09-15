function YIN_Pitch = YIN(Input, Original_VAD, Fs, Frame_length, Frame_shift)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Function outline
% This function estimates a pitch using YIN algorithm.
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Relatied works
%  A. e Cheveign\'{e} and H. Kawahara,
%  ``YIN, a fundamental frequency estimator for speech and music,"
%  J. Acoust. Soc. Am., vol.111, no.4, pp.1917-1930, 2002.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
%   Parameter   %
%%%%%%%%%%%%%%%%%

% Number of the analyzed frames
Number_frame = floor(length(Input)/Frame_shift - (Frame_length/Frame_shift)) + 1;

% Initialize output
YIN_Pitch = zeros(length(Original_VAD), 1);

%%%%%%%%%%%%%%%%%%%%%
%   YIN algorithm   %
%%%%%%%%%%%%%%%%%%%%%

for Frame_index = 1:Number_frame
    
    % Input to be analyzed
    Input_analyzed = Input((Frame_index-1)*Frame_shift+1:(Frame_index-1)*Frame_shift+Frame_length);
            
    % Difference function
    Difference_function = zeros(1, Frame_length/2);
    
    for Lag_index = 0:Frame_length/2-1
        for Sample_index = 1:Frame_length/2
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
    Normalized_Difference_function_Localmin = (Normalized_Difference_function(islocalmin(Normalized_Difference_function)));
    Lag_index_Localmin = (Lag_index_all(islocalmin(Normalized_Difference_function)));
            
    % Ristrict pitch between 50-400 Hz
    Normalized_Difference_function_Localmin = Normalized_Difference_function_Localmin(Lag_index_Localmin>19);
    Lag_index_Localmin = Lag_index_Localmin(Lag_index_Localmin>19);
        
    % Select the smallest lag index of the locla minima
    if isempty(Lag_index_Localmin)
        YIN_Pitch(Frame_index) = 0;
    else
        if Original_VAD(Frame_index)>0
            [~, min_Lag_index_Localmin] = min(Normalized_Difference_function_Localmin);
            YIN_Pitch(Frame_index) = Fs./Lag_index_Localmin(min_Lag_index_Localmin);
        else
            YIN_Pitch(Frame_index) = 0;
        end
    end
       
end

