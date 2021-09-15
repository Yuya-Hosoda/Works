function Ans = Evaluation(truePitch, EstPitch)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              
% Function outline
% This function calculates Gross pitch error rate (GPE) and octave erorr rate (OER).
%                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
%   Parameter   %
%%%%%%%%%%%%%%%%%

% Lowest pitch range
Pitch_low = 50;

% Acceptable error[%]
AC_error = 0.20;

% Calcurate GPE
Number_Voice_frame = 0;
Number_Error_frame = 0;

for frame_index = 1:length(truePitch)
    if truePitch(frame_index)>Pitch_low
        Number_Voice_frame = Number_Voice_frame +1;
        if abs(truePitch(frame_index)-EstPitch(frame_index))>truePitch(frame_index)*AC_error
            Number_Error_frame = Number_Error_frame+1;
        end
    end
end

GPE = Number_Error_frame/Number_Voice_frame;


% Calculate OER
Number_Voice_frame = 0;
Number_Error_frame = 0;

for frame_index=1:length(truePitch)
    if truePitch(frame_index)>Pitch_low
        Number_Voice_frame = Number_Voice_frame +1;
        if abs(log2(truePitch(frame_index)/EstPitch(frame_index)))>=0.99
            Number_Error_frame = Number_Error_frame+1;
        end
    end
end

OER = Number_Error_frame/Number_Voice_frame;


Ans = [GPE, OER];