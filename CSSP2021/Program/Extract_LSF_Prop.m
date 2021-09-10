function [UB_Bit] = Extract_LSF_Prop(UB_Sp, frame, VQ_n, Bit_length)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Global variable Setting   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The Number of an OverLapping
OL_num = 1;

% The number of frames
frame_L = floor(length(UB_Sp)/(frame/OL_num)-OL_num+1);

% hanning window
win = hann(frame+frame/2);

% LPC pole
LPC_Pole = 10;

time = 1;
UB_Bit = char(zeros(frame_L, VQ_n));
UB_Sp = [zeros(frame/4, 1); UB_Sp; zeros(frame/4, 1)];

for frame_n = 1:frame_L
    
    % UB Speech
    now = UB_Sp(time:time+frame+frame/2-1).*win;
    input = now(frame/4+1:frame/4+frame);
    
    if sum(abs(input))==0
        input = input + 0.0001*randn(frame, 1);
    end
    
    % UB LPC (10 pole)    
    AR = lpc(input, LPC_Pole);
         
    % UB spectral envelope (VQ)
    LSF = poly2lsf(AR);
    
    % UB spectral envelope -> VQ
    [UB_Bit(frame_n, :), LSF_VQ] = VQ_converter_LSF(LSF, VQ_n, Bit_length);
    
    time = time+frame/OL_num;
end
