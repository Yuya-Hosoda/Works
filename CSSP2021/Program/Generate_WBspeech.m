function WBspeech = Generate_WBspeech(NBspeech, UBspeech)



UBspeech_US = upsample(UBspeech, 2);
UBspeech_US = lowpass(UBspeech_US, 3400, 16000);

time = 1:length(UBspeech_US);
UBspeech_US = 2*UBspeech_US.*cos(-2*pi*3400*time'/16000);
UBspeech_US = highpass(UBspeech_US, 3400, 16000);

NBspeech_US = upsample(NBspeech, 2);
NBspeech_US = lowpass(NBspeech_US, 3400, 16000);

WBspeech = NBspeech_US + UBspeech_US;