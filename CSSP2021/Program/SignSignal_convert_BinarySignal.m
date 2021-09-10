function BinarySignal = SignSignal_convert_BinarySignal(SSignSignal, BinarySignal_length)

BinarySignal=[];

for Bit_n = 1:BinarySignal_length
    if SSignSignal(Bit_n)>0
        BinarySignal = [BinarySignal, '1'];
    else
        BinarySignal = [BinarySignal, '0'];
    end
end
        