function [d2v2,kern] = getD2v2(sig, order, deltaT)
    [kern, ~] = myYuleWalker2(sig, order);
    [KLc, ~] = getKLc(kern, 1);
    d2v2 = KLc/deltaT*log2(exp(1));
    %output is in units of bits/s
end
