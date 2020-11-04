% function [f,alpha,CS] = CycSpecFft(x,MapN,fs,M)
% Calculate the Cyclic Spectrum by FFT
% ----------------------------------------------
% Input:
%   x:      signal sequence
%   MapN:   the size of output cyclic spectrum
%   fs:     sample rate
%   M:      smooth length in frequency spectrum, should be even
% Output:
%   f:      frequency meshgrid of cyclic spectrum (size: (MapN+1)*(MapN+1))
%   alpha:  cyclic frequency meshgrid of cyclic spectrum (size: (MapN+1)*(MapN+1))
%   CS:     cyclic spectrum (size: (MapN+1)*(MapN+1))
%  ----------------------------------------------
% Version 2020-11-1
%
function [f,alpha,CS] = CycSpecFft(x,MapN,fs,M)
% X = fft(x,2*MapN).'/sqrt(MapN);
span = lcm(2*MapN,2*length(x));
spanRate = span/(2*MapN);
X0 = fft(x,span)/sqrt(length(x));
X = zeros(1,2*MapN);
for i = 1:2*MapN
    X(i) = X0(1+(i-1)*spanRate);
end
[f,alpha]=meshgrid(-fs/2:fs/MapN:fs/2,-fs/2:fs/MapN:fs/2);
CS = zeros(MapN+1,MapN+1);
for index_a = 0:MapN
    if ((M>1) && (mod(M,2) == 0))
        for m = -M:2:M
            X1 = [X(mod(m+round(MapN/2)+index_a,2*MapN)+1:end),X(1:mod(m+round(MapN/2)+index_a,2*MapN))];
            X2 = [X(mod(m-round(MapN/2)-index_a,2*MapN)+1:end),X(1:mod(m-round(MapN/2)-index_a,2*MapN))];
            X1_down_rate = zeros(1,MapN+1);
            X2_down_rate = zeros(1,MapN+1);
            for i = 1:MapN
                X1_down_rate(i) = X1(2*i-1);
                X2_down_rate(i) = X2(2*i-1);
            end
            CS(index_a+1,:) = CS(index_a+1,:) + 1/(M+1) * X1_down_rate.*conj(X2_down_rate);
        end
    else
        X1 = [X(mod(round(MapN/2)+index_a,2*MapN)+1:end),X(1:mod(round(MapN/2)+index_a,2*MapN))];
        X2 = [X(mod(-round(MapN/2)-index_a,2*MapN)+1:end),X(1:mod(-round(MapN/2)-index_a,2*MapN))];
        X1_down_rate = zeros(1,MapN+1);
        X2_down_rate = zeros(1,MapN+1);
        for i = 1:MapN
            X1_down_rate(i) = X1(2*i-1);
            X2_down_rate(i) = X2(2*i-1);
        end
        CS(index_a+1,:) = CS(index_a+1,:) + X1_down_rate.*conj(X2_down_rate);
    end
end
