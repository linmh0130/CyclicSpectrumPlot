clear
tic
code0 = ones(1,1);
N_data = 1024;
upRate = 8; % 过采样倍数
% shape_h =  rcosdesign(0.5,6,upRate,'sqrt'); % 根升余弦成形
shape_h = ones(1,upRate); % 方窗成形

testN = 100; % 实验次数，后面做平均
% SNR_dB = -10:1:10;
SNR_dB = 10; 
nB_dB = 0; % 噪声的波动性，可以仿真增加噪声的波动
MapN = 1024; % 画循环谱的大小
cumulative_cyc_spec = zeros(MapN+1,1);
cumulative_CS_H1 = zeros(MapN+1, MapN+1);

SNR = 10^(SNR_dB/10);
amp_noise_dB = (2*rand()-1)*nB_dB;
amp_noise = sqrt(10^(amp_noise_dB/10));

for testTime = 1:testN

data = randi([0,1],1,N_data)*2-1; % 随机数据
data_DSSS = zeros(1,length(code0)*N_data); % 可扩频，但这里其实没有扩频
for i = 1:N_data
    for j = 1:length(code0)
        data_DSSS(j+(i-1)*length(code0)) = code0(j)*data(i);
    end
end
data_insertZero = zeros(1,upRate*length(data_DSSS));

for i = 1:length(data_DSSS) % 内插0
    data_insertZero(upRate*(i-1)+1) = data_DSSS(i);
end
sample_upRate = conv(shape_h,data_insertZero); % 直接过成形滤波器
noise = randn(1,length(sample_upRate))/sqrt(2)/sqrt(SNR);
x = sample_upRate; % 无噪波形
% x = sample_upRate + amp_noise*noise; % 波形加噪

% 计算循环谱
fs = 200;
fixLength = N_data*upRate*length(code0);
[f,alpha,CS_H1] = CycSpecFft(x(1:fixLength),MapN,fs,32);
z = abs(CS_H1(:,MapN/2+1));
cumulative_CS_H1 = cumulative_CS_H1 + CS_H1;
% p = CS_H1(MapN/2+1-MapN/upRate/length(code0),MapN/2+1);
end
cumulative_CS_H1 = cumulative_CS_H1 /testN;
alpha_line=-fs/2:fs/MapN:fs/2;

mesh(f,alpha,abs(cumulative_CS_H1))
xlabel('f / Hz')
ylabel('\alpha / Hz')
toc