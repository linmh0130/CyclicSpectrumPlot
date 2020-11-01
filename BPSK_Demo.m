clear
% PN1_obj = comm.PNSequence('Polynomial',[3 1 0],'InitialConditions',[0 0 1],'SamplesPerFrame',7);
% PN1_obj = comm.PNSequence('Polynomial',[4 1 0],'InitialConditions',[0 0 1 0],'SamplesPerFrame',15);
% chip1 = PN1*2-1CycSpecFft;
chip1 = ones(8,1);

upRate = 8; % up sampling rate
shape_h =  rcosdesign(0.5,4,upRate,'sqrt');
% shape_h = ones(upRate,1);
chip1_insertZero = zeros(length(chip1)*upRate,1);
for i = 1:length(chip1)
    chip1_insertZero(1+(i-1)*upRate) = chip1(i);
end
chip1_upRate = filter(shape_h,1,chip1_insertZero);
Lchip1 = length(chip1_upRate);

% data = [ones(5,1);-ones(5,1)];
data = [ones(10,1)]; % the same 10 BPSK samples
% data = randi([0,1],16,1)*2-1;
N_data = length(data);
x1 = zeros(N_data*Lchip1,1);
for i = 1:N_data
    x1(1+(i-1)*Lchip1:i*Lchip1) = data(i)*chip1_upRate;
end

t =(0:length(x1)-1)';
x1 = x1.*cos(pi/4*t); % cosine modulation
N1 = length(x1);
fs = 200;
MapN = 200;
[f,alpha,CS1] = CycSpecFft(x1,MapN,fs,8);
toc
figure
mesh(f,alpha,abs(CS1))
xlabel('f/Hz')
ylabel('\alpha/Hz')

figure
f_line=-fs/2:fs/MapN:fs/2;
z = CS1(MapN/2+1,:);
plot(f_line,abs(z))
xlabel('f Hz')
grid on

figure
z = CS1(:,MapN/2+1);
alpha_line=-fs/2:fs/MapN:fs/2;
plot(alpha_line,abs(z)); 
xlabel('\alpha Hz')
grid on
