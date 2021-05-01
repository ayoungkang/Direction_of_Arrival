clc
clear
close all


doas = [-5 20]*pi/180; %DOAs in radian
M = 4;     % number of sensors
K = 200;   % number of time snapshots
d = 0.5;   % sensor spacing
P = ones(1, length(doas))*10; %Power: 10dB
noise_var = 1;   % variance of noise
r = length(doas);  

angles=(-90:0.05:90);
A=exp(-1i*2*pi*d*(0:M-1)'*sin([doas(:).'])); % steering vector

% signal and noise generation
sig = round(rand(r,K))*2-1; % BPSK
noise = sqrt(noise_var/2)*(randn(M,K)+1i*randn(M,K));  
        
a1=exp(-1i*2*pi*d*(0:M-1)'*sin([angles(:).']*pi/180));
X = A*diag(sqrt(P))*sig+noise;   % data matrix
R = (X*X')/K;  % covariance matrix

% CBF
for k=1:length(angles)
    cbf(k)=(a1(:,k)'*R*a1(:,k));
end
figure;
plot(angles,10*log10(abs(cbf)), 'linewidth',2.0)
xlabel('Angle (deg)')
ylabel ('Power (dB)')
title('Conventional Beamformer')
%title('Comparison of resolution')
grid on

%MVDR
IR=inv(R); %Inverse of covariance matrix
for k=1:length(angles)
    mvdr(k)=1/(a1(:,k)'*IR*a1(:,k));
end
figure;
plot(angles,10*log10(abs(mvdr)), 'linewidth',2.0)
xlabel('Angle (deg)')
ylabel ('Power (dB)')
title('MVDR')
grid on


[Q, D] = eig(R);   % eigendecomposition of covariance matrix
[D, I] = sort(diag(D),1,'descend');   
Q = Q(:,I);   
Qs = Q(:,1:r);   
Qn = Q(:,r+1:M);

% MUSIC 
angles=(-90:0.05:90);


for k=1:length(angles) 
    music_spectrum(k)=1/norm(a1(:,k)'*Qn)^2;
end

figure;
plot(angles,10*log10(abs(music_spectrum)), 'linewidth',2.0)
title('MUSIC Spectrum')
xlabel('Angle (deg)')
ylabel ('Power (dB)')
grid on

%ESPRIT Algorithm
phi= linsolve(Qs(1:M-1,:),Qs(2:M,:));
ESPRIT_doas=asin(-angle(eig(phi))/(2*pi*d))*180/pi;
E_ind = round(ESPRIT_doas./0.05)+ceil(length(angles)/2);
ESPRIT = zeros(1, length(angles));
ESPRIT(E_ind) = 1;
figure;
plot(angles, ESPRIT, '-o')
title('ESPRIT Spectrum')
xlabel('Angle (deg)')
ylabel ('Normalized Power (dB)')
grid on




