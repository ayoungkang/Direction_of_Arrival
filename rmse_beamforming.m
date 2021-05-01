clc
clear 
close all

trial = 2000;
SNRdB = 0:5:25;
SNR=10.^(SNRdB/10);
doas = [0 25]*pi/180; %DOAs in radian
M = 4;     % number of sensors
K = 200;   % number of time snapshots
d = 0.5;   % sensor spacing
noise_var = 1;   % variance of noise
r = length(doas);  

angles= -90:0.05:90;
A=exp(-1i*2*pi*d*(0:M-1)'*sin([doas(:).'])); % steering vector

rmse_music = zeros(1, length(SNR));
rmse_esprit = zeros(1, length(SNR));
rmse_mvdr = zeros(1, length(SNR));
rmse_cbf = zeros(1, length(SNR));

for j = 1 : length(SNR)
    P = ones(1, length(doas))*SNR(j);
    music_est = zeros(trial, length(doas));
    esprit_est = zeros(trial, length(doas));
    mvdr_est = zeros(trial, length(doas));
    cbf_est = zeros(trial, length(doas));
    
    for i = 1 : trial
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
        [pks, loc] = findpeaks(abs(cbf), angles);
        [~, ind] = sort(pks, 'descend');
        cbf_est(i,:) = sort(loc(ind(1:length(doas))));
        
        % MVDR
        IR=inv(R); %Inverse of covariance matrix
        for k=1:length(angles)
            mvdr(k)=1/(a1(:,k)'*IR*a1(:,k));
        end
        
        [pks, loc] = findpeaks(abs(mvdr), angles);
        [~, ind] = sort(pks, 'descend');
        mvdr_est(i,:) = sort(loc(ind(1:length(doas))));
        
        % MUSIC
        [Q, D] = eig(R);   % eigendecomposition of covariance matrix
        [D, I] = sort(diag(D),1,'descend');   
        Q = Q(:,I);   
        Qs = Q(:,1:r);   
        Qn = Q(:,r+1:M);   

        for k=1:length(angles)  
            music_spectrum(k)=(a1(:,k)'*a1(:,k))/(a1(:,k)'*Qn*Qn'*a1(:,k));
        end

        [pks, loc] = findpeaks(abs(music_spectrum), angles);
        [~, ind] = sort(pks, 'descend');
        music_est(i,:) = sort(loc(ind(1:length(doas))));
    

        %ESPRIT Algorithm
        phi= linsolve(Qs(1:M-1,:),Qs(2:M,:));
        ESPRIT_doas=asin(-angle(eig(phi))/(2*pi*d))*180/pi;
        esprit_est(i,:) = sort(ESPRIT_doas);
        E_ind = round(ESPRIT_doas./0.05)+ceil(length(angles)/2);
        ESPRIT = zeros(1, length(angles));
        ESPRIT(E_ind) = 1;
        
    end
    % music
    rmse = zeros(trial, 2);
    for i = 1 : trial
        rmse(i, 1) = (doas(1)*180/pi - music_est(i,1))^2;
        rmse(i, 2) = (doas(2)*180/pi - music_est(i,2))^2;
    end
    rmse_music(j) = 0.5*(sum(rmse(:,1)) + sum(rmse(:,2)))/trial;
    
    % esprit
    rmse = zeros(trial, 2);
    for i = 1 : trial
        rmse(i, 1) = (doas(1)*180/pi - esprit_est(i,1))^2;
        rmse(i, 2) = (doas(2)*180/pi - esprit_est(i,2))^2;
    end
    rmse_esprit(j) = 0.5*(sum(rmse(:,1)) + sum(rmse(:,2)))/trial;
    
    rmse = zeros(trial, 2);
    for i = 1 : trial
        rmse(i, 1) = (doas(1)*180/pi - mvdr_est(i,1))^2;
        rmse(i, 2) = (doas(2)*180/pi - mvdr_est(i,2))^2;
    end
    rmse_mvdr(j) = 0.5*(sum(rmse(:,1)) + sum(rmse(:,2)))/trial;
    
    rmse = zeros(trial, 2);
    for i = 1 : trial
        rmse(i, 1) = (doas(1)*180/pi - cbf_est(i,1))^2;
        rmse(i, 2) = (doas(2)*180/pi - cbf_est(i,2))^2;
    end
    rmse_cbf(j) = 0.5*(sum(rmse(:,1)) + sum(rmse(:,2)))/trial;
end

figure;
semilogy(SNRdB, rmse_music, '-o', 'linewidth',2.0, 'MarkerSize', 10.0)
hold on
semilogy(SNRdB, rmse_esprit, '-s', 'linewidth',2.0, 'MarkerSize', 10.0)
semilogy(SNRdB, rmse_mvdr, '-d', 'linewidth',2.0, 'MarkerSize', 10.0)
semilogy(SNRdB, rmse_cbf, '-*', 'linewidth',2.0, 'MarkerSize', 10.0)
legend('MUSIC', 'ESPRIT', 'MVDR', 'CBF')
grid on
hold off
xlabel('SNR (dB)')
ylabel('RMSE (degree)')
title('Performance of DOA Estimation')
hold off
