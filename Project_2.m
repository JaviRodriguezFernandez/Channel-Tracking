% Project ATE

clc, clear all, close all;
rng(1);
SNR = [-10,-5,0];
SNR = 30;
snr = 10.^(SNR/10);

NMC = 100;
SIM_MODE = 'OFF-GRID';

SimParams = struct;
SimParams.U = 1;
SimParams.Nms_v = 16*ones(1,SimParams.U);
SimParams.Nbs = 256;
SimParams.Lms_v = 2*ones(1,SimParams.U);
SimParams.Lbs = 8;
SimParams.Gms_v = 8*ones(1,SimParams.U);
SimParams.Gbs = 8;
SimParams.Ktotal = 16;
SimParams.K = SimParams.Ktotal;%min(32,SimParams.Ktotal/4);
SimParams.Kdata = SimParams.Ktotal*(1- 1/8); % 1/8 pilot subcarriers
SimParams.Zp = min(64,ceil(SimParams.K/4));
% SimParams.Zp = SimParams.Ktotal/4;
SimParams.L_v = ones(1,SimParams.U);
SimParams.vel_v = 30*ones(1,SimParams.U);
SimParams.D_v = ones(1,SimParams.U);
SimParams.Ns = 1;
SimParams.A = 'A1A1';
fs = 2560e6;
SimParams.rolloff = 0.25;
SimParams.B = (1+SimParams.rolloff)*fs;
SimParams.Ts = 0.509e-9;
fs = 1/SimParams.Ts;
SimParams.Pwr = 1/SimParams.Ns/SimParams.U*ones(SimParams.K,SimParams.U); %per-user per-subcarrier power
SimParams.Ntr = 4*(SimParams.Ktotal + SimParams.Zp);%1024;
SimParams.Mtck = 2;
SimParams.Bdata = 200;
SimParams.f_c = 28e9;
SimParams.c = 3e8;
SimParams.lambda_c = SimParams.c/SimParams.f_c;
var_n = 1/snr;
SimParams.var_n = var_n;

a_v_Tx = 15;
b_v_Tx = 25;
a_v_Rx = 20;
b_v_Rx = 30;

% Positions and velocities
for u = 1:SimParams.U
    Tx_pos = 70*randn(3,1);
    Tx_pos(3) = abs(Tx_pos(3));
    Rx_pos = 100*randn(3,1);
    Rx_pos(3) = abs(Rx_pos(3));
    SimParams.Tx_pos = Tx_pos;
    SimParams.Rx_pos(:,u) = Rx_pos;
    SimParams.v_Rx(:,u) = sqrt(b_v_Rx-a_v_Rx)*randn(3,1);% + a_v_Rx*sign(randn(3,1));
    SimParams.u = u;
end
SimParams.v_Tx = sqrt(b_v_Tx-a_v_Tx)*randn(3,1);% + a_v_Tx*sign(randn(3,1));

% Tracking time = time_tracking_frames + time_synch_tracking_frames +
% time_synch_data + time_trans_data
SimParams.Ttck = SimParams.Ts*((SimParams.Ktotal + SimParams.Zp + SimParams.Ntr)*SimParams.Mtck + (SimParams.Ktotal + SimParams.Zp + SimParams.Ntr)*SimParams.Bdata);
SimParams.Overhead_Ttck = SimParams.Mtck*(SimParams.Ktotal + SimParams.Zp + SimParams.Ntr);
SimParams.Total_Ttck = SimParams.Overhead_Ttck + SimParams.Bdata*(SimParams.Ktotal + SimParams.Zp + SimParams.Ntr);
SimParams.Efficiency_Ttck = 1 - SimParams.Overhead_Ttck/SimParams.Total_Ttck;
SimParams.Efficiency_Ttck2 = SimParams.Bdata*(SimParams.Ktotal)/(SimParams.Bdata*(SimParams.Ktotal + SimParams.Zp + SimParams.Ntr) + SimParams.Mtck*(SimParams.Ktotal + SimParams.Zp +SimParams.Ntr));

% Space evolution
a_v = -20;
b_v = 20;
for u = 1:SimParams.U
    diff_v(:,u) = (b_v - a_v)*rand(3,1) + a_v;
    diff_v(:,u) = diff_v(:,u).*sign(randn(3,1));
    dx(u) = norm(SimParams.Ttck*diff_v(:,u),2);
end
SimParams.diff_v = diff_v;

% Correlation coefficients
for u = 1:SimParams.U
    SimParams.rho(u) = cos(0.45*dx(u)/SimParams.lambda_c)*exp(-0.1*dx(u)/SimParams.lambda_c);
%     SimParams.rho_NLOS(u) = exp(-0.26*dx(u)/SimParams.lambda_c);
%     fd = SimParams.f_c*norm(diff_v(:,u))/SimParams.c;
%     Tblock = SimParams.Ttck;
%     SimParams.rho(u) = besselj(0,2*pi*fd*Tblock);
end

SimParams.Rho = zeros(SimParams.U*SimParams.K,SimParams.U*SimParams.K);
% Correlation matrix
for u = 1:SimParams.U
    for v = 1:SimParams.U
        for k = 1:SimParams.K
            for l = 1:SimParams.K
                SimParams.Rho((u-1)*SimParams.K + k,(v-1)*SimParams.K + l) = sqrt(SimParams.rho(u)*SimParams.rho(v)*(u == v))*exp(-abs((k-l)*SimParams.B/SimParams.K*SimParams.Ttck));    
            end
        end
    end
end

F = 1/sqrt(SimParams.K)*dftmtx(SimParams.K);
F1 = F(:,1:SimParams.Zp);

NMC = 100;
Nslots = 100;
SimParams.Time = (1:Nslots)*SimParams.Ttck;

P_hat_MC = zeros(SimParams.U*SimParams.K,SimParams.U*SimParams.K,length(SimParams.Time),NMC);
Gains_hat_MC = zeros(SimParams.U*SimParams.K,length(SimParams.Time),NMC);
Gains_MC = zeros(SimParams.U*SimParams.K,length(SimParams.Time),NMC);

for nmc = 1:NMC
    Percentage_MC = nmc/NMC*100
    % Generation of the truth
    SimParams.nmc = nmc;
    SimParams.L = 1;
    SimParams.theta_AoA = zeros(SimParams.U,length(SimParams.Time));
    SimParams.phi_AoA = zeros(size(SimParams.theta_AoA));
    SimParams.theta_AoD = zeros(size(SimParams.theta_AoA));
    SimParams.phi_AoD = zeros(size(SimParams.theta_AoA));
    SimParams.d = zeros(size(SimParams.theta_AoA));
    SimParams.Gains = zeros(SimParams.U,SimParams.K,length(SimParams.Time));
    SimParams.Txp = zeros(3,SimParams.L,length(SimParams.Time));
    SimParams.Rxp = zeros(3,SimParams.U,SimParams.L,length(SimParams.Time));
    SimParams.Txp(:,:,1) = Tx_pos;
    SimParams.Rxp(:,1,:,1) = Rx_pos(:,1);
    Rel_pos = [25,60,15].';
    a_p = -60;
    b_p = 60;
    tau = sort(rand(SimParams.U,1))*(SimParams.Zp - 1)*SimParams.Ts;
    Mfilter = 1;
    q = 50;
    % q = 0.01;
    AS_azi = pi/180;
    AS_ele = pi/180;
    Gt_eff = [];
    
    for u = 1:SimParams.U
        SimParams.u = u;
        Rel_pos(:,u) = (b_p - a_p)*rand(3,1) + a_p;
        Rel_pos(:,u) = Rel_pos(:,u);
        [d_1u,theta_1u,phi_1u] = TransfC2S(Rel_pos(1,u),Rel_pos(2,u),Rel_pos(3,u));
        SimParams.theta_AoA_ini(u) = theta_1u;
        SimParams.theta_AoD_ini(u) = pi - theta_1u;
        SimParams.phi_AoA_ini(u) = phi_1u;
        SimParams.phi_AoD_ini(u) = pi/2 - phi_1u;
        SimParams.d_1u = d_1u;
        
        SimParams.theta_AoA(u,1) = theta_1u + AS_azi/sqrt(2)*randl(1,1);
        SimParams.phi_AoA(u,1) = phi_1u + AS_ele/sqrt(2)*randl(1,1);
        
        std_angle(u) = q*norm(diff_v(:,u),2)^2*SimParams.Ttck/d_1u^2;
        SimParams.std_angle(u) = std_angle(u);
        
        [d_k,theta_k,phi_k,x_k,y_k,z_k] = propagation(SimParams.d_1u,SimParams.theta_AoA(u,1),SimParams.phi_AoA(u,1),SimParams.Ttck,diff_v(1,u),diff_v(2,u),diff_v(3,u));
        SimParams.theta_AoA(u,1) = theta_k + std_angle(u)/sqrt(2)*randl(1,1);
        SimParams.phi_AoA(u,1) = phi_k + std_angle(u)/sqrt(2)*randl(1,1);
        SimParams.phi_AoD(u,1) = pi/2 - SimParams.phi_AoA(u,1);
        SimParams.theta_AoD(u,1) = pi - SimParams.theta_AoA(u,1);
        SimParams.d(u,1) = d_k;
        
        for d = 0:SimParams.Zp-1
            Gt(d+1,u) = sinc((d*SimParams.Ts-tau(u))/Mfilter/SimParams.Ts)*cos(pi*SimParams.rolloff*(d*SimParams.Ts-tau(u))/Mfilter/SimParams.Ts)/(1-(2*SimParams.rolloff*(d*SimParams.Ts-tau(u))/Mfilter/SimParams.Ts)^2);
        end
        F1 = F(:,1:SimParams.Zp)*diag(Gt(:,u));
        SimParams.Pbar_0_Gains(1:SimParams.K,1:SimParams.K,u) = eye(SimParams.K);%F1(:,1:SimParams.Zp)*F1(:,1:SimParams.Zp)'*eye(SimParams.K);%*SimParams.Nbs*SimParams.Nms_v(u);
        SimParams.m_gains_0(u,:,1) = F(:,1:SimParams.Zp)*(Gt(:,u).*(randn(SimParams.Zp,1)+1i*randn(SimParams.Zp,1))*sqrt(1/2));
        Gt_eff = blkdiag(Gt_eff,Gt(:,u));
    end
    
    SimParams.Q_Gains = eye(SimParams.K*SimParams.U) - SimParams.Rho*SimParams.Rho';
    SimParams.Q_Gains = SimParams.Q_Gains;%*kron(eye(SimParams.K),diag(SimParams.Nbs*SimParams.Nms_v));
    
    SimParams.m_gains_0_vect = zeros(SimParams.K*SimParams.U,1);
    SimParams.Pbar_Angles = [];
    SimParams.Q_Angles = [];
    SimParams.Pbar_Gains = [];
    SimParams.Pbar_Gains2 = [];
    for u = 1:SimParams.U
        SimParams.Pbar_Gains = 1e-10*blkdiag(SimParams.Pbar_Gains,SimParams.Pbar_0_Gains(:,:,u));
        SimParams.Pbar_Gains2 = 1e-10*blkdiag(SimParams.Pbar_Gains2,F(:,1:SimParams.Zp)*F(:,1:SimParams.Zp)');%*SimParams.Nbs*SimParams.Nms_v(u));
        covar_u = kron(eye(2),[std_angle(u),0;0,std_angle(u)]);
        SimParams.Pbar_Angles = [AS_azi,0; 0, AS_ele];% blkdiag(SimParams.Pbar_Angles,covar_u);
        covar_u2 = eye(4)*std_angle(u)^2;
        SimParams.Q_Angles = blkdiag(SimParams.Q_Angles,covar_u2);
    end
    SimParams.R_noise = var_n*eye(SimParams.Lbs*SimParams.Mtck*SimParams.K);
    SimParams.Gains_vect(:,1) = SimParams.m_gains_0_vect + sqrtm(SimParams.Pbar_Gains/2)*(randn(SimParams.K*SimParams.U,1) + 1i*randn(SimParams.K*SimParams.U,1));
    SimParams.Gains_vect(:,1) = SimParams.Rho*SimParams.Gains_vect(:,1) + chol(SimParams.Q_Gains/2)*(randn(SimParams.K*SimParams.U,1) + 1i*randn(SimParams.K*SimParams.U,1));
    Channel_gains_energy(1) = norm(SimParams.Gains_vect(:,1),2)^2;
    SimParams = gen_Multiuser_Dict(SimParams,SimParams.theta_AoA(:,1),SimParams.phi_AoA(:,1),SimParams.theta_AoD(:,1),SimParams.phi_AoD(:,1),1);
    SimParams.Channel(:,1) = kron(eye(SimParams.K),SimParams.Psi)*SimParams.Gains_vect(:,1);
    Norm_Energy(1) = norm(SimParams.Channel(:,1),2)^2;
    
    for t = 2:length(SimParams.Time)
        for u = 1:SimParams.U
            [d_k,theta_k,phi_k,x_k,y_k,z_k] = propagation(SimParams.d(u,t-1),SimParams.theta_AoA(u,t-1),SimParams.phi_AoA(u,t-1),SimParams.Ttck,diff_v(1,u),diff_v(2,u),diff_v(3,u));
            SimParams.d(u,t) = d_k;
            SimParams.theta_AoA(u,t) = theta_k + std_angle(u)/sqrt(2)*randl(1,1);%std_angle(u)*randn(1,1);
            SimParams.phi_AoA(u,t) = phi_k + std_angle(u)/sqrt(2)*randl(1,1);%std_angle(u)*randn(1,1);
            SimParams.phi_AoD(u,t) = pi/2 - SimParams.phi_AoA(u,t);
            SimParams.theta_AoD(u,t) = pi - SimParams.theta_AoA(u,t);
        end
        SimParams.Gains_vect(:,t) = SimParams.Rho*SimParams.Gains_vect(:,t-1) + chol(SimParams.Q_Gains/2)*(randn(SimParams.K*SimParams.U,1)+1i*randn(SimParams.K*SimParams.U,1));
        SimParams = gen_Multiuser_Dict(SimParams,SimParams.theta_AoA(:,t),SimParams.phi_AoA(:,t),SimParams.theta_AoD(:,t),SimParams.phi_AoD(:,t),t);
        SimParams.Channel(:,t) = kron(eye(SimParams.K),SimParams.Psi)*SimParams.Gains_vect(:,t);
        Norm_Energy(t) = norm(SimParams.Channel(:,t),2)^2;
    end
    
%     Scaling_factor = SimParams.Nbs*sum(SimParams.Nms_v)*SimParams.K/mean(Norm_Energy);
%     SimParams.Gains_vect = sqrt(Scaling_factor)*SimParams.Gains_vect;
%     SimParams.Channel = sqrt(Scaling_factor)*SimParams.Channel;
    
%     figure, plot(SimParams.Time,squeeze(SimParams.theta_AoA).'*180/pi);
%     figure, plot(SimParams.Time,squeeze(SimParams.phi_AoA).'*180/pi);
    
    
    %%%%%%% DEFINE SENSING MATRIX %%%%%%%%%
    [SimParams]=CSsensingMatrix_MUChannelEstimation(SimParams);
    C_noise = SimParams.W_2'*SimParams.W_2;
    % Combine noise vector
    for t = 1:length(SimParams.Time)
        SimParams.v(:,t) = chol(SimParams.R_noise)*1/sqrt(2)*(randn(SimParams.Mtck*SimParams.Lbs*SimParams.K,1) + 1i*randn(SimParams.Mtck*SimParams.Lbs*SimParams.K,1));
        SimParams.Phi_w(:,:,t) = chol(inv(C_noise))*SimParams.Phi;
        SimParams = gen_Multiuser_Dict(SimParams,SimParams.theta_AoA(:,t),SimParams.phi_AoA(:,t),SimParams.theta_AoD(:,t),SimParams.phi_AoD(:,t),t);
        A = kron(eye(SimParams.K),SimParams.Phi_w(:,:,t)*SimParams.Psi);
%         SimParams.y(:,t) = A*SimParams.Gains_vect(:,t) + SimParams.v(:,t);
        Signal_energy(t) = norm(A*SimParams.Gains_vect(:,t),2)^2;
        SimParams.var_n = Signal_energy(1)/size(A,1)/snr;
        SimParams.R_noise = SimParams.var_n*eye(SimParams.Mtck*SimParams.Lbs*SimParams.K);
        SimParams.v(:,t) = chol(SimParams.R_noise)*1/sqrt(2)*(randn(SimParams.Mtck*SimParams.Lbs*SimParams.K,1) + 1i*randn(SimParams.Mtck*SimParams.Lbs*SimParams.K,1));
        SimParams.y(:,t) = A*SimParams.Gains_vect(:,t) + SimParams.v(:,t);
        Noise_energy(t) = norm(SimParams.v(:,t),2)^2;
        SimParams.SNR(t,nmc) = Signal_energy(t)/Noise_energy(t);
    end
    
%     SNR
    SNR_av = 10*log10(mean(SimParams.SNR))
    %%%%%%% Rao - Blackwellized Particle Filter %%%%%%%%%
    SimParams.N = 200;
    
    SimParams = RBPF_Project(SimParams);
    hat_theta_AoA = SimParams.hat_Theta_AoA;
    hat_theta_AoD = SimParams.hat_Theta_AoD;
    hat_phi_AoA = SimParams.hat_Phi_AoA;
    hat_phi_AoD = SimParams.hat_Phi_AoD;
    hat_Gains = SimParams.hat_Gains;
    P_hat = SimParams.P_hat_Gains;
    Gains_hat_MC(:,:,nmc) = hat_Gains;
    Gains_MC(:,:,nmc) = SimParams.Gains_vect;
    P_hat_MC(:,:,:,nmc) = P_hat;
    
    NMSE(nmc) = norm(SimParams.Channel - SimParams.hat_Channel,'fro')^2/norm(SimParams.Channel,'fro')^2
    
    close all;
    
end

Error = Gains_MC - Gains_hat_MC;

sample_mean = 1/NMC*sum(Error,3);
sample_cov = zeros(SimParams.K*SimParams.U,SimParams.K*SimParams.U,length(SimParams.Time));
for t = 1:length(SimParams.Time)
    for nmc = 1:NMC
        sample_cov(:,:,t) = sample_cov(:,:,t) + 1/NMC*(Error(:,t,nmc) - sample_mean(:,t))*(Error(:,t,nmc) - sample_mean(:,t))';
    end
end

figure(1), hold on,
plot1(1:NMC) = plot(SimParams.Time,squeeze(Error(1,:,:)),'b-','Linewidth',2.0,'Markersize',14);
plot2(1:NMC) = plot(SimParams.Time,3*squeeze(sqrt(P_hat_MC(1,1,:,:))),'r-','Linewidth',2.0,'Markersize',14);
plot3(1:NMC) = plot(SimParams.Time,-3*squeeze(sqrt(P_hat_MC(1,1,:,:))),'r-','Linewidth',2.0,'Markersize',14);
plot4(1) = plot(SimParams.Time,squeeze(sample_mean(1,:)),'k*-','Linewidth',2.0,'Markersize',14);
plot5(1) = plot(SimParams.Time,3*squeeze(sqrt(sample_cov(1,1,:))),'m*-','Linewidth',2.0,'Markersize',14);
plot6(1) = plot(SimParams.Time,-3*squeeze(sqrt(sample_cov(1,1,:))),'m*-','Linewidth',2.0,'Markersize',14);
legend([plot1(1),plot2(1),plot4(1),plot5(1)],{'Estimation error','\pm 3 x Predicted Standard Deviation','Sample mean','\pm 3 x Sample standard deviation'});
xlabel('t')
ylabel('Statistics of estimation error')
title('Component # 1')
axis tight
set(gca,...
  'Box'         , 'on',...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...  
    'LineWidth',2,...
    'FontSize'    , 14);
grid on;

hold off;












