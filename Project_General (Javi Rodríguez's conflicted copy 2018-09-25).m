% Project ATE

clc, clear all, close all;
rng(1);
SNR = [-10,-5,0];
SNR = 0;
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
SimParams.Ktotal = 256;
SimParams.K = min(32,SimParams.Ktotal/4);
SimParams.Kdata = SimParams.Ktotal*(1- 1/8); % 1/8 pilot subcarriers
SimParams.Zp = min(64,ceil(SimParams.K/4));
SimParams.L_v = 4*ones(1,SimParams.U);
SimParams.vel_v = 30*ones(1,SimParams.U);
SimParams.D_v = ones(1,SimParams.U);
SimParams.Ns = 1;
SimParams.A = 'A1A1';
SimParams.Ts = 0.509e-9;
SimParams.rolloff = 0.25;
SimParams.fs = 1/SimParams.Ts;
SimParams.B = (1+SimParams.rolloff)*SimParams.fs/4;
SimParams.Pwr = 1/SimParams.Ns/SimParams.U*ones(SimParams.K,SimParams.U); %per-user per-subcarrier power
SimParams.Ntr = 1024;
SimParams.Mtck = 4;
SimParams.Bdata = 400;
SimParams.f_c = 60e9;
SimParams.c = 3e8;
SimParams.lambda_c = SimParams.c/SimParams.f_c;
var_n = 1/SimParams.U/snr;%/SimParams.Nbs/SimParams.Nms_v(1);
SimParams.var_n = var_n;
SimParams.K_factor = [4,2,1,3];
SimParams.Heights = [60,45,20,15];

a_v_Tx = -25;
b_v_Tx = 25;
a_v_Rx = -30;
b_v_Rx = 30;

SimParams.K_factor_nat = 10.^(SimParams.K_factor/10);

% Positions and velocities
 Tx_pos = 80*randn(3,1);
 Tx_pos(3) = max(abs(Tx_pos(3)),10);
for u = 1:SimParams.U
    Rx_pos = 50*randn(3,1);
    Rx_pos(3) = SimParams.Heights(u);%abs(Rx_pos(3));
    SimParams.Tx_pos = Tx_pos;
    SimParams.Rx_pos(:,u) = Rx_pos;
    SimParams.v_Rx(:,u) = sqrt(b_v_Rx-a_v_Rx)*randn(3,1);% + a_v_Rx*sign(randn(3,1));
    SimParams.u = u;
%     SimParams.K_factor_nat(u) = SimParams.K_factor_nat(u)*Rx_pos(3)/10;
end
SimParams.v_Tx = sqrt(b_v_Tx-a_v_Tx)*randn(3,1);% + a_v_Tx*sign(randn(3,1));

% Tracking time = time_tracking_frames + time_synch_tracking_frames +
% time_synch_data + time_trans_data
SimParams.Ttck = SimParams.Ts*((SimParams.Ktotal + SimParams.Zp + SimParams.Ntr)*SimParams.Mtck + (SimParams.Ktotal + SimParams.Zp + SimParams.Ntr)*SimParams.Bdata);
Training_time = SimParams.Mtck*(SimParams.Zp + SimParams.Ktotal)*SimParams.Ts;
Sync_training_time = SimParams.Mtck*(SimParams.Zp + SimParams.Ntr)*SimParams.Ts;
Sync_data_time = (SimParams.Zp + SimParams.Ntr)*SimParams.Ts;
Data_time = SimParams.Bdata*(SimParams.Zp + SimParams.Ktotal)*SimParams.Ts;
SimParams.Ttck = Sync_training_time + Training_time + Sync_data_time + Data_time;

SimParams.Efficiency_Ttck = Data_time/(Sync_training_time + Training_time + Data_time + Sync_data_time);

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
    for l = 1:SimParams.L_v(1)
%         SimParams.rho(u,l) = cos(0.45*dx(u)/SimParams.lambda_c)*exp(-0.1*dx(u)/SimParams.lambda_c);
        SimParams.rho(u,l) = 0.99*exp(-dx(u)*2.05)*(l == 1) + (0.9*exp(-dx(u)*1.05) - 0.1)*(l > 1);
%         SimParams.rho(u,l) = sqrt(SimParams.Nbs*SimParams.Nms_v(1)/SimParams.L_v(1))*SimParams.rho(u,l);
%     SimParams.rho_NLOS(u) = exp(-0.26*dx(u)/SimParams.lambda_c);
%     fd = SimParams.f_c*norm(diff_v(:,u))/SimParams.c;
%     Tblock = SimParams.Ttck;
%     SimParams.rho(u) = besselj(0,2*pi*fd*Tblock);
    end
end
SimParams.rho_reshaped = diag(reshape(SimParams.rho.',[],1));
SimParams.Rho = zeros(SimParams.U*SimParams.K*SimParams.L_v(1),SimParams.U*SimParams.K*SimParams.L_v(1));
% Correlation matrix
k = 0:(SimParams.K-1);
l = 0:(SimParams.K-1);
% Organization: Users x Subcarriers
SimParams.Rho = kron(exp(-abs(k'-l)*SimParams.B/SimParams.K*SimParams.Ttck),sqrt(SimParams.rho_reshaped));
SimParams.Rho = SimParams.rho_reshaped;
% Organization: Paths x Users x Subcarriers

% for u = 1:SimParams.U
%     for v = 1:SimParams.U
%         for k = 1:SimParams.K
%             for l = 1:SimParams.K
%                 for m = 1:SimParams.L_v(1)
%                     for n = 1:SimParams.L_v(1) 
%                         SimParams.Rho((u-1)*SimParams.K + k,(v-1)*SimParams.K + l) = sqrt(SimParams.rho(u)*SimParams.rho(v)*(u == v))*exp(-abs((k-l)*SimParams.B/SimParams.K*SimParams.Ttck));    
%                     end
%                 end
%             end
%         end
%     end
% end

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
    SimParams.theta_AoA = zeros(SimParams.U,SimParams.L_v(1),length(SimParams.Time));
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
    a_p = -80;
    b_p = 80;
    tau = sort(rand(SimParams.U,1))*(SimParams.Zp - 1)*SimParams.Ts;
    Mfilter = 1;
    q = 500;
%     q = 0.01;
    AS_azi = (pi/180);
    AS_ele = (pi/180);
    SimParams.AS_azi = AS_azi;
    SimParams.AS_ele = AS_ele;
    Gt_eff = [];
    
    for u = 1:SimParams.U
        for l = 1:SimParams.L_v(u)
            SimParams.u = u;
            Rel_pos(:,u,l) = (b_p - a_p)*rand(3,1) + a_p;
            Rel_pos(3,u,l) = (l == 1)*SimParams.Heights(u) + (l > 1)*abs(Rel_pos(3,u,l))*10;
            [d_1u,theta_1u,phi_1u] = TransfC2S(Rel_pos(1,u,l),Rel_pos(2,u,l),Rel_pos(3,u,l));
            SimParams.theta_AoA_ini(u,l) = theta_1u;
            SimParams.theta_AoD_ini(u,l) = pi - theta_1u;
            SimParams.phi_AoA_ini(u,l) = phi_1u;
            SimParams.phi_AoD_ini(u,l) = pi/2 - phi_1u;
            SimParams.d_1u(u,l) = d_1u;
            SimParams.tau_1u(u,l) = d_1u/3e8;
            
            SimParams.theta_AoA(u,l,1) = theta_1u + AS_azi/sqrt(2)*randl(1,1,1);
            SimParams.phi_AoA(u,l,1) = phi_1u + AS_ele/sqrt(2)*randl(1,1,1);
            
            std_angle(u) = q*norm(diff_v(:,u),2)^2*SimParams.Ttck/d_1u^2;
            SimParams.std_angle(u) = std_angle(u);
        end
    end
    
    SimParams.tau_1u = SimParams.tau_1u/max(max(SimParams.tau_1u));
    SimParams.tau_1u = SimParams.tau_1u*(SimParams.Zp - 1)*SimParams.Ts;
    
    for u = 1:SimParams.U
        for l = 1:SimParams.L_v(u)
            [d_k,theta_k,phi_k,x_k,y_k,z_k] = propagation(SimParams.d_1u(u,l),SimParams.theta_AoA(u,l,1),SimParams.phi_AoA(u,l,1),SimParams.Ttck,diff_v(1,u),diff_v(2,u),diff_v(3,u));
            SimParams.theta_AoA(u,l,1) = theta_k + std_angle(u)/sqrt(2)*randl(1,1,1);
            SimParams.phi_AoA(u,l,1) = phi_k + std_angle(u)/sqrt(2)*randl(1,1,1);
            SimParams.phi_AoD(u,l,1) = pi/2 - SimParams.phi_AoA(u,l,1);
            SimParams.theta_AoD(u,l,1) = pi - SimParams.theta_AoA(u,l,1);
            SimParams.d(u,l,1) = d_k;
            SimParams.tau(u,l,1) = d_k/3e8;
            
            for d = 0:SimParams.Zp-1
                Gt(d+1,u) = sinc((d*SimParams.Ts-SimParams.tau_1u(u,l,1))/Mfilter/SimParams.Ts)*cos(pi*SimParams.rolloff*(d*SimParams.Ts-SimParams.tau_1u(u,l,1))/Mfilter/SimParams.Ts)/(1-(2*SimParams.rolloff*(d*SimParams.Ts-SimParams.tau_1u(u,l,1))/Mfilter/SimParams.Ts)^2);
            end
            F1 = F(:,1:SimParams.Zp)*diag(Gt(:,u));
            Scaling_factor = eye(SimParams.Zp);%*diag([1,1/SimParams.K_factor_nat(u)*ones(1,SimParams.Zp-1)]);
            SimParams.Pbar_0_Gains(1:SimParams.K,1:SimParams.K,u,l) = F1*Scaling_factor*F1';%*eye(SimParams.K,SimParams.K)*SimParams.Nbs*SimParams.Nms_v(u);
            SimParams.Pbar_0_Gains2(1:SimParams.K,1:SimParams.K,u,l) = F(:,1:SimParams.Zp)*Scaling_factor*F(:,1:SimParams.Zp)';%*SimParams.Nbs*SimParams.Nms_v(u);
            SimParams.m_gains_0(u,:,1) = F(:,1:SimParams.Zp)*(Gt(:,u).*(randn(SimParams.Zp,1)+1i*randn(SimParams.Zp,1))*sqrt(1/2));
            Gt_eff = blkdiag(Gt_eff,Gt(:,u));
            SimParams.Q_Gains(u,l) = 1 - SimParams.rho(u,l)*SimParams.rho(u,l);%SimParams.Nbs*SimParams.Nms_v(1)/SimParams.L_v(1)*(1 - SimParams.rho(u,l)*SimParams.rho(u,l));
        end
    end
    
%     SimParams.Q_Gains = SimParams.Q_Gains*kron(eye(SimParams.K*SimParams.L_v(1)),diag(SimParams.Nbs*SimParams.Nms_v));%/SimParams.L_v(1));
    
    SimParams.m_gains_0_vect = zeros(SimParams.K,1);%  1/sqrt(2)*(randn(SimParams.K,1) + 1i*randn(SimParams.K,1));%ones(SimParams.K,1);
    SimParams.Pbar_Angles = [];
    SimParams.Q_Angles = [];
    SimParams.Pbar_Gains = [];
    SimParams.Pbar_Gains2 = [];
    
    for u = 1:SimParams.U
        for l = 1:SimParams.L_v(1)
            SimParams.Pbar_Gains = blkdiag(SimParams.Pbar_Gains,SimParams.Pbar_0_Gains(:,:,u,l));
            SimParams.Pbar_Gains2 = blkdiag(SimParams.Pbar_Gains2,SimParams.Pbar_0_Gains2(:,:,u,l));
            covar_u = kron(eye(2),[std_angle(u),0;0,std_angle(u)]);
            SimParams.Pbar_Angles = blkdiag(SimParams.Pbar_Angles,covar_u);
            covar_u2 = eye(4)*std_angle(u)^2;
            SimParams.Q_Angles = blkdiag(SimParams.Q_Angles,covar_u2);
            SimParams.Gains_vect(:,u,l,1) = SimParams.m_gains_0_vect + sqrtm(SimParams.Pbar_0_Gains(:,:,u,l)/2)*(randn(SimParams.K,1) + 1i*randn(SimParams.K,1));
            SimParams.Gains_vect(:,u,l,1) = SimParams.rho(u,l)*SimParams.Gains_vect(:,u,l,1) + sqrtm(SimParams.Q_Gains(u,l)/2)*(randn(SimParams.K,1) + 1i*randn(SimParams.K,1));
        end
    end
    SimParams.R_noise = var_n*eye(SimParams.Lbs*SimParams.Mtck*SimParams.K);

    Channel_gains_energy(1) = sum(sum(sum(abs(SimParams.Gains_vect(:,:,:,1)).^2)));
    SimParams = gen_Multiuser_Dict(SimParams,SimParams.theta_AoA(:,:,1),SimParams.phi_AoA(:,:,1),SimParams.theta_AoD(:,:,1),SimParams.phi_AoD(:,:,1),1);
    for u = 1:SimParams.U
        SimParams.u = u;
        Hk = gen_channel(SimParams,SimParams.theta_AoA(:,:,1),SimParams.theta_AoD(:,:,1),SimParams.phi_AoA(:,:,1),SimParams.phi_AoD(:,:,1),SimParams.Gains_vect(:,:,:,1));
        SimParams.Hk_MU{u,1} = Hk;
    end
    
    SimParams.Q_Gains_w = [];
    SimParams.Q_Gains_w2 = [];
    Energy_gains = zeros(length(SimParams.Time),1);
    Energy_channel = zeros(length(SimParams.Time),1);
    
    for t = 2:length(SimParams.Time)
        Gt = [];
        Q_1 = [];
        Q_2 = [];
        for u = 1:SimParams.U
            for l = 1:SimParams.L_v(1)
                [d_k,theta_k,phi_k,x_k,y_k,z_k] = propagation(SimParams.d(u,l,t-1),SimParams.theta_AoA(u,l,t-1),SimParams.phi_AoA(u,l,t-1),SimParams.Ttck,diff_v(1,u),diff_v(2,u),diff_v(3,u));
                SimParams.d(u,l,t) = d_k;
                SimParams.theta_AoA(u,l,t) = theta_k + std_angle(u)/sqrt(2)*randl(1,1,1);%std_angle(u)*randn(1,1);
                SimParams.phi_AoA(u,l,t) = phi_k + std_angle(u)/sqrt(2)*randl(1,1,1);%std_angle(u)*randn(1,1);
                SimParams.phi_AoD(u,l,t) = pi/2 - SimParams.phi_AoA(u,l,t);
                SimParams.theta_AoD(u,l,t) = pi - SimParams.theta_AoA(u,l,t);
                SimParams.tau(u,l,t) = d_k/3e8;
                for d = 0:SimParams.Zp-1
                    Gt(d+1,u) = sinc((d*SimParams.Ts-SimParams.tau(u,l,t))/Mfilter/SimParams.Ts)*cos(pi*SimParams.rolloff*(d*SimParams.Ts-SimParams.tau(u,l,t))/Mfilter/SimParams.Ts)/(1-(2*SimParams.rolloff*(d*SimParams.Ts-SimParams.tau(u,l,t))/Mfilter/SimParams.Ts)^2);
                end
                Scaling_factor = eye(SimParams.Zp)*diag([1,1/SimParams.K_factor_nat(u)*ones(1,SimParams.Zp-1)]);
                F1 = F(:,1:SimParams.Zp)*diag(Gt(:,u));
                SimParams.Q_Gains_u(1:SimParams.K,1:SimParams.K,u,l) = F1*Scaling_factor*F1';%*sqrt(SimParams.Nbs*SimParams.Nms_v(u));
                SimParams.Q_Gains_u2(1:SimParams.K,1:SimParams.K,u,l) = F(:,1:SimParams.Zp)*Scaling_factor*F(:,1:SimParams.Zp)';%*sqrt(SimParams.Nbs*SimParams.Nms_v(u));
                Q_1 = blkdiag(Q_1,SimParams.Q_Gains_u(:,:,u,l));
                Q_2 = blkdiag(Q_2,SimParams.Q_Gains_u2(:,:,u,l));
                SimParams.Gains_vect(:,u,l,t) = SimParams.rho(u,l)*SimParams.Gains_vect(:,u,l,t-1) + sqrtm(SimParams.Q_Gains(u,l)/2)*(randn(SimParams.K,1) + 1i*randn(SimParams.K,1));
                Energy_gains(t) = sum(sum(sum(abs(SimParams.Gains_vect(:,:,:,t)).^2)));
            end
        end
%         SimParams.Q_Gains_final = SimParams.Q_Gains*Q_1*SimParams.Q_Gains';
%         SimParams.Q_Gains_final2 = SimParams.Q_Gains*Q_2*SimParams.Q_Gains';
%         SimParams.Gains_vect(:,t) = SimParams.Rho*SimParams.Gains_vect(:,t-1) + sqrtm(SimParams.Q_Gains_final/2)*(randn(SimParams.K*SimParams.U*SimParams.L_v(1),1)+1i*randn(SimParams.K*SimParams.U*SimParams.L_v(1),1));
%         SimParams.Channel(:,t) = kron(eye(SimParams.K),SimParams.Psi)*SimParams.Gains_vect(:,t);

    end
    SimParams.Gains_vect = sqrt(SimParams.Nbs*SimParams.Nms_v(1)/SimParams.L_v(1))*SimParams.Gains_vect;
    for t = 1:length(SimParams.Time)
        for u = 1:SimParams.U
                SimParams.u = u;
            Hk = gen_channel(SimParams,SimParams.theta_AoA(:,:,t),SimParams.theta_AoD(:,:,t),SimParams.phi_AoA(:,:,t),SimParams.phi_AoD(:,:,t),SimParams.Gains_vect(:,:,:,t));
            %         SimParams.Hk_MU{u,1} = Hk;
            SimParams.Hk_MU{u,t} = Hk;
            Energy_channel(t) = Energy_channel(t) + sum(sum(sum(abs(Hk).^2)));
        end
    end
    
%     Scaling_factor = SimParams.Nbs*sum(SimParams.Nms_v)*SimParams.K/mean(Norm_Energy);
%     SimParams.Gains_vect = sqrt(Scaling_factor)*SimParams.Gains_vect;
%     SimParams.Channel = sqrt(Scaling_factor)*SimParams.Channel;
    
    figure, plot(SimParams.Time,squeeze(SimParams.theta_AoA(:,1,:)).'*180/pi);
    figure, plot(SimParams.Time,squeeze(SimParams.phi_AoA(:,1,:)).'*180/pi);
    close all
    
    %%%%%%% DEFINE SENSING MATRIX %%%%%%%%%
    [SimParams]=CSsensingMatrix_MUChannelEstimation(SimParams);
    C_noise = SimParams.W_2'*SimParams.W_2;
    SimParams.C_noise = C_noise;
    SimParams.R_noise = SimParams.var_n*eye(SimParams.Mtck*SimParams.Lbs);
    % Combine noise vector
    for t = 1:length(SimParams.Time)
        
        for k = 1:SimParams.K
            H_k = [];
            for u = 1:SimParams.U
                H_dummy = SimParams.Hk_MU{u,t};
                H_k = [H_k H_dummy(:,:,k)];
            end
            
            SimParams.v(:,k,t) = chol(SimParams.R_noise)*1/sqrt(2)*(randn(SimParams.Mtck*SimParams.Lbs,1) + 1i*randn(SimParams.Mtck*SimParams.Lbs,1));
            SimParams.Phi_w(:,:,t) = chol(inv(C_noise))*SimParams.Phi;
            SimParams = gen_Multiuser_Dict(SimParams,SimParams.theta_AoA(:,:,t),SimParams.phi_AoA(:,:,t),SimParams.theta_AoD(:,:,t),SimParams.phi_AoD(:,:,t),t);
            SimParams.y(:,k,t) = SimParams.Phi_w(:,:,t)*reshape(H_k,[],1) + SimParams.v(:,k,t);
            Signal_energy(k,t) = norm(SimParams.y(:,k,t),2)^2;
            Noise_energy(k,t) = norm(SimParams.v(:,k,t),2)^2;
            SNR(k,t,nmc) = Signal_energy(k,t)/Noise_energy(k,t);
        end
    end
    
%     SNR
    SNR_av = 10*log10(mean(mean(SNR,1),2))
    %%%%%%% Rao - Blackwellized Particle Filter %%%%%%%%%
    SimParams.N = 120;
    
    SimParams = RBPF_Project_General(SimParams);
    hat_theta_AoA = SimParams.hat_theta_AoA;
    hat_theta_AoD = SimParams.hat_theta_AoD;
    hat_phi_AoA = SimParams.hat_phi_AoA;
    hat_phi_AoD = SimParams.hat_phi_AoD;
    hat_Gains = SimParams.hat_Gains;
    P_hat = SimParams.P_hat_Gains;
    Gains_hat_MC(:,:,nmc) = hat_Gains;
    Gains_MC(:,:,nmc) = SimParams.Gains_vect;
    P_hat_MC(:,:,:,nmc) = P_hat;
    
    NMSE(nmc) = norm(SimParams.Channel - SimParams.hat_Channel,'fro')^2/norm(SimParams.Channel,'fro')^2
    figure, plot(real(SimParams.theta_AoA.'),'b*-')
    hold on
    plot(real(SimParams.hat_theta_AoA.'),'ks-')
    close all;
    
    figure, plot(real(SimParams.phi_AoA.'),'b*-')
    hold on
    plot(real(SimParams.hat_phi_AoA.'),'ks-')
    
    figure,
    plot(real(SimParams.SE_DL_PCSI),'ms-')
    hold on
    plot(real(SimParams.SE_DL_Est))

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












