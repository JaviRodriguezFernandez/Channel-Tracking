% Project ATE

clc, clear all, close all;
rng(1);
SNR = [-10,-5,0];
SNR = -10;
snr = 10.^(SNR/10);

NMC = 100;
SIM_MODE = 'OFF-GRID';

SimParams = struct;
SimParams.U = 1;
SimParams.Nms_v = 64*ones(1,SimParams.U);
SimParams.Nbs = 32;
SimParams.Lms_v = 4*ones(1,SimParams.U);
SimParams.Lbs = 4;
SimParams.Gms_v = 8*ones(1,SimParams.U);
SimParams.Gbs = 8;
SimParams.Ktotal = 256;
SimParams.K = min(8,SimParams.Ktotal/4);
SimParams.Kdata = SimParams.Ktotal*(1- 1/8); % 1/8 pilot subcarriers
SimParams.Zp = min(64,ceil(SimParams.K/4));
SimParams.Zp = SimParams.K;
SimParams.L_v = 10*ones(1,SimParams.U);
SimParams.vel_v = 30*ones(1,SimParams.U);
SimParams.D_v = ones(1,SimParams.U);
SimParams.Ns = 2;
SimParams.A = 'A1A1';
SimParams.Ts = 0.509e-9;
SimParams.rolloff = 0.3;
SimParams.fs = 1/SimParams.Ts;
SimParams.B = (1+SimParams.rolloff)*SimParams.fs/4;
SimParams.Pwr = 1/SimParams.Ns/SimParams.U*ones(SimParams.K,SimParams.U); %per-user per-subcarrier power
SimParams.Ntr = 1024;
SimParams.Mtck = 8;
SimParams.Bdata = 200;
SimParams.f_c = 28e9;
SimParams.c = 3e8;
SimParams.lambda_c = SimParams.c/SimParams.f_c;
var_n = 1/snr/SimParams.U;
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
 Rx_pos = 50*randn(3,SimParams.U,SimParams.L_v(1));
for u = 1:SimParams.U
    for l = 1:SimParams.L_v(u)
        Rx_pos(3,u,l) = SimParams.Heights(u)*(l == 1) + abs(Rx_pos(3,u,l))*(l > 1);
        SimParams.Tx_pos = Tx_pos;
        SimParams.Rx_pos(:,u,l) = Rx_pos(:,u,l);
        SimParams.v_Rx(:,u) = sqrt(b_v_Rx-a_v_Rx)*randn(3,1);% + a_v_Rx*sign(randn(3,1));
        SimParams.u = u;
%     SimParams.K_factor_nat(u) = SimParams.K_factor_nat(u)*Rx_pos(3)/10;
    end
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
a_v = -10;
b_v = 10;
for u = 1:SimParams.U
    diff_v(:,u) = a_v*sign(randn(3,1));%(b_v - a_v)*rand(3,1) + a_v;
%     diff_v(:,u) = diff_v(:,u).*sign(randn(3,1));
    dx(u) = norm(SimParams.Ttck*diff_v(:,u),2);
end
SimParams.diff_v = diff_v;

% Correlation coefficients
for u = 1:SimParams.U
    for l = 1:SimParams.L_v(1)
%         SimParams.rho(u,l) = cos(0.45*dx(u)/SimParams.lambda_c)*exp(-0.1*dx(u)/SimParams.lambda_c);
        factor_LOS = SimParams.K_factor_nat(u)/(SimParams.K_factor_nat(u)+1);
        factor_NLOS = 1/(SimParams.K_factor_nat(u)+1);
        SimParams.rho(u,l) = factor_LOS*0.99*exp(-dx(u)*2.05)*(l == 1) + factor_NLOS*(0.9*exp(-dx(u)*1.05) - 0.1)*(l > 1);
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
% SimParams.Rho = kron(exp(-abs(k'-l)*SimParams.B/SimParams.K*SimParams.Ttck),SimParams.rho_reshaped);

% Organization: Paths x Users x Subcarriers

for k = 1:SimParams.K
    for j = 1:SimParams.K
        SimParams.Rho((k-1)*SimParams.L_v(1)*SimParams.U + (1:SimParams.U*SimParams.L_v(1)), (j-1)*SimParams.L_v(1)*SimParams.U + (1:SimParams.U*SimParams.L_v(1))) = exp(-abs((k-j)*SimParams.B/SimParams.K*SimParams.Ttck))*SimParams.rho_reshaped;
    end
end

F = 1/sqrt(SimParams.K)*dftmtx(SimParams.K);
F1 = F(:,1:SimParams.Zp);

NMC = 100;
Nslots = 20;
SimParams.Time = (1:Nslots)*SimParams.Ttck;

P_hat_MC = zeros(SimParams.U*SimParams.K*SimParams.L_v(1),SimParams.U*SimParams.K*SimParams.L_v(1),length(SimParams.Time),NMC);
Gains_hat_MC = zeros(SimParams.U*SimParams.K*SimParams.L_v(1),length(SimParams.Time),NMC);
Gains_MC = zeros(SimParams.U*SimParams.K*SimParams.L_v(1),length(SimParams.Time),NMC);

% Generation of the truth

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
Rel_pos = zeros(3,SimParams.U,SimParams.L_v(1));
Rel_pos2 = zeros(size(Rel_pos));
a_p = -40;
b_p = 40;
tau = sort(rand(SimParams.U,1))*(SimParams.Zp - 1)*SimParams.Ts;
Mfilter = 1;
q = 2;
%     q = 0.01;
AS_azi = 1*(pi/180);
AS_ele = 1*(pi/180);
SimParams.AS_azi = AS_azi;
SimParams.AS_ele = AS_ele;
Gt_eff = [];

for u = 1:SimParams.U
    for l = 1:SimParams.L_v(u)
        SimParams.u = u;
        Rel_pos(:,u,l) = SimParams.Tx_pos - SimParams.Rx_pos(:,u,l);% (b_p - a_p)*rand(3,1) + a_p;
        Rel_pos2(:,u,l) = (SimParams.Rx_pos(:,u,1) - SimParams.Tx_pos)*(l == 1) + (SimParams.Rx_pos(:,u,1) - SimParams.Rx_pos(:,u,l))*(l > 1);
        [d_1u,theta_1u,phi_1u] = TransfC2S(Rel_pos(1,u,l),Rel_pos(2,u,l),Rel_pos(3,u,l));
        [d_1u2,theta_1u2,phi_1u2] = TransfC2S(Rel_pos2(1,u,l),Rel_pos2(2,u,l),Rel_pos2(3,u,l));
        
        SimParams.theta_AoA_ini(u,l) = rand(1,1)*pi;%theta_1u;
        SimParams.theta_AoD_ini(u,l) = rand(1,1)*pi;%theta_1u2;
        SimParams.phi_AoA_ini(u,l) = rand(1,1)*pi/2;%phi_1u;
        SimParams.phi_AoD_ini(u,l) = rand(1,1)*pi/2;%phi_1u2;
        %             SimParams.theta_AoA_ini(u,l) = theta_1u;
        %             SimParams.theta_AoD_ini(u,l) = theta_1u2;
        %             SimParams.phi_AoA_ini(u,l) = phi_1u;
        %             SimParams.phi_AoD_ini(u,l) = phi_1u2;
        SimParams.d_1u(u,l) = d_1u;
        SimParams.d_1u2(u,l) = d_1u2;
        SimParams.tau_1u(u,l) = (d_1u + d_1u2)/3e8;
        
        SimParams.theta_AoA(u,l,1) = SimParams.theta_AoA_ini(u,l) + AS_azi/sqrt(2)*randl(1,1,1);
        SimParams.phi_AoA(u,l,1) = SimParams.phi_AoA_ini(u,l) + AS_ele/sqrt(2)*randl(1,1,1);
        SimParams.theta_AoD(u,l,1) = SimParams.theta_AoD_ini(u,l) + AS_azi/sqrt(2)*randl(1,1,1);
        SimParams.phi_AoD(u,l,1) = SimParams.phi_AoD_ini(u,l) + AS_ele/sqrt(2)*randl(1,1,1);
        
        std_angle(u) = sqrt(q*norm(diff_v(:,u),2)^2*SimParams.Ttck/d_1u^2);
        SimParams.std_angle(u) = std_angle(u);
        
        [d_k,theta_k,phi_k,x_k,y_k,z_k] = propagation(SimParams.d_1u(u,l),SimParams.theta_AoA(u,l,1),SimParams.phi_AoA(u,l,1),SimParams.Ttck,diff_v(1,u),diff_v(2,u),diff_v(3,u));
        SimParams.theta_AoA(u,l,1) = theta_k + std_angle(u)/sqrt(2)*randl(1,1,1);
        SimParams.phi_AoA(u,l,1) = phi_k + std_angle(u)/sqrt(2)*randl(1,1,1);
        SimParams.d1(u,l,1) = d_k;
        [d_k,theta_k,phi_k,x_k,y_k,z_k] = propagation(SimParams.d_1u2(u,l),SimParams.theta_AoD(u,l,1),SimParams.phi_AoD(u,l,1),SimParams.Ttck,diff_v(1,u),diff_v(2,u),diff_v(3,u));
        SimParams.theta_AoD(u,l,1) = theta_k + std_angle(u)/sqrt(2)*randl(1,1,1);
        SimParams.phi_AoD(u,l,1) = phi_k + std_angle(u)/sqrt(2)*randl(1,1,1);
        SimParams.d2(u,l,1) = d_k;
        SimParams.tau(u,l,1) = (SimParams.d1(u,l,1) + SimParams.d2(u,l,1))/3e8;
    end
end

SimParams.tau = SimParams.tau/max(max(SimParams.tau));
SimParams.tau = sort(SimParams.tau,2);
SimParams.tau = SimParams.tau*(SimParams.Zp - 1)*SimParams.Ts;

for u = 1:SimParams.U
    for l = 1:SimParams.L_v(u)
        for d = 0:SimParams.Zp-1
            Gt(d+1,u) = sinc((d*SimParams.Ts-SimParams.tau(u,l,1))/Mfilter/SimParams.Ts)*cos(pi*SimParams.rolloff*(d*SimParams.Ts-SimParams.tau(u,l,1))/Mfilter/SimParams.Ts)/(1-(2*SimParams.rolloff*(d*SimParams.Ts-SimParams.tau(u,l,1))/Mfilter/SimParams.Ts)^2);
        end
        F1 = F(:,1:SimParams.Zp)*diag(Gt(:,u));
        Scaling_factor = eye(SimParams.Zp);%diag([1,1/SimParams.K_factor_nat(u)*ones(1,SimParams.Zp-1)]);
        SimParams.Pbar_0_Gains(1:SimParams.K,1:SimParams.K,u,l) = F1*Scaling_factor*F1';%*eye(SimParams.K,SimParams.K)*SimParams.Nbs*SimParams.Nms_v(u);
        SimParams.Pbar_0_Gains2(1:SimParams.K,1:SimParams.K,u,l) = F(:,1:SimParams.Zp)*Scaling_factor*F(:,1:SimParams.Zp)';%*SimParams.Nbs*SimParams.Nms_v(u);
        SimParams.m_gains_0(u,:,1) = F(:,1:SimParams.Zp)*(Gt(:,u).*(randn(SimParams.Zp,1)+1i*randn(SimParams.Zp,1))*sqrt(1/2));
        Gt_eff = blkdiag(Gt_eff,Gt(:,u));
        factor_LOS = SimParams.K_factor_nat(u)/(SimParams.K_factor_nat(u)+1);
        factor_NLOS = 1/(SimParams.K_factor_nat(u)+1);
        SimParams.Q_Gains(:,:,u,l) = eye(SimParams.K)*(1 - SimParams.rho(u,l)*SimParams.rho(u,l));
        SimParams.Q_Gains(:,:,u,l) = factor_LOS*SimParams.Q_Gains(:,:,u,l)*(l == 1) + factor_NLOS*SimParams.Q_Gains(:,:,u,l)*(l > 1);
    end
end
%     SimParams.Q_Gains = SimParams.Q_Gains*kron(eye(SimParams.K*SimParams.L_v(1)),diag(SimParams.Nbs*SimParams.Nms_v));%/SimParams.L_v(1));

for nmc = 1:NMC
    SimParams.nmc = nmc;
    SimParams.m_gains_0_vect = zeros(SimParams.K*SimParams.U*SimParams.L_v(1),1);
    Gt = [];
    for u = 1:SimParams.U
        gains_u = sqrt(1/2)*(randn(SimParams.L_v(u),1) + 1i*randn(SimParams.L_v(u),1));
        KF_u = SimParams.K_factor_nat(u);
        Pow_LOS = abs(gains_u(1))^2;
        Pow_NLOS = abs(gains_u(2:end)).^2;
        current_KF = 10*log10(Pow_LOS/sum(Pow_NLOS))
        gain_ref = gains_u(1);
        Pow_scale = abs(gain_ref)^2/sum(abs(gains_u(2:end)).^2);
        Scale = KF_u/Pow_scale;
        gains_u(1) = gain_ref*sqrt(Scale);
        Check_KF = 10*log10(abs(gains_u(1))^2/sum(abs(gains_u(2:end)).^2));
        for p = 1:SimParams.L_v(u)
            for d = 0:SimParams.Zp-1
                Gt(d+1,p) = gains_u(p)*sinc((d*SimParams.Ts-SimParams.tau(u,l,1))/Mfilter/SimParams.Ts)*cos(pi*SimParams.rolloff*(d*SimParams.Ts-SimParams.tau(u,l,1))/Mfilter/SimParams.Ts)/(1-(2*SimParams.rolloff*(d*SimParams.Ts-SimParams.tau(u,l,1))/Mfilter/SimParams.Ts)^2);
            end
        end
        gains_u_fd(:,:,u) = F(:,1:SimParams.Zp)*Gt;
    end
%     for k = 1:SimParams.K
%         SimParams.m_gains_0_vect((k-1)*SimParams.L_v(1)*SimParams.U + (1:SimParams.L_v(1)*SimParams.U)) = reshape(gains_u_fd(k,:,:),[],1);
%     end
        
%     SimParams.m_gains_0_vect = 1/sqrt(2)*(randn(SimParams.K*SimParams.U*SimParams.L_v(1),1) + 1i*randn(SimParams.K*SimParams.U*SimParams.L_v(1),1));%ones(SimParams.K,1);
    SimParams.Pbar_Angles = [];
    SimParams.Q_Angles = [];
    SimParams.Pbar_Gains = [];
    SimParams.Pbar_Gains2 = [];
    
    Pbar_Gains = [];
    Pbar_Gains2 = [];
    Q_Gains = [];
    for k = 1:SimParams.K
        for j = 1:SimParams.K
            Pbar_current = squeeze(SimParams.Pbar_0_Gains(k,j,:,:));
            Pbar_Gains((k-1)*SimParams.U*SimParams.L_v(1) + (1:SimParams.U*SimParams.L_v(1)),(j-1)*SimParams.U*SimParams.L_v(1)+(1:SimParams.U*SimParams.L_v(1))) = diag(reshape(Pbar_current.',[],1));
            Pbar_current2 = squeeze(SimParams.Pbar_0_Gains2(k,j,:,:));
            Pbar_Gains2((k-1)*SimParams.U*SimParams.L_v(1) + (1:SimParams.U*SimParams.L_v(1)),(j-1)*SimParams.U*SimParams.L_v(1)+(1:SimParams.U*SimParams.L_v(1))) = diag(reshape(Pbar_current2.',[],1));
            Q_current = squeeze(SimParams.Q_Gains(k,j,:,:));
            Q_Gains((k-1)*SimParams.U*SimParams.L_v(1) + (1:SimParams.U*SimParams.L_v(1)),(j-1)*SimParams.U*SimParams.L_v(1)+(1:SimParams.U*SimParams.L_v(1))) = diag(reshape(Q_current.',[],1));
        end
    end
    
    SimParams.Pbar_0_Gains = Pbar_Gains;
    SimParams.Pbar_0_Gains2 = Pbar_Gains2;
    SimParams.Q_Gains = Q_Gains;
    SimParams.Q_Gains2 = Q_Gains;
    
    for u = 1:SimParams.U
        for l = 1:SimParams.L_v(1)
            covar_u = kron(eye(2),[std_angle(u),0;0,std_angle(u)]);
            SimParams.Pbar_Angles = blkdiag(SimParams.Pbar_Angles,covar_u);
            covar_u2 = eye(4)*std_angle(u)^2;
            SimParams.Q_Angles = blkdiag(SimParams.Q_Angles,covar_u2);
        end
    end
    
    SimParams.Gains_vect(:,1) = SimParams.m_gains_0_vect + sqrtm(SimParams.Pbar_0_Gains/2)*(randn(SimParams.K*SimParams.U*SimParams.L_v(1),1) + 1i*randn(SimParams.K*SimParams.U*SimParams.L_v(1),1));
    SimParams.Gains_vect(:,1) = SimParams.Rho*SimParams.Gains_vect(:,1) + sqrtm(SimParams.Q_Gains/2)*(randn(SimParams.K*SimParams.U*SimParams.L_v(1),1) + 1i*randn(SimParams.K*SimParams.U*SimParams.L_v(1),1));
    
    SimParams.R_noise = var_n*eye(SimParams.Lbs*SimParams.Mtck*SimParams.K);

    Channel_gains_energy(1) = sum(sum(sum(abs(SimParams.Gains_vect(:,:,:,1)).^2)));
    SimParams = gen_Multiuser_Dict(SimParams,SimParams.theta_AoA(:,:,1),SimParams.phi_AoA(:,:,1),SimParams.theta_AoD(:,:,1),SimParams.phi_AoD(:,:,1),1);
%     Gains_vect_1 = reshape(SimParams.Gains_vect(:,1),SimParams.L_v(1),SimParams.U,SimParams.K);
%     for u = 1:SimParams.U
%         SimParams.u = u;
%         Hk = gen_channel(SimParams,SimParams.theta_AoA(:,:,1),SimParams.theta_AoD(:,:,1),SimParams.phi_AoA(:,:,1),SimParams.phi_AoD(:,:,1),Gains_vect_1);
%         SimParams.Hk_MU{u,1} = Hk;
%     end
    
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
                [d_k,theta_k,phi_k,x_k,y_k,z_k] = propagation(SimParams.d1(u,l,t-1),SimParams.theta_AoA(u,l,t-1),SimParams.phi_AoA(u,l,t-1),SimParams.Ttck,diff_v(1,u),diff_v(2,u),diff_v(3,u));
                SimParams.d1(u,l,t) = d_k;
                SimParams.theta_AoA(u,l,t) = theta_k + std_angle(u)/sqrt(2)*randl(1,1,1);%std_angle(u)*randn(1,1);
                SimParams.phi_AoA(u,l,t) = phi_k + std_angle(u)/sqrt(2)*randl(1,1,1);%std_angle(u)*randn(1,1);
                [d_k,theta_k,phi_k,x_k,y_k,z_k] = propagation(SimParams.d2(u,l,t-1),SimParams.theta_AoD(u,l,t-1),SimParams.phi_AoD(u,l,t-1),SimParams.Ttck,diff_v(1,u),diff_v(2,u),diff_v(3,u));
                SimParams.d2(u,l,t) = d_k;
                SimParams.theta_AoD(u,l,t) = theta_k + std_angle(u)/sqrt(2)*randl(1,1,1);%std_angle(u)*randn(1,1);
                SimParams.phi_AoD(u,l,t) = phi_k + std_angle(u)/sqrt(2)*randl(1,1,1);%std_angle(u)*randn(1,1);
                SimParams.tau(u,l,t) = (SimParams.d1(u,l,t) + SimParams.d2(u,l,t))/3e8;
            end
        end
        SimParams.tau(:,:,t) = SimParams.tau(:,:,t)/max(max(SimParams.tau(:,:,t)));
        SimParams.tau(:,:,t) = SimParams.tau(:,:,t)*(SimParams.Zp - 1)*SimParams.Ts;
        for u = 1:SimParams.U
            for l = 1:SimParams.L_v(u)
                for d = 0:SimParams.Zp-1
                    Gt(d+1,u) = sinc((d*SimParams.Ts-SimParams.tau(u,l,t))/Mfilter/SimParams.Ts)*cos(pi*SimParams.rolloff*(d*SimParams.Ts-SimParams.tau(u,l,t))/Mfilter/SimParams.Ts)/(1-(2*SimParams.rolloff*(d*SimParams.Ts-SimParams.tau(u,l,t))/Mfilter/SimParams.Ts)^2);
                end
                Scaling_factor = eye(SimParams.Zp);%*diag([1,1/SimParams.K_factor_nat(u)*ones(1,SimParams.Zp-1)]);
                F1 = F(:,1:SimParams.Zp)*diag(Gt(:,u)) + 0.1*eye(SimParams.K);
                SimParams.Q_Gains_u(1:SimParams.K,1:SimParams.K,u,l) = F1*Scaling_factor*F1';%*sqrt(SimParams.Nbs*SimParams.Nms_v(u));
                SimParams.Q_Gains_u2(1:SimParams.K,1:SimParams.K,u,l) = F(:,1:SimParams.Zp)*Scaling_factor*F(:,1:SimParams.Zp)';%*sqrt(SimParams.Nbs*SimParams.Nms_v(u));
                Q_1 = blkdiag(Q_1,SimParams.Q_Gains_u(:,:,u,l));
                Q_2 = blkdiag(Q_2,SimParams.Q_Gains_u2(:,:,u,l));
            end
        end
        
        Q_Gains_ext = [];
        for k = 1:SimParams.K
            for j = 1:SimParams.K
                Q_current = squeeze(SimParams.Q_Gains_u(k,j,:,:));
                Q_Gains_ext((k-1)*SimParams.U*SimParams.L_v(1) + (1:SimParams.U*SimParams.L_v(1)),(j-1)*SimParams.U*SimParams.L_v(1)+(1:SimParams.U*SimParams.L_v(1))) = diag(reshape(Q_current.',[],1));
            end
        end
        
        SimParams.Q_Gains = Q_Gains;% 10*sqrtm(Q_Gains)*Q_Gains_ext*sqrtm(Q_Gains)';
        
        SimParams.Gains_vect(:,t) = SimParams.Rho*SimParams.Gains_vect(:,t-1) + sqrtm(SimParams.Q_Gains/2)*(randn(SimParams.K*SimParams.U*SimParams.L_v(1),1) + 1i*randn(SimParams.K*SimParams.U*SimParams.L_v(1),1));
        Energy_gains(t) = sum(abs(SimParams.Gains_vect(:,t)).^2);
        
%         SimParams.Q_Gains_final = SimParams.Q_Gains*Q_1*SimParams.Q_Gains';
%         SimParams.Q_Gains_final2 = SimParams.Q_Gains*Q_2*SimParams.Q_Gains';
%         SimParams.Gains_vect(:,t) = SimParams.Rho*SimParams.Gains_vect(:,t-1) + sqrtm(SimParams.Q_Gains_final/2)*(randn(SimParams.K*SimParams.U*SimParams.L_v(1),1)+1i*randn(SimParams.K*SimParams.U*SimParams.L_v(1),1));
%         SimParams.Channel(:,t) = kron(eye(SimParams.K),SimParams.Psi)*SimParams.Gains_vect(:,t);

    end
    SimParams.Gains_vect = sqrt(SimParams.Nbs*SimParams.Nms_v(1)/SimParams.L_v(1))*SimParams.Gains_vect;
    for t = 1:length(SimParams.Time)
        Gains_vect_t = reshape(SimParams.Gains_vect(:,t),SimParams.L_v(1),SimParams.U,SimParams.K);
        for u = 1:SimParams.U
                SimParams.u = u;
            [Hk,norm_factor] = gen_channel(SimParams,SimParams.theta_AoA(:,:,t),SimParams.theta_AoD(:,:,t),SimParams.phi_AoA(:,:,t),SimParams.phi_AoD(:,:,t),Gains_vect_t);
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
        
%         for u = 1:SimParams.U
%             SimParams.Hk_MU_current{u} = SimParams.Hk_MU{u,t};
%         end
%         SimParams.Flag_Use_Estimate = 0;
%         SimParams = SVDMMSE_v2 (SimParams);
%         [r] = achievableRate (SimParams.Ns,SimParams.U,SimParams.K,SimParams.var_n*eye(SimParams.Nms_v(1)),SimParams.Hk_MU,SimParams.P,SimParams.Q);
% 
%         [r] = achievableRate (SimParams.Ns,SimParams.U,SimParams.K,SimParams.var_n*eye(SimParams.Nms_v(1)),SimParams.Hk_MU_current,SimParams.P,SimParams.Q);
%         SimParams.SE_DL_PCSI(t,SimParams.nmc) = sum(r);
%         r_PCSI = sum(r)
    end
%     SNR
    SNR_av = 10*log10(mean(mean(SNR,1),2))
    %%%%%%% Rao - Blackwellized Particle Filter %%%%%%%%%
    SimParams.N = 1000;
    SimParams.nmc = nmc;
    SimParams = RBPF_Project_General_v3(SimParams);
    hat_theta_AoA = SimParams.hat_theta_AoA;
    hat_theta_AoD = SimParams.hat_theta_AoD;
    hat_phi_AoA = SimParams.hat_phi_AoA;
    hat_phi_AoD = SimParams.hat_phi_AoD;
    hat_Gains = SimParams.hat_Gains;
    P_hat = SimParams.P_hat_Gains;
    Gains_hat_MC(:,:,nmc) = hat_Gains;
    Gains_MC(:,:,nmc) = SimParams.Gains_vect;
    P_hat_MC(:,:,:,nmc) = P_hat;

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












