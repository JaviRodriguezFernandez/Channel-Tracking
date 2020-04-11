function [Hk,norm_factor] = gen_channel_ULA(SimParams,theta_AoA,theta_AoD,Gains)

u = SimParams.u;
Nbs = SimParams.Nbs;
Nms_v = SimParams.Nms_v;
K = SimParams.K;
Lu = size(theta_AoA,2);

% Generation of the multi-user dictionary matrix
SimParams.Psi = [];

epsilon_azi_BS=[0:Nbs-1]';
epsilon_azi_MS=[0:Nms_v(u)-1]';

SimParams.Abs_azi{u} = ([exp(-1i*pi*cos(theta_AoA(u,:).')*epsilon_azi_BS')]/sqrt(Nbs)).';
SimParams.Ams_azi{u} = ([exp(-1i*pi*cos(theta_AoD(u,:).')*epsilon_azi_MS')]/sqrt(Nms_v(u))).';

for l = 1:size(SimParams.Abs_azi{u},2)
    Ams(:,l) = SimParams.Ams_azi{u}(:,l);
    Abs(:,l) = SimParams.Abs_azi{u}(:,l);
end

% for i=1:Lu
%     Ar_azi(:,i) = [exp(-1i*pi*cos(AoA_azi_u(i))*epsilon_azi_BS)]/sqrt(Nbsx);
%     At_azi(:,i) = [exp(-1i*pi*cos(AoD_azi_u(i))*epsilon_azi_MS)]/sqrt(Nmsux);
%     Ar_ele(:,i) = [exp(-1i*pi*sin(AoA_ele_u(i))*sin(AoA_azi_u(i))*epsilon_ele_BS)]/sqrt(Nbsy);
%     At_ele(:,i) = [exp(-1i*pi*sin(AoD_ele_u(i))*sin(AoD_azi_u(i))*epsilon_ele_MS)]/sqrt(Nmsuy);
%     Ar(:,i) = kron(Ar_azi(:,i),Ar_ele(:,i));
%     At(:,i) = kron(At_azi(:,i),At_ele(:,i));
%     Psi_u(:,i) = kron(conj(At(:,i)),Ar(:,i));
% end

normHk = zeros(K,1);

Hk = zeros(Nbs,Nms_v(u),K);
for k = 1:K
    Hv(:,:,k) = diag(squeeze(Gains(:,u,k)));
    Hk(:,:,k) = Abs*Hv(:,:,k)*Ams';
    normHk(k) = norm(Hk(:,:,k),'fro')^2;
end

rho = Nbs*Nms_v(1)/mean(normHk);
norm_factor = rho;
Hk = Hk*sqrt(rho);

% Hk = sqrt(SimParams.Nbs*SimParams.Nms_v(u)/SimParams.L_v(u))*Hk;





