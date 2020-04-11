function SimParams = Frequency_Domain_MU_Channel_Estimation_UL(SimParams)

%%%%%%%% Multiuser UL Wideband Channel Estimation exploiting common suppport

Nbs = SimParams.Nbs;
Nms_v = SimParams.Nms_v;
Lbs = SimParams.Lbs;
Lms_v = SimParams.Lms_v;
Gbs = SimParams.Gbs;
Gms_v = SimParams.Gms_v;
tsteps = SimParams.tsteps;
K = SimParams.K;
Q = SimParams.W;
C_noise = zeros(tsteps*Lbs,tsteps*Lbs);
var_n = SimParams.var_n;
Phi = SimParams.Phi;
Hk = SimParams.Hk;
U = SimParams.U;

n_c = zeros(tsteps*Lbs,K);

for k = 1:K
    for m = 1:tsteps
        noise = sqrt(var_n/2)*(randn(Nbs,1)+1i*randn(Nbs,1));
        n_c((m-1)*Lbs+(1:Lbs),k) = Q(:,(m-1)*Lbs+(1:Lbs))'*noise;
        C_noise((m-1)*Lbs+(1:Lbs),(m-1)*Lbs+(1:Lbs)) = Q(:,(m-1)*Lbs+(1:Lbs))'*Q(:,(m-1)*Lbs+(1:Lbs));
    end
end

for k = 1:K
    H_k = [];
    for u = 1:U
        H_dummy = Hk{u};
        H_k = [H_k H_dummy(:,:,k)];
    end
    y(:,k) = Phi*reshape(H_k,[],1);
    r(:,k) = y(:,k) + n_c(:,k);
    SNR(k) = y(:,k)'*y(:,k)/(n_c(:,k)'*n_c(:,k));
end

SNR_av = 10*log10(mean(SNR))
Dw = chol(C_noise)';
SimParams.Dw = Dw;
SimParams.r = r;

SimParams = SC_SS_SW_IPM(SimParams);

for u = 1:SimParams.U
    Hk_MU_est_u = zeros(Nbs,Nms_v(u),K);
    Psi_u = SimParams.hatPsi{u};
    for k = 1:K
        Hk_MU_est_u(:,:,k) = reshape(Psi_u*SimParams.hatG(u,k),[Nbs,Nms_v(u)]);
        SimParams.Hk_MU_est{u} = Hk_MU_est_u;
    end
end

num = 0;
den = 0;
for u = 1:SimParams.U
    num = num + sum(sum(sum(abs(SimParams.Hk_MU_est{u} - Hk{u}).^2)));
    den = den + sum(sum(sum(abs(Hk{u}).^2)));
end

NMSE_Freq_SWOMP_UL = num/den;

SimParams.NMSE_Freq_SWOMP_UL(SimParams.slot) = NMSE_Freq_SWOMP_UL;
% NMSE_Freq_SWOMP_UL = sum(sum(sum(sum(abs(Hk_MU_est(:,:,:,:) - Hk_MU(:,:,:,:)).^2))))/sum(sum(sum(sum(abs(Hk_MU(:,:,:,:)).^2))));

