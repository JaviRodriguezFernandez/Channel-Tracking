function SimParams=CSsensingMatrix(SimParams);
%Nt=64; Nr=16; tsteps=2*64*16; Lt=1; Lr=1; A='A1A2';

% input parameters
% M: number of training steps
% Nt: number of Tx antennas
% Nr: number or Rx antennas
% Lt: number of Tx RF chains
% Lr: number of Rx RF chains
% A: architecture, e.g. 'A1A2' A1 at Bs and A2 at MS
% method: random, deterministic

% tsteps=mt*mr; snapshots
% measurements = mt*mr*Lr 

Nt = SimParams.Nms_v;
Nr = SimParams.Nbs;
Lt = SimParams.Lms_v;
Lr = SimParams.Lbs;
A = SimParams.A;
Ns = SimParams.Ns;
M = SimParams.M_ini;


Nbits = 4;
Nres = 2^Nbits;

switch (A)
    case 'A1A1'    % Fully-connected network of phase-shifters at both Tx and Rx    
        tt=randi(Nres,[Nt M*Lt]);
        for i = 1:Nres
            tt(tt==i) = exp(1i*2*pi*(i-1)/Nres);
        end
        P=tt;  
        
        tt=randi(Nres,[Nr M*Lr]);
        for i = 1:Nres
            tt(tt==i) = exp(1i*2*pi*(i-1)/Nres);
        end
        Q=tt;    
        P = P/sqrt(Nt);
        Q = Q/sqrt(Nr);
        
        for i=1:M,
            signal = sqrt(1/2/Lt)*(sign(randn(Lt,1))+1i*sign(randn(Lt,1)));
            SimParams.signal(:,i) = signal; 
            phi((i-1)*Lr+(1:Lr),:)=kron(signal.'*P(:,(i-1)*Lt+(1:Lt)).',Q(:,(i-1)*Lr+(1:Lr))');
        end
end

[tsteps  ~]=size(phi);
tsteps=tsteps/Lr;

SimParams.Phi = phi;
SimParams.P = P;
SimParams.Q = Q;
