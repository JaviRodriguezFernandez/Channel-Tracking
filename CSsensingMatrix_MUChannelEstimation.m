function [SimParams]=CSsensingMatrix_MUChannelEstimation(SimParams)

tsteps = SimParams.Mtck;
Nms_v = SimParams.Nms_v;
Nbs = SimParams.Nbs;
Lms_v = SimParams.Lms_v;
Lbs = SimParams.Lbs;
A = SimParams.A;
U = SimParams.U;

%Nt=64; Nr=16; tsteps=2*64*16; Lt=1; Lr=1; A='A1A2';

% input
% tsteps: number of training steps
% Nt: number of Tx antennas
% Nr: number or Rx antennas
% Lt: number of Tx RF chains
% Lr: number of Rx RF chains
% A: architecture, e.g. 'A1A2' A1 at Bs and A2 at MS
% method: random, deterministic

% tsteps=mt*mr; snapshots
% measurements = mt*mr*Lr 
% tsteps=mr*mt;
Nbits = 4;
Nres = 2^Nbits;

switch (A)
    case 'A1A1'    % phase shifters     
        
        % Generate pseudorandom precoding matrices for every user
        W = [];
        F = [];
        Phi = [];
        W_2 = [];
        for m = 1:tsteps
            f_m = [];
            for u = 1:U
                tt_u=randi(Nres,[Nms_v(u) Lms_v(u)]);
                q_u = sqrt(1/2/Lms_v(u)/U)*(sign(randn(Lms_v(u),1))+1i*sign(randn(Lms_v(u),1)));
                tt_u = reshape(tt_u,[],1);
                for i = 1:Nres
                    tt_u(tt_u==i) = 1/sqrt(Nms_v(u))*exp(1i*2*pi*(i-1)/Nres);
                end
                
                tt_u = reshape(tt_u,[Nms_v(u) Lms_v(u)]);
                f_u = tt_u*q_u;
                f_u = f_u/norm(f_u)*sqrt(size(f_u,2));
                f_m = [f_m; f_u/sqrt(U)];
            end
            F = [F f_m];
            tt=randi(Nres,[Nbs Lbs]);
            tt = reshape(tt,[],1);
            for i = 1:Nres
                tt(tt==i) = 1/sqrt(Nbs)*exp(1i*2*pi*(i-1)/Nres);
            end
            tt = reshape(tt,[Nbs Lbs]);
            Phi = [Phi; kron(f_m.',tt')];
            W = [W tt];
            W_2 = blkdiag(W_2,tt);
        end
        SimParams.Phi = Phi;
        SimParams.W = W;
        SimParams.W_2 = W_2;
        SimParams.F = F;
        
    case 'A1A2'         
        
        tt=randi(16,[Nt tsteps]);
        for i = 1:15
            tt(tt==i) = exp(1i*2*pi*(i-1)/16);
        end
        P=tt;  
        P=P/sqrt(Nt);                          
        
        Q=INr(1:Nr,randi(Nr,[tsteps*Lr,1])); 
        
        for i=1:tsteps,
            phi((i-1)*Lr+(1:Lr),:)=kron(P(:,i).',Q(:,(i-1)*Lr+(1:Lr))');
        end            
        
    case 'A1A3'   % phase shifters to a subset of antennas
        tt=randi(4,[Nt tsteps]);
        tt(tt==2)=-1i;
        tt(tt==3)=1i;
        tt(tt==4)=-1;
        P=tt;  
        P=P/sqrt(Nt);                          
        
        Q=[];
        for i=1:Lr,
            tt=randi(4,[Mr tsteps]);
            tt(tt==2)=-1i;
            tt(tt==3)=1i;
            tt(tt==4)=-1;
            Qf=tt;                      
            Q=blkdiag(Q,Qf);
        end                
        
        for i=1:tsteps,
            phi((i-1)*Lr+(1:Lr),:)=kron(P(:,i).',Q(:,(i-1)*Lr+(1:Lr))');
        end
        
     case 'A1A4' % switches combined signal from a subset of antennas         
        tt=randi(4,[Nt tsteps]);
        tt(tt==2)=-1i;
        tt(tt==3)=1i;
        tt(tt==4)=-1;
        P=tt;  
        P=P/sqrt(Nt);                          
        
        Q=[];
        for i=1:Lr,
            Qf=IMr(1:Mr,randi(Mr,[tsteps,1])); 
            Q=blkdiag(Q,Qf);
        end                
        
        for i=1:tsteps,
            phi((i-1)*Lr+(1:Lr),:)=kron(P(:,i).',Q(:,(i-1)*Lr+(1:Lr))');
        end
        
end

