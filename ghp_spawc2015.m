function [Frf,Fbb]=ghp_spawc2015(Y,p,Qbits)
% "Dictionary-free hybrid precoders and combiners for mmWave MIMO systems"
% Spawc2015
% input:
%  Y  : Fopt optimum precoder
%  p  : number of Rf chains
%------------------------
% size(Y) = [m n]
% n num of Streams
% m num Antennas
% Fopt  = Frf Fbb
% Y = AB

deltaQ = 2*pi/2^Qbits;
% 1) Precoder
[m n] = size(Y);
A = zeros(m,p);
B = zeros(p,n);

R=Y;
Res=R;
A(:,1:n)=Y(:,1:n)./abs(Y(:,1:n));    
for k=1:p-n,
    B(1:k+n-1,:)=A(:,1:k+n-1)\Y;           
    R=Res-A*B;    
    [u s v] = svd(R);    
    A(:,k+n)=u(:,1)./abs(u(:,1));            
end
AngleM = angle(A);
AngleQ = round(AngleM/deltaQ)*deltaQ;
AngleQ = AngleM;
A = exp(1i*AngleQ);
B=A\Y;

Fbb=B;
Frf=A;
% normalize power constraint ||Fopt||^2_F=NS
Cs=sqrt(n)/norm((Frf*Fbb),'fro');
Fbb=Fbb*Cs;