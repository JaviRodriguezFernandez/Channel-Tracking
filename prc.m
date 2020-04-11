function val = prc(SimParams,tau,d)

rolloff = 0.8;
Ts = SimParams.Ts;

val = sinc((d*Ts-tau)/Mfilter/Ts)*cos(pi*rolloff*(d*Ts-tau)/Mfilter/Ts)/(1-(2*rolloff*(d*Ts-tau)/Mfilter/Ts)^2);