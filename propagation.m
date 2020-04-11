function [d_k,theta_k,phi_k,x_k,y_k,z_k] = propagation(d_prev,theta_prev,phi_prev,dt,vx,vy,vz)

x_k = d_prev*sin(phi_prev)*cos(theta_prev) + vx*dt;
y_k = d_prev*sin(phi_prev)*sin(theta_prev) + vy*dt;
z_k = d_prev*cos(phi_prev) + vz*dt;
d_k = sqrt(x_k^2 + y_k^2 + z_k^2);
theta_k = atan2(y_k,x_k);
if abs(theta_k - theta_prev) > pi/2
    theta_k = mod(theta_k,2*pi)-pi;
end
phi_k = atan2(sqrt(x_k^2 + y_k^2),z_k);
if abs(phi_k - phi_prev) > pi/4
    phi_k = mod(phi_k,pi)-pi/2;
end