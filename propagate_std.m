function [xe, Pe, PHI_std,G_std] = propagate_std(xe,Pe,dt,v_m,omega_m,sigma_v,sigma_w)


J = [0 -1; 1 0];
Q = [sigma_v^2  0;   0    sigma_w^2];

xR_k_k = xe(1:3,1);


% % propagate state
xe(1:3,1) = [ xe(1) + v_m*dt*cos(xe(3,1));
    xe(2) + v_m*dt*sin(xe(3,1));
    pi_to_pi(xe(3) + omega_m*dt) ];


% % jacobians: evaluated at k+1|k, k|k
PHI = [ eye(2)   J*(xe(1:2,1)-xR_k_k(1:2,1));
    zeros(1,2)  1 ];


G = [dt*cos(xR_k_k(3,1))   0;
    dt*sin(xR_k_k(3,1))   0;
    0     dt];


Qprime = G*Q*G';
Qprime = compute_Q(xR_k_k,dt,v_m,omega_m,sigma_v,sigma_w);


% propagate covariance
Pe(1:3,1:3) = PHI*Pe(1:3,1:3)*PHI' + Qprime;
if size(Pe,1)>3
    Pe(1:3,4:end) = PHI*Pe(1:3,4:end);
    Pe(4:end,1:3) = Pe(1:3,4:end)';
end
Pe = 0.5*(Pe+Pe');

% %
PHI_std = blkdiag(PHI,eye(size(xe,1)-3));
G_std = [G; zeros(size(xe,1)-3,2)];

