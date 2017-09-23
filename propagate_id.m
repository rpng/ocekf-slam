function [xe,Pe,PHI_id,G_id] = propagate_id(xe,Pe,dt,v_m,omega_m,sigma_v,sigma_w,xR_true,v_true,omega_true)


J = [0 -1; 1 0];
Q = [sigma_v^2  0;   0    sigma_w^2];

xR_k_k = xe(1:3,1);


% % propagate state
xe(1:3,1) = [ xe(1) + v_m*dt*cos(xe(3));
    xe(2) + v_m*dt*sin(xe(3));
    pi_to_pi(xe(3) + omega_m*dt) ];

% % jacobians evaluated at true values of the state
% % 
xR_true_2 = [ xR_true(1) + v_true*dt*cos(xR_true(3));
        xR_true(2) + v_true*dt*sin(xR_true(3));
        pi_to_pi(xR_true(3) + omega_true*dt) ];


PHI = [ eye(2)   J*(xR_true_2(1:2,1)-xR_true(1:2,1));
    zeros(1,2)  1 ];



G = [dt*cos(xR_true(3)) 0;
    dt*sin(xR_true(3)) 0;
    0    dt ];


Qprime = G*Q*G';
Qprime = compute_Q(xR_true,dt,v_true,omega_true,sigma_v,sigma_w);


% % propagate covariance
Pe(1:3,1:3) = PHI*Pe(1:3,1:3)*PHI' + Qprime;
if size(Pe,1)>3
    Pe(1:3,4:end) = PHI*Pe(1:3,4:end);
    Pe(4:end,1:3) = Pe(1:3,4:end)';
end
Pe = 0.5*(Pe+Pe');


PHI_id = blkdiag(PHI,eye(size(xe,1)-3)); %PHI%
G_id = [G; zeros(size(xe,1)-3,2)]; %G

