function [Q] = compute_Q(x,dt,v,omega,sigma_v,sigma_w)

x1k = x(1:3,1);
P1k = zeros(3);

nddt = 3;
ddt = dt/nddt;

J = [0 -1; 1 0];
Q = [sigma_v^2  0;  0  sigma_w^2]  / ddt; %quatanion tech report, eq.101


for i=1:nddt
    
    xk = [ x1k(1) + v*ddt*cos(x1k(3));
        x1k(2) + v*ddt*sin(x1k(3));
        pi_to_pi(x1k(3) + omega*ddt) ];
    
    F1k = [1 0 -v*ddt*sin(x1k(3));
        0 1  v*ddt*cos(x1k(3));
        0 0   1];

    
    Gk = [ddt*cos(x1k(3))   0;
        ddt*sin(x1k(3))   0;
        0     ddt];
    

    P1k = F1k*P1k*F1k'+Gk*Q*Gk'; %Lyapunov eq.

    x1k = xk;
    
end

Q = P1k;
