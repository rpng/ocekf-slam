function [xe, Pe, PHI,G,PHI_mult, xR_k_k1,dpR_star_prev,pR_star_prev,lambda] = ...
    propagate_ocekf(xe,Pe,dt,v_m,omega_m,sigma_v,sigma_w, PHI_mult, xR_k_k1,dpR_star_prev,pR_star_prev,xL_1,lambda, lm_seq,z)

J = [0 -1; 1 0];

xR_k_k = xe(1:3,1);


% propagate state
xe(1:3,1) = [ xe(1) + v_m*dt*cos(xe(3,1));
    xe(2) + v_m*dt*sin(xe(3,1));
    pi_to_pi(xe(3) + omega_m*dt) ];


% compute lagrangian multiplier: lambda
observed_land_id = [];
lenz = length(find(z(3,:)>0));
for i= 1:lenz
    %data association (based on landmark id)
    is_exist = ~(lm_seq - z(3,i));
    idx = find(is_exist);
    if ~isempty(idx)
        observed_land_id = [observed_land_id;idx];
    end
end
nland = length(observed_land_id);

if nland>0
    A = kron(ones(nland)+eye(nland), eye(2));
    for i=1:nland
        ii = 2*i+(-1:0);
        idx = observed_land_id(i);
        fpos = 3+idx*2-1;
        lpos = idx*2-1;
        b(ii,1) = 2*( (xe(fpos:fpos+1,1)-xL_1(lpos:lpos+1,1)) - (xR_k_k(1:2,1)-xR_k_k1(1:2,1)-dpR_star_prev) );
    end
    lambda = A\b;
    lambda = reshape(lambda, 2,nland);
        
    % linearization point for pR
    pR_star = xR_k_k(1:2,1) + sum(lambda,2)/2;
else
    lambda = 0;
    pR_star = xR_k_k(1:2,1);
end

PHI_R = [ eye(2)   J*(xe(1:2,1)-pR_star);  
    zeros(1,2)  1 ];


G_R = [dt*cos(xR_k_k(3,1))   0;
    dt*sin(xR_k_k(3,1))   0;
    0     dt];


% odometry noise cov
Q = [sigma_v^2  0;   0    sigma_w^2];

Qprime = G_R*Q*G_R';


% propagate covariance
Pe(1:3,1:3) = PHI_R*Pe(1:3,1:3)*PHI_R' + Qprime;
if size(Pe,1)>3
    Pe(1:3,4:end) = PHI_R*Pe(1:3,4:end);
    Pe(4:end,1:3) = Pe(1:3,4:end)';
end
Pe = 0.5*(Pe+Pe');


% %
PHI = blkdiag(PHI_R,eye(size(xe,1)-3));
G = [G_R; zeros(size(xe,1)-3,2)];


% % % PHI_mult=PHI_R(k+m-1)*...*PHI_R(k+1)*PHI_R(k), used in computing nullspace of H
PHI_mult = PHI_R*PHI_mult;


%%%%%
dpR_star_prev = dpR_star_prev + xR_k_k1(1:2,1)-pR_star ;
pR_star_prev = pR_star;
xR_k_k1 = xe(1:3,1);

