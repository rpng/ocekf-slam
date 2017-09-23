function [zhat, H] = measurement_model_ocekf_1(xe,idx,V_1,xL_1,dpR,lambda_i)

global gDISTBEAR

fpos = 3+idx*2-1;
J = [0 -1; 1 0];
C = [cos(xe(3)) -sin(xe(3));  sin(xe(3))  cos(xe(3)) ];

k_xL = C'*(xe(fpos:fpos+1)-xe(1:2));

if gDISTBEAR    
    rho = norm(k_xL);
    th = atan2(k_xL(2),k_xL(1));
    zhat = [rho; th];
else
    zhat = k_xL;
end



% % Jacobian
H = zeros(2,size(xe,1));

pL_star = xe(fpos:fpos+1,1)-lambda_i/2;
k_xL = C'*(pL_star-xe(1:2));

if gDISTBEAR
    H_Lk = [ 1/norm(k_xL)*k_xL'; 1/norm(k_xL)^2*k_xL'*J' ];    
    H(:,1:3) = -  H_Lk*C'*[eye(2)  J*(pL_star-xe(1:2,1))];
    H(:,fpos:fpos+1) =  H_Lk*C';
else
    H(:,1:3) = - C'*[eye(2)  J*(pL_star-xe(1:2,1))];
    H(:,fpos:fpos+1) = C';
end

