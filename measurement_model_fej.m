function [zhat, H] = measurement_model_fej(xe,idx,xL_1)


global gDISTBEAR


fpos = 3+idx*2-1;

C = [cos(xe(3)) -sin(xe(3));  sin(xe(3))  cos(xe(3)) ];

k_xL = C'*(xe(fpos:fpos+1)-xe(1:2));

if gDISTBEAR
    rho = norm(k_xL);
    th = atan2(k_xL(2),k_xL(1));
    zhat = [rho; th];
else
    zhat = k_xL;
end


% jacobian H evaluated at robot estimate of k|k-1 and landmark
% estimate of first time step

J = [0 -1; 1 0];

C = [cos(xe(3)) -sin(xe(3));   sin(xe(3)) cos(xe(3))];

H = zeros(2,size(xe,1));

k_xL = C'*(xL_1(idx*2-1:idx*2)-xe(1:2));

if gDISTBEAR
    H_Lk = [ 1/norm(k_xL)*k_xL'; 1/norm(k_xL)^2*k_xL'*J' ];    
    H(:,1:3) = - H_Lk* C'*[eye(2)  J*(xL_1(idx*2-1:idx*2)-xe(1:2,1))];
    H(:,fpos:fpos+1) = H_Lk* C';
else
    H(:,1:3) = -  C'*[eye(2)  J*(xL_1(idx*2-1:idx*2)-xe(1:2,1))];
    H(:,fpos:fpos+1) = C';
end
