function [zhat, H] = measurement_model_std(xe,idx)
%
% relative position meausrement model
%

global gDISTBEAR

fpos = 3+idx*2-1;

J = [0 -1; 1 0];
C = [cos(xe(3,1)) -sin(xe(3,1));
    sin(xe(3,1)) cos(xe(3,1))];

k_xL = C'*(xe(fpos:fpos+1)-xe(1:2));

if gDISTBEAR
    rho = norm(k_xL);
    th = atan2(k_xL(2),k_xL(1));
    zhat = [rho; th];
else
    zhat = k_xL;
end

H = zeros(2,size(xe,1));

if gDISTBEAR
    H_Lk = [ 1/norm(k_xL)*k_xL'; 1/norm(k_xL)^2*k_xL'*J' ];
    H(:,1:3) = - H_Lk* C'*[eye(2)  J*(xe(fpos:fpos+1,1)-xe(1:2,1))];
    H(:,fpos:fpos+1)= H_Lk* C';
    
else
    H(:,1:3) = - C'*[eye(2)  J*(xe(fpos:fpos+1,1)-xe(1:2,1))];
    H(:,fpos:fpos+1)= C';
end

