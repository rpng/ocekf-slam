function [zhat, H] = measurement_model_id(xe,idx,idm,xR_true,xL_true)
% 
% relative position meausrement model
% 
% idx: landmark id in state vector
% idm: landmark id in map
% 

global gDISTBEAR

fpos = 3+idx*2-1;

C= [cos(xe(3)) -sin(xe(3));  sin(xe(3)) cos(xe(3)) ];

k_xL = C'*(xe(fpos:fpos+1)-xe(1:2));

if gDISTBEAR
    rho = norm(k_xL);
    th = atan2(k_xL(2),k_xL(1));
    zhat = [rho; th];
else
    zhat = k_xL;
end


% % true jacobians
J= [0 -1; 1 0];
C= [cos(xR_true(3)) -sin(xR_true(3));   sin(xR_true(3)) cos(xR_true(3))];

H = zeros(2,size(xe,1));

k_xL = C'*(xL_true(:,idm)-xR_true(1:2,1));


if gDISTBEAR
    H_Lk = [ 1/norm(k_xL)*k_xL'; 1/norm(k_xL)^2*k_xL'*J' ];    
    H(:,1:3) = - H_Lk* C'*[eye(2)  J*(xL_true(:,idm)-xR_true(1:2,1))];
    H(:,fpos:fpos+1)= H_Lk* C';
else
    H(:,1:3) = - C'*[eye(2)  J*(xL_true(:,idm)-xR_true(1:2,1))];
    H(:,fpos:fpos+1)= C';
end
