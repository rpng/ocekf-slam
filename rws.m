function [v_m,omega_m,v,omega,xR_true,z,R] = rws(timesteps,dt,v_true,omega_true,sigma_v,sigma_w,sigma_r,sigma_th,sigma_p,xL_true,max_range,min_range)
% % REAL WORLD SLAM SIMULATION

% %true odometry
v = v_true*ones(1,timesteps);
omega = omega_true*ones(1,timesteps);  %loop
% omega = ((rand(1,timesteps)-.5)) *.25;  %random

% %odometry measurements
v_m = v + sigma_v*randn(size(v)) ;
omega_m = omega + sigma_w*randn(size(omega)) ;


% %
nL = size(xL_true,2);
z = zeros(3,nL,timesteps);

for k= 1:timesteps
    R{k}=[];
end


for k= 1:timesteps
    
    % % constant velocity motion model
    if k==1
        xR_true(:,k) = zeros(3,1);
    else
        xR_true(1,k) = xR_true(1,k-1) + v(k-1)*dt*cos(xR_true(3,k-1));
        xR_true(2,k) = xR_true(2,k-1) + v(k-1)*dt*sin(xR_true(3,k-1));
        xR_true(3,k) = xR_true(3,k-1) + omega(k-1)*dt;
    end
    
    curr_meas_num = 0;
    
    for i= 1:nL
        % measurement wrt *global* frame
        [th,r] = cart2pol(xL_true(1,i)-xR_true(1,k), xL_true(2,i)-xR_true(2,k));
        %measurement wrt robot
        th = pi_to_pi(th-xR_true(3,k));
        sigma_r = sigma_p*r; %sigma_p percentage of range
        %use measurement only if landmark is closer than max_range
        if r<max_range && r>min_range
            curr_meas_num = curr_meas_num+1;
            %store measurement, and landmark id
            global gDISTBEAR
            if gDISTBEAR==1 || gDISTBEAR==-1
                %distance-bearing meausement
                Rii = diag([sigma_r^2,sigma_th^2]);
                R{k} = blkdiag(R{k},Rii);
                noise = mvnrnd([0;0],Rii);
                r = r + sigma_r*randn; %noise(1);%
                th = th + sigma_th*randn; %noise(2);%
                
                z(1:2,curr_meas_num,k) = [r;th];%
            else
                %relative position measurement
                [dx,dy] = pol2cart(th,r);
                %add noise
                sig_d= sigma_p*r; %sigma_p percentage of range
                sig_d = sigma_p; % or constant std
                sig_o= 0;
                Rii= [sig_d^2  sig_o^2; sig_o^2  sig_d^2];
                R{k} = blkdiag(R{k},Rii);
                noise = mvnrnd([0;0],Rii);
                dx = dx + noise(1); %sig_d*randn; %
                dy = dy + noise(2); %sig_d*randn; %
                z(1:2,curr_meas_num,k) = [dx;dy];
            end
            
            z(3,curr_meas_num,k) = i;
            
        end
        
    end%nL
    
end%timesteps

