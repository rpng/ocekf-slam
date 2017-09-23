%% -----------------------------------------------------------------------------------------------
% %  Simulation script: 2D OC-EKF SLAM SIMULATION
% %  Reference:
% %  - Guoquan Huang, Anastasios I. Mourikis and Stergios I. Roumeliotis.
% %    Observability-based Rules for Designing Consistent EKF SLAM Estimators.
% %    International Journal of Robotics Research, vol. 29, no. 5, pp. 502-528, April 2010.
% %
% %  Copyright (C) Robot Perception and Navigation Group (RPNG) - All Rights Reserved
% %  Unauthorized copying of this file, via any medium is strictly prohibited
% %  Proprietary and confidential material are included.
% %
% %  Written by:
% %  Guoquan (Paul) Huang <ghuang@udel.edu>
%% -----------------------------------------------------------------------------------------------


clear all
close all
clc


%%
randomselect='repeatlast'; % 'random' or 'repeatlast';
switch randomselect
    case 'repeatlast'
        load randstates;
        rand('state',randstate);
        randn('state',randnstate);
    case 'random'
        % do nothing
    otherwise
        rand('state',randomselect);
        randn('state',randomselect);
end
% %
randstate=rand('state');
randnstate=randn('state');
save randstates.mat randstate randnstate;


%% Simulation Parameters

global gDISTBEAR
gDISTBEAR = 0; %1 - range and bearing;  0 - relative position;

dt = 1;
v_true = .25;
omega_true = .025;
sigma = .05*v_true;
sigma_v = sigma/sqrt(2); %.1*v_true; %
sigma_w = 2*sqrt(2)*sigma; %.1*omega_true; 1*pi/180; %
Q = diag([sigma_v^2 sigma_w^2]);

sigma_p = .1; %noise is the percentagae of distance measurement, BUT double check rws.m since sometimes we use this as a constant absolute sigma (like in invariant EKF tro submission)
sigma_r = 1; %range measmnt noise
sigma_th = 10*pi/180;%bearing measuremnt noise

nL = 20; %number of landmarks
nSteps = 2500; %nubmer of time steps
nRuns = 5; %number of monte carlo runs
if nL==1
    max_range = 200;%always observe this landmark
else
    max_range = 5;%env_size/10;%.5*v_true/omega_true;%
end
min_range = .5;
init_steps = 0;%3
max_delay = 10;%for delayed initial
NIncr = 0; %increment number of incremental MAP: 0 - not runing

global gNPERIODIC gREORDERING gRELINEARIZATION
gNPERIODIC = 100; %in isam, periodic batch update..
gREORDERING = 1; %isam- reordering
gRELINEARIZATION = 1; %iSAM: relinearization and redo qr: check if it is matlab qr or davis' cs_qr



%% preallocate memory for saving resutls

% Ideal EKF
xRest_id = zeros(3,nSteps,nRuns); %estimated traj
xRerr_id = zeros(3,nSteps,nRuns); %all err state
Prr_id = zeros(3,nSteps,nRuns); %actually diag of Prr
neesR_id = zeros(1,nSteps,nRuns); %nees (or mahalanobis distance)
rmsRp_id =  zeros(1,nSteps,nRuns); %rms of robot position
rmsRth_id = zeros(1,nSteps,nRuns); %rms of robot orientation
xLest_id = zeros(2,nL,nSteps,nRuns);
xLerr_id = zeros(2,nL,nSteps,nRuns);
Pll_id = zeros(2,nL,nSteps,nRuns);
neesL_id = zeros(1,nL,nSteps,nRuns);
rmsL_id = zeros(1,nL,nSteps,nRuns);
nees_id = zeros(1,nSteps,nRuns); %nees for the whole state

% Standard EKF
xRest_std = zeros(3,nSteps,nRuns); %estimated traj
xRerr_std = zeros(3,nSteps,nRuns); %all err state
Prr_std = zeros(3,nSteps,nRuns); %actually diag of Prr
neesR_std = zeros(1,nSteps,nRuns); %nees (or mahalanobis distance)
rmsRp_std =  zeros(1,nSteps,nRuns); %rms of robot position
rmsRth_std = zeros(1,nSteps,nRuns); %rms of robot orientation
xLest_std = zeros(2,nL,nSteps,nRuns);
xLerr_std = zeros(2,nL,nSteps,nRuns);
Pll_std = zeros(2,nL,nSteps,nRuns);
neesL_std = zeros(1,nL,nSteps,nRuns);
rmsL_std = zeros(1,nL,nSteps,nRuns);
nees_std = zeros(1,nSteps,nRuns); %nees for the whole state
kld_std = zeros(1,nSteps,nRuns); % KLD

% FEJ-EKF
xRest_fej = zeros(3,nSteps,nRuns); %estimated traj
xRerr_fej = zeros(3,nSteps,nRuns); %all err state
Prr_fej = zeros(3,nSteps,nRuns); %actually diag of Prr
neesR_fej = zeros(1,nSteps,nRuns); %nees (or mahalanobis distance)
rmsRp_fej =  zeros(1,nSteps,nRuns); %rms of robot position
rmsRth_fej = zeros(1,nSteps,nRuns); %rms of robot orientation
xLest_fej = zeros(2,nL,nSteps,nRuns);
xLerr_fej = zeros(2,nL,nSteps,nRuns);
Pll_fej = zeros(2,nL,nSteps,nRuns);
neesL_fej = zeros(1,nL,nSteps,nRuns);
rmsL_fej = zeros(1,nL,nSteps,nRuns);
nees_fej = zeros(1,nSteps,nRuns); %nees for the whole state
kld_fej = zeros(1,nSteps,nRuns); % KLD

% OC-EKF
xRest_ocekf_1 = zeros(3,nSteps,nRuns); %estimated traj
xRerr_ocekf_1 = zeros(3,nSteps,nRuns); %all err state
Prr_ocekf_1 = zeros(3,nSteps,nRuns); %actually diag of Prr
neesR_ocekf_1 = zeros(1,nSteps,nRuns); %nees (or mahalanobis distance)
rmsRp_ocekf_1 =  zeros(1,nSteps,nRuns); %rms of robot position
rmsRth_ocekf_1 = zeros(1,nSteps,nRuns); %rms of robot orientation
xLest_ocekf_1 = zeros(2,nL,nSteps,nRuns);
xLerr_ocekf_1 = zeros(2,nL,nSteps,nRuns);
Pll_ocekf_1 = zeros(2,nL,nSteps,nRuns);
neesL_ocekf_1 = zeros(1,nL,nSteps,nRuns);
rmsL_ocekf_1 = zeros(1,nL,nSteps,nRuns);
nees_ocekf_1 = zeros(1,nSteps,nRuns); %nees for the whole state
kld_ocekf_1 = zeros(1,nSteps,nRuns); % KLD


%% LANDMARK GENERATION: same landmarks in each run
if nL==1
    xL_true_fixed = [0;v_true/omega_true];
elseif nL==2
    xL_true_fixed = [5 -5; 5 15];
else
    xL_true_fixed = gen_map(nL,v_true,omega_true,min_range, max_range, nSteps,dt);%max_range=5
end


%% Monte Carlo Simulations

tic
for kk = 1:nRuns
    
    kk
    
    % % real world simulation data % %
    xL_true(:,:,kk) = xL_true_fixed;
    [v_m,omega_m, v_true_all,omega_true_all, xR_true(:,:,kk), z,R] = rws(nSteps, dt,v_true,omega_true,sigma_v,sigma_w,sigma_r,sigma_th,sigma_p,xL_true(:,:,kk),max_range,min_range);
    
    
    % % INITIALIZATION
    x0 = zeros(3,1);
    if init_steps
        P0 = zeros(3);
    else
        P0 = diag([(.0001)^2,(.0001)^2,(.0001)^2]);
    end
    
    for k=1:init_steps
        [x0,P0] = propagate_std(x0,P0,dt,v_m(k),omega_m(k),sigma_v,sigma_w);
    end
    
    
    % Ideal EKF
    xe_id = x0;
    Pe_id = P0;
    V_id = [];
    
    % Standard EKF
    xe_std = x0;
    Pe_std = P0;
    V_std = [];
    
    
    % FEJ-EKF
    xe_fej = x0;
    Pe_fej = P0;
    xL_fej_1 = [];
    xR_fej_k_k1 = xe_fej(1:3,1);
    PHI_mult_fej = eye(1); %propagation jacobian product
    V_fej = [];
    
    % OC-EKF
    xe_ocekf_1 = x0;
    Pe_ocekf_1 = P0;
    xL_ocekf_1_1 = [];
    xR_oc_k_k1_1 = xe_ocekf_1(1:3,1);
    dpR_star_prev_1 = zeros(2,1);
    pR_star_prev = x0(1:2,1);
    dpR_ocekf_1 = zeros(2,1);
    V_ocekf_1 = [];
    PHI_mult_ocekf_1 = eye(1);
    lambda_1 = zeros(2,nL);
    
    
    
    % list of landmark ids that sequentially appear in the state vector
    lm_seq_id = [];
    lm_seq_std = [];
    lm_seq_fej = [];
    lm_seq_ocekf_1 = [];
    
    
    
    for k= init_steps+1:nSteps-1
        %first init_steps for ekf propagation to produce nonzero init cov
        
        timestep = k+1
        
        
        % % PROPAGATE: k+1|k
        [xe_id,Pe_id,PHI_id,G_id] = propagate_id(xe_id,Pe_id,dt,v_m(k),omega_m(k),sigma_v,sigma_w,xR_true(:,k,kk),v_true_all(k),omega_true_all(k));
        [xe_std,Pe_std,PHI_std,G_std] = propagate_std(xe_std,Pe_std,dt,v_m(k),omega_m(k),sigma_v,sigma_w);
        [xe_fej,Pe_fej,xR_fej_k_k1,PHI_fej,G_fej] = propagate_fej(xe_fej,Pe_fej,dt,v_m(k),omega_m(k),sigma_v,sigma_w,xR_fej_k_k1, PHI_mult_fej);
        [xe_ocekf_1,Pe_ocekf_1,PHI_ocekf_1,G_ocekf_1,PHI_mult_ocekf_1,  xR_oc_k_k1_1,dpR_star_prev_1,pR_star_prev,lambda_1] = propagate_ocekf_1(xe_ocekf_1,Pe_ocekf_1,dt,v_m(k),omega_m(k),sigma_v,sigma_w,PHI_mult_ocekf_1, xR_oc_k_k1_1,dpR_star_prev_1,pR_star_prev,xL_ocekf_1_1,lambda_1, lm_seq_ocekf_1,z(:,:,k+1));
        
        
        % % UPDATE: k+1|k+1
        [xe_id,Pe_id,lm_seq_id,V_id] = update_id(xe_id,Pe_id,lm_seq_id,z(:,:,k+1),R{k+1},xR_true(:,k+1,kk),xL_true(:,:,kk), V_id,PHI_id);
        [xe_std,Pe_std,lm_seq_std] = update_std(xe_std,Pe_std,lm_seq_std,z(:,:,k+1),R{k+1});
        [xe_fej,Pe_fej,xL_fej_1,lm_seq_fej,V_fej] = update_fej(xe_fej,Pe_fej,xL_fej_1,lm_seq_fej,z(:,:,k+1),R{k+1}, V_fej);
        [xe_ocekf_1,Pe_ocekf_1,xL_ocekf_1_1,lm_seq_ocekf_1, V_ocekf_1,dpR_ocekf_1,lambda_1] = update_ocekf_1(xe_ocekf_1,Pe_ocekf_1,xL_ocekf_1_1,lm_seq_ocekf_1,z(:,:,k+1),R{k+1}, PHI_mult_ocekf_1,V_ocekf_1,dpR_ocekf_1, lambda_1, dpR_star_prev_1);
        
        
        % % SAVE RESULTS        
        % Ideal EKF
        xRest_id(:,k+1,kk) = xe_id(1:3,1);
        err = xR_true(1:3,k+1,kk) - xe_id(1:3,1);
        err(3) = pi_to_pi(err(3));
        xRerr_id(:,k+1,kk) = err;
        Prr_id(:,k+1,kk) = diag(Pe_id(1:3,1:3));
        neesR_id(:,k+1,kk) = err'*inv(Pe_id(1:3,1:3))*err;
        rmsRp_id(:,k+1,kk) = err(1:2,1)'*err(1:2,1);
        rmsRth_id(:,k+1,kk) = err(3,1)'*err(3,1);
        
        for i=1:length(lm_seq_id)
            xLest_id(:,lm_seq_id(i),k+1,kk) = xe_id(3+2*i-1:3+2*i,1);
            err = xL_true(:,lm_seq_id(i),kk) - xLest_id(:,lm_seq_id(i),k+1,kk);
            xLerr_id(:,lm_seq_id(i),k+1,kk) = err;
            Pll_id(:,lm_seq_id(i),k+1,kk) = diag(Pe_id(3+2*i-1:3+2*i,3+2*i-1:3+2*i));
            
            neesL_id(:,lm_seq_id(i),k+1,kk) = err'*inv(Pe_id(3+2*i-1:3+2*i,3+2*i-1:3+2*i))*err;
            rmsL_id(:,lm_seq_id(i),k+1,kk) =  err'*err;
        end
        
        
        % Standard EKF
        xRest_std(:,k+1,kk) = xe_std(1:3,1);
        err = xR_true(1:3,k+1,kk) - xe_std(1:3,1);
        err(3) = pi_to_pi(err(3));
        xRerr_std(:,k+1,kk) = err;
        Prr_std(:,k+1,kk) = diag(Pe_std(1:3,1:3));
        neesR_std(:,k+1,kk) = err'*inv(Pe_std(1:3,1:3))*err;
        rmsRp_std(:,k+1,kk) = err(1:2,1)'*err(1:2,1);
        rmsRth_std(:,k+1,kk) = err(3,1)'*err(3,1);
        
        for i=1:length(lm_seq_std)
            xLest_std(:,lm_seq_std(i),k+1,kk) = xe_std(3+2*i-1:3+2*i,1);
            err = xL_true(:,lm_seq_std(i),kk) - xLest_std(:,lm_seq_std(i),k+1,kk);
            xLerr_std(:,lm_seq_std(i),k+1,kk) = err;
            Pll_std(:,lm_seq_std(i),k+1,kk) = diag(Pe_std(3+2*i-1:3+2*i,3+2*i-1:3+2*i));
            
            neesL_std(:,lm_seq_std(i),k+1,kk) = err'*inv(Pe_std(3+2*i-1:3+2*i,3+2*i-1:3+2*i))*err;
            rmsL_std(:,lm_seq_std(i),k+1,kk) =  err'*err;
        end
        
        
        % FEJ-EKF
        xRest_fej(:,k+1,kk) = xe_fej(1:3,1);
        err = xR_true(1:3,k+1,kk) - xe_fej(1:3,1);
        err(3) = pi_to_pi(err(3));
        xRerr_fej(:,k+1,kk) = err;
        Prr_fej(:,k+1,kk) = diag(Pe_fej(1:3,1:3));
        neesR_fej(:,k+1,kk) = err'*inv(Pe_fej(1:3,1:3))*err;
        rmsRp_fej(:,k+1,kk) = err(1:2,1)'*err(1:2,1);
        rmsRth_fej(:,k+1,kk) = err(3,1)'*err(3,1);
        
        for i=1:length(lm_seq_fej)
            xLest_fej(:,lm_seq_fej(i),k+1,kk) = xe_fej(3+2*i-1:3+2*i,1);
            err = xL_true(:,lm_seq_fej(i),kk) - xLest_fej(:,lm_seq_fej(i),k+1,kk);
            xLerr_fej(:,lm_seq_fej(i),k+1,kk) = err;
            Pll_fej(:,lm_seq_fej(i),k+1,kk) = diag(Pe_fej(3+2*i-1:3+2*i,3+2*i-1:3+2*i));
            
            neesL_fej(:,lm_seq_fej(i),k+1,kk) = err'*inv(Pe_fej(3+2*i-1:3+2*i,3+2*i-1:3+2*i))*err;
            rmsL_fej(:,lm_seq_fej(i),k+1,kk) =  err'*err;
        end
        
        % OC-EKF
        xRest_ocekf_1(:,k+1,kk) = xe_ocekf_1(1:3,1);
        err = xR_true(1:3,k+1,kk) - xe_ocekf_1(1:3,1);
        err(3) = pi_to_pi(err(3));
        xRerr_ocekf_1(:,k+1,kk) = err;
        Prr_ocekf_1(:,k+1,kk) = diag(Pe_ocekf_1(1:3,1:3));
        neesR_ocekf_1(:,k+1,kk) = err'*inv(Pe_ocekf_1(1:3,1:3))*err;
        rmsRp_ocekf_1(:,k+1,kk) = err(1:2,1)'*err(1:2,1);
        rmsRth_ocekf_1(:,k+1,kk) = err(3,1)'*err(3,1);
        
        for i=1:length(lm_seq_ocekf_1)
            xLest_ocekf_1(:,lm_seq_ocekf_1(i),k+1,kk) = xe_ocekf_1(3+2*i-1:3+2*i,1);
            err = xL_true(:,lm_seq_ocekf_1(i),kk) - xLest_ocekf_1(:,lm_seq_ocekf_1(i),k+1,kk);
            xLerr_ocekf_1(:,lm_seq_ocekf_1(i),k+1,kk) = err;
            Pll_ocekf_1(:,lm_seq_ocekf_1(i),k+1,kk) = diag(Pe_ocekf_1(3+2*i-1:3+2*i,3+2*i-1:3+2*i));
            
            neesL_ocekf_1(:,lm_seq_ocekf_1(i),k+1,kk) = err'*inv(Pe_ocekf_1(3+2*i-1:3+2*i,3+2*i-1:3+2*i))*err;
            rmsL_ocekf_1(:,lm_seq_ocekf_1(i),k+1,kk) =  err'*err;
        end
        
    end%end of all nSteps
    
    
end%end of monte carlo runs


%% Monte Carlo Results
% % average nees and rms of robot pose over all runs
neesR_avg_std = sum(neesR_std,3)/nRuns;
neesR_avg_id = sum(neesR_id,3)/nRuns;
neesR_avg_fej = sum(neesR_fej,3)/nRuns;
neesR_avg_ocekf_1 = sum(neesR_ocekf_1,3)/nRuns;

rmsRp_avg_std = sqrt(sum(rmsRp_std,3)/nRuns);
rmsRp_avg_id = sqrt(sum(rmsRp_id,3)/nRuns);
rmsRp_avg_fej = sqrt(sum(rmsRp_fej,3)/nRuns);
rmsRp_avg_ocekf_1 = sqrt(sum(rmsRp_ocekf_1,3)/nRuns);

rmsRth_avg_std = sqrt(sum(rmsRth_std,3)/nRuns);
rmsRth_avg_id = sqrt(sum(rmsRth_id,3)/nRuns);
rmsRth_avg_fej = sqrt(sum(rmsRth_fej,3)/nRuns);
rmsRth_avg_ocekf_1 = sqrt(sum(rmsRth_ocekf_1,3)/nRuns);

% % average nees and rms over landmarks and over time and runs
neesL_avg_std = sum(mean(mean(neesL_std(:,lm_seq_std,:,:),2),3),4)/nRuns;
neesL_avg_id = sum(mean(mean(neesL_id(:,lm_seq_id,:,:),2),3),4)/nRuns;
neesL_avg_fej = sum(mean(mean(neesL_fej(:,lm_seq_fej,:,:),2),3),4)/nRuns;
neesL_avg_ocekf_1 = sum(mean(mean(neesL_ocekf_1(:,lm_seq_ocekf_1,:,:),2),3),4)/nRuns;

rmsL_avg_std = sqrt(sum(mean(mean(rmsL_std(:,lm_seq_std,:,:),2),3),4)/nRuns);
rmsL_avg_id = sqrt(sum(mean(mean(rmsL_id(:,lm_seq_id,:,:),2),3),4)/nRuns);
rmsL_avg_fej = sqrt(sum(mean(mean(rmsL_fej(:,lm_seq_fej,:,:),2),3),4)/nRuns);
rmsL_avg_ocekf_1 = sqrt(sum(mean(mean(rmsL_ocekf_1(:,lm_seq_ocekf_1,:,:),2),3),4)/nRuns);

% % average robot pose err w/ cov
xRerr_avg_std = sum(xRerr_std,3)/nRuns;
xRerr_avg_id = sum(xRerr_id,3)/nRuns;
xRerr_avg_fej = sum(xRerr_fej,3)/nRuns;
xRerr_avg_ocekf_1 = sum(xRerr_ocekf_1,3)/nRuns;

% % avg. robot cov
Prr_avg_std = sum(Prr_std,3)/nRuns;
Prr_avg_id = sum(Prr_id,3)/nRuns;
Prr_avg_fej = sum(Prr_fej,3)/nRuns;
Prr_avg_ocekf_1 = sum(Prr_ocekf_1,3)/nRuns;

% % avg. landmark cov
%one landmark
Pll_avg_id = sum(squeeze(Pll_id),3)/nRuns;
Pll_avg_std = sum(squeeze(Pll_std),3)/nRuns;
Pll_avg_fej = sum(squeeze(Pll_fej),3)/nRuns;
Pll_avg_ocekf_1 = sum(squeeze(Pll_ocekf_1),3)/nRuns;


%% plot figures
plot_all

