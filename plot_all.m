close all

% % setup 
for nr= 1%1:nRuns
    figure
    plot(xR_true(1,:,nr), xR_true(2,:,nr), 'b-.','Linewidth',1), hold on
    plot(xL_true(1,:,nr), xL_true(2,:,nr), 'b*','Linewidth',1)
    title('Robot trajectory and landmarks')
    legend('True Trajectory','True Landmarks')
    xlabel('x (m)','FontWeight','bold'), ylabel('y (m)','FontWeight','bold')
    axis equal
    print('-dpng','setup_sim')
end



% % robot % %

start= 1;  incr= 5;

% Robot NEES
figure, hold on
plot([start:incr:nSteps],neesR_avg_id(start:incr:end)','g:','Linewidth',3);
plot([start:incr:nSteps],neesR_avg_std(start:incr:end)','b-o','Linewidth',1);
plot([start:incr:nSteps],neesR_avg_fej(start:incr:end)','m-.','Linewidth',2);
plot([start:incr:nSteps],neesR_avg_ocekf_1(start:incr:end)','k--','Linewidth',2);

xlabel('Time (sec)','FontWeight','bold'), ylabel('Robot pose NEES','FontWeight','bold')
legend( 'Ideal-EKF','Std-EKF','FEJ-EKF','OC-EKF')%,'OC-EKF3','Robocentric')
print('-dpng','slam_robot_nees_oc3_sim')

% Robot RMSE
figure
subplot(2,1,1), hold on
plot([start:incr:nSteps],rmsRp_avg_id(start:incr:end)','g:','Linewidth',3);
plot([start:incr:nSteps],rmsRp_avg_std(start:incr:end)','b-o','Linewidth',1);
plot([start:incr:nSteps],rmsRp_avg_fej(start:incr:end)','m-.','Linewidth',2);
plot([start:incr:nSteps],rmsRp_avg_ocekf_1(start:incr:end)','k--','Linewidth',2);
ylabel('Position RMSE (m)','FontWeight','bold')
legend( 'Ideal-EKF','Std-EKF','FEJ-EKF','OC-EKF')%,'OC-EKF3','Robocentric') %#ok<*LEGINTPAR>
subplot(2,1,2), hold on
plot([start:incr:nSteps],rmsRth_avg_id(start:incr:end)','g:','Linewidth',3);
plot([start:incr:nSteps],rmsRth_avg_std(start:incr:end)','b-o','Linewidth',1);
plot([start:incr:nSteps],rmsRth_avg_fej(start:incr:end)','m-.','Linewidth',2);
plot([start:incr:nSteps],rmsRth_avg_ocekf_1(start:incr:end)','k--','Linewidth',2);
xlabel('Time (sec)','FontWeight','bold'),
ylabel('Heading RMSE (rad)','FontWeight','bold')
print('-dpng','slam_robot_rms_oc3_sim')


NEES_Robot = [mean(neesR_avg_id(1:end)),mean(neesR_avg_std(1:end)), mean(neesR_avg_fej(1:end)), mean(neesR_avg_ocekf_1(1:end)) ]
RMS_Position = [mean(rmsRp_avg_id(1:end)),mean(rmsRp_avg_std(1:end)), mean(rmsRp_avg_fej(1:end)),   mean(rmsRp_avg_ocekf_1(1:end)) ]
RMS_Heading = [mean(rmsRth_avg_id(1:end)),mean(rmsRth_avg_std(1:end)), mean(rmsRth_avg_fej(1:end)), mean(rmsRth_avg_ocekf_1(1:end))]



% % landmarks % %

% NEES
figure
bar([neesL_avg_id,...
    neesL_avg_std,...
    neesL_avg_fej,...
    neesL_avg_ocekf_1 ]);
ylabel('Avg. Landmark NEES','FontWeight','bold')
print('-dpng','lm_nees')


% RMSE
figure
bar([rmsL_avg_id,...
    rmsL_avg_std,...
    rmsL_avg_fej,...
    rmsL_avg_ocekf_1  ]);
ylabel('Avg. Landmark RMSE','FontWeight','bold')
print('-dpng','lm_rms')


NEES_L = [neesL_avg_id,neesL_avg_std,neesL_avg_fej, neesL_avg_ocekf_1]
RMS_L = [rmsL_avg_id,rmsL_avg_std,rmsL_avg_fej, rmsL_avg_ocekf_1 ]


