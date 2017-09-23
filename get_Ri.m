function Ri = get_Ri(r_i, theta_i, sigma_r, sigma_theta)
% computes the covariance of the point p_i = r_i*[cos(theta_i); sin(theta_i)
% in cartesian coordinates


% first we compute the matrix square root of R, which consists of the
% Jacobian * sqrt(R_in_polar_coords)
Ri = [cos(theta_i) -r_i*sin(theta_i)
      sin(theta_i)  r_i*cos(theta_i) ] * diag([sigma_r sigma_theta]); 

% and then we get the entire Ri   
Ri = Ri * Ri';