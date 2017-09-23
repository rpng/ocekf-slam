function [xe,Pe,xL_1,lm_seq,V_fej] = update_fej(xe,Pe,xL_1,lm_seq,z,R,V_fej)


lenz = length(find(z(3,:)>0));
lenx= size(xe,1);

global gDISTBEAR

xR_k_k1 = xe(1:3,1);

nf = 0;

% % update
for i= 1:lenz
    % data association (based on landmark id)
    is_exist = ~(lm_seq - z(3,i));
    idx = find(is_exist);
    
    ii = 2*i+(-1:0);
    
    % update: already in the state vecor
    if ~isempty(idx)
        nf = nf+1;
        jj = 2*nf+(-1:0);
        if gDISTBEAR
            [zhat,Hii] = measurement_model_fej(xe,idx,xL_1);   %Hfej=Hii
            H(jj,:) = Hii;
            r(jj,1) = [ z(1,i)-zhat(1,1);  pi_to_pi(z(2,i)-zhat(2,1)) ];
        else
            [zhat,Hii] = measurement_model_fej(xe,idx,xL_1);   %Hfej=Hii
            H(jj,:) = Hii;
            r(jj,1) = [ z(1,i)-zhat(1,1);  (z(2,i)-zhat(2,1)) ];
        end
        Rf(jj,jj) = R(ii,ii);
    end
end

if nf~=0
    S = H*Pe*H'+ Rf;
    S = (S+S')*0.5;
    
    
    if isspd(S) 
        
        K = Pe*H'/S;
        xe = xe + K*r;
        Pe = (eye(length(Pe)) - K*H) * Pe *(eye(length(Pe)) - K*H)' + K*Rf*K';
        
    end
end



% % augment
for i= 1:lenz
    % data association (known)
    is_exist = ~(lm_seq - z(3,i));
    idx = find(is_exist);
    
    lenx= size(xe,1);
    ii = 2*i + (-1:0);
    
    % add the new landmark into the state vector
    if isempty(idx)
        lm_seq = [lm_seq; z(3,i)];
        
        
        if gDISTBEAR
            % augment state
            d = z(1,i);
            th = z(2,i);
            k_xL = [d*cos(th); d*sin(th)];
            
            x_L = xe(1:2,1) + [ d*cos(th+xe(3)); d*sin(th+xe(3)) ];
            xe = [xe; x_L];
            
            xL_1 = [xL_1; x_L]; %
            
            
            % jacobians
            J = [0 -1; 1 0];
            C = [cos(xR_k_k1(3))  -sin(xR_k_k1(3));  sin(xR_k_k1(3))  cos(xR_k_k1(3)) ];
            
            H_Lk = [ 1/norm(k_xL)*k_xL'; 1/norm(k_xL)^2*k_xL'*J' ];
            
            HR = - H_Lk* C'*[eye(2)  J*(x_L-xe(1:2,1))];
            HL = H_Lk* C';
        else
            % augment state
            k_xL = z(1:2,i);
            
            C = [cos(xe(3)) -sin(xe(3));  sin(xe(3)) cos(xe(3)) ];
            x_L = xe(1:2,1) + C*k_xL;
            
            xe = [xe; x_L*1];
            xL_1 = [xL_1; x_L ]; %
            
            
            % jacobians
            J = [0 -1; 1 0];
            C = [cos(xR_k_k1(3))  -sin(xR_k_k1(3));  sin(xR_k_k1(3))  cos(xR_k_k1(3)) ];
            
            HR = - C'*[eye(2)  J*(x_L-xe(1:2,1))];
            HL = C';
        end
        
        V_fej = null([HR,HL]);
        
        
        % augment covariance
        rng= lenx+1:lenx+2;
        Pe(rng,rng)= inv(HL)*HR*Pe(1:3,1:3)*HR'*inv(HL)' + inv(HL)*R(ii,ii)*inv(HL)'; % landmark cov
        Pe(rng,1:3)= -inv(HL)*HR*Pe(1:3,1:3); % landmark-robot xcorr
        Pe(1:3,rng)= Pe(rng,1:3)';
        if lenx>3
            rnm= 4:lenx;
            Pe(rng,rnm)= -inv(HL)*HR*Pe(1:3,rnm);
            Pe(rnm,rng)= Pe(rng,rnm)';
        end
        Pe = 0.5*(Pe+Pe');
        
    end
end

