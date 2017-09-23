function [D] = cal_kld(xR1,xR2,Pr,th1,th2,Rz)
%
% compute kld or relative entropy D(f||g)
% see Bailey ICRA 2003
%


nSamples = 100; % modify if necessary


% step 1
mu = [xR1;xR2;th1;th2]'; %original state w/ bearing measuremnts
P =  blkdiag(Pr,Rz);

spd = isspd(P);
if ~spd
    D = -1;
    return
end

X = mvnrnd(mu,P,nSamples);
for k = 1: nSamples
    w(k) = mvnpdf(X(k,:),mu,P);
end


% step 2-3
for k = 1: nSamples
    x1 = X(k,1:3)';
    x2 = X(k,4:6)';

    psi1 = pi_to_pi(x1(3)+th1);
    s1 = sin(psi1);
    c1 = cos(psi1);
    psi2 = pi_to_pi(x2(3)+th2);
    s2 = sin(psi2);
    c2 = cos(psi2);

    xf = 1/(s1*c2-s2*c1) * ...
        [ x1(1)*s1*c2 - x2(1)*s2*c1 + (x2(2)-x1(2))*c1*c2;
        x2(2)*s1*c2 - x1(2)*s2*c1 + (x1(1)-x2(1))*s1*s2 ];

    Xt(k,:) = [x1;x2;xf]'; %transformed samples with landmark position

    detHf = ( x1(1)*(xf(2)-x2(2)) + x2(1)*(x1(2)-xf(2)) + xf(1)*(x2(2)-x1(2)) ) / ( norm(xf-x1(1:2))^2 * norm(xf-x2(1:2))^2 );
    wt(k) = w(k) * abs(detHf);
end


% step 4
x1 = xR1;
x2 = xR2;

psi1 = pi_to_pi(x1(3)+th1);
s1 = sin(psi1);
c1 = cos(psi1);
psi2 = pi_to_pi(x2(3)+th2);
s2 = sin(psi2);
c2 = cos(psi2);

Hx1 = 1/(s1*c2-s2*c1) * [ s1*c2; s1*s2 ];
Hy1 = 1/(s1*c2-s2*c1) * [ -c1*c2; -s2*c1 ];
Hphi1 = -(c1*c2+s2*s1)/(s1*c2-s2*c1)^2 * [ x1(1)*s1*c2-x2(1)*s2*c1+(x2(2)-x1(2))*c1*c2; x2(2)*s1*c2-x1(2)*s2*c1+(x1(1)-x2(1))*s1*s2 ] + ...
    1/(s1*c2-s2*c1) * [ x1(1)*c1*c2+x2(1)*s2*s1-(x2(2)-x1(2)*s1*c2); x2(2)*c1*c2+x1(2)*s2*s1+(x1(1)-x2(1))*c1*s2 ];

HR1 = [Hx1, Hy1, Hphi1];

Hx2 = 1/(s1*c2-s2*c1) * [ -s2*c1; -s1*s2 ];
Hy2 = 1/(s1*c2-s2*c1) * [ c1*c2; s1*c2 ];
Hphi2 = (s1*s2+c1*c2)/(s1*c2-s2*c1)^2 * [ x1(1)*s1*c2-x2(1)*s2*c1+(x2(2)-x1(2)*c1*c2); x2(2)*s1*c2-x1(1)*s2*c1+(x1(1)-x2(1))*s1*s2; ] + ...
    1/(s1*c2-s2*c1) * [ -x1(1)*s1*s2-x2(1)*c2*c1-(x2(2)-x1(2))*c1*s2; -x2(2)*s1*s2-x1(2)*c2*c1+(x1(1)-x2(1))*s1*c2 ];

HR2 = [Hx2, Hy2, Hphi2];

Hr = [HR1, HR2];

Hz1 = Hphi1;
Hz2 = Hphi2;

Hz = [Hz1, Hz2];

F = eye(8);
F(end-1:end,:) = [Hr Hz];

Pt = F*P*F'; %covariance of transformed samples

spd = isspd(Pt);
if ~spd
    D = -1;
    return
end


xf = 1/(s1*c2-s2*c1) * ...
    [ x1(1)*s1*c2 - x2(1)*s2*c1 + (x2(2)-x1(2))*c1*c2;
    x2(2)*s1*c2 - x1(2)*s2*c1 + (x1(1)-x2(1))*s1*s2 ];

mut = [x1;x2;xf]'; % mean of the tranformed samples
% mut = mean(Xt,1);



% step 5
for k = 1: nSamples
    v(k) = mvnpdf(Xt(k,:),mut,Pt);
end


% step 6
D = sum(log(w)-log(v)) / nSamples;

