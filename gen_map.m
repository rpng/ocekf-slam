function xL = gen_map(nL,v,omega,rmin,rmax, nSteps,dt)


if omega~=0

    % radius
    r = v/omega;

    for i=1:nL

        if rand>.5
            rho = r + rmin + (rmax-rmin)*rand;
            th = 2*pi*rand;
        else
            rho = r - rmax + (rmax-rmin)*rand;
            th = 2*pi*rand;
        end

        [x,y] = pol2cart(th,rho);
        xL(:,i) = [x;y+r]; %shift y w/ r since robot starts at (0,0)

    end
    
end



if omega==0 %exploration case
       
    for i=1:nL
        if round(rand)
            rho = rmin + (rmax-rmin)*rand;
            th = 2*pi*rand;
        else
            rho = - rmax + (rmax-rmin)*rand;
            th = 2*pi*rand;
        end
        [x,y] = pol2cart(th,rho);
        yy(1,i) = rho;%y;
    end
    
    dtraj = v*dt*nSteps;
    
    xL = [dtraj*rand(1,nL); yy ];
    
    xL = [1:dtraj/nL:dtraj; yy];

    xL = [dtraj*rand(1,nL); dtraj/2*(rand(1,nL)-0.5) ];

end

