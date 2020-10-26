function coo_ = get_coordinates(res,wb)
    waitbar(0,wb,'Getting coordinate system...');
    coo = res;
    n = coo.resolution;
    coo.epsilon = 2*pi/n;   
    
    
    coo.ax.theta = linspace(0,pi,n);
    coo.ax.phi = linspace(0,2*pi,2*n);    
    [coo.phi,coo.theta] = meshgrid(coo.ax.phi,coo.ax.theta);
    
    coo.d.theta = n/pi;
    coo.d.phi = n*diag(1./sin(coo.ax.theta));    
    
    coo_ = coo;