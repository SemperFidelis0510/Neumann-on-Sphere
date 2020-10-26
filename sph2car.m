function [x,y,z] = sph2car(phi,theta,r)
    [x,y,z] = sph2cart(phi,(pi/2)-theta,r);
end