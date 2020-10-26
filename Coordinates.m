classdef Coordinates
    properties
        resolution %grid resolution
        angularResolution %resolution for circular sampling around a point
        colorResolution %resolution for plot colors
        pointSize %size of extrema points in plot
        uncertainty %mostly an unsuccessful tool to help prevent bugs
        grid = cell(2,1)
        g
        ginv
    end
    methods
        function obj = Coordinates(res,ang_res,clr_res,pt_sz)
            %just orgenize the inputs from the GUI
            obj.resolution = res; 
            obj.angularResolution = ang_res; 
            obj.colorResolution = clr_res; 
            obj.pointSize = pt_sz; 
            obj.uncertainty = pi/res; 
            
            %creating the grid itself
            [obj.grid{1},obj.grid{2}] = meshgrid(linspace(0,2*pi,2*res),linspace(0,pi,res));
            
            obj.g = @(theta) [sin(theta).^2,1]; %the riemannian metric
            obj.ginv = @(theta) [1./sin(theta).^2,1]; %its inverse
        end
    end
end