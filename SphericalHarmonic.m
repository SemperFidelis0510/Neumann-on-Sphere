classdef SphericalHarmonic
    properties %function identity
        coeffs %trial function coefficients
        nodal_doms = 0; %# of nodal domains    
        maxDeg = [] %maxima degrees
        minDeg = [] %minima degrees      
    end
    properties %function elements
        func = 0 %the handle for the trial function itself, as an anonymous function        
        grad = cell(2,1) %the trial function's gradient      
        hessian = cell(2,2) %the hessian of the trial function
        values = [] %function values on the grid points
        nodal = [] %coordinates of the nodal points
        maxima = [] %coordinates of the maxima       
        minima = [] %coordinates of the minima
        saddles = [] %coordinates of the saddles
        neumann = cell(0) %array of the nodal domains
    end
    methods %calculation methods
        function obj = SphericalHarmonic(coeffs) %just creates the object 
            %input: coefficients structure. Must be as described in
            %"neumann_on_sphere.m"
            obj.coeffs = coeffs;            
        end
        function obj = get_function(obj, wb) %gets function handle, gradient and hessian
            t = sym('t', 'real'); %theta symbolic variable
            p = sym('p', 'real'); %phi symbolic variable
            C = obj.coeffs.vals;
            l = obj.coeffs.l;
            phi_part = sym(zeros(1,l+1));
            dphi_part = sym(zeros(1,l+1));

            %creates the spherical harmonics and their gradients and then
            %the hessian simultaniously
            for m=0:l
                coef_m = ((-1)^m)*(2^(-l))*sqrt(factorial(l-m)/factorial(l+m));
                if ~m
                    phi_part(1) = coef_m*C(l+1);
                    dphi_part(1) = 0;
                else
                    phi_part(m+1) = simplify(coef_m*(C(l+m+1)*cos(m*p)+C(l-m+1)*sin(m*p)));
                    dphi_part(m+1) = simplify(m*coef_m*(-C(l+m+1)*sin(m*p)+C(l-m+1)*cos(m*p))/sin(t));
                end
                if wb~=0
                    waitbar(0.05+(0.05*m/l),wb,'Getting function coefficients...');
                end
            end

            legendre = simplify(diff((t^2-1)^l, l));    
            theta_part = simplify(subs(legendre(1), t, cos(t)));
            dtheta_part = simplify(diff(theta_part, t));

            for m=1:l
                legendre = [legendre, simplify(diff(legendre(end),t))]; %#ok<*AGROW>
                theta_part = [theta_part, simplify((sin(t)^m)*subs(legendre(end),t,cos(t)))];
                dtheta_part = [dtheta_part, simplify(diff(theta_part(end),t))];
                if wb~=0
                    waitbar(0.1+(0.3*m/l),wb,sprintf('Getting function m components... %d/%d',2*m+1,2*l+1));
                end
            end 
            if wb~=0
                waitbar(0.4,wb,'Getting function...');
            end
            obj.func = matlabFunction(dot(theta_part, phi_part), 'vars', {'p', 't'});    
            if wb~=0
                waitbar(0.45,wb,'Getting gradient...');
            end
            obj.grad{1} = matlabFunction(dot(theta_part, dphi_part), 'vars', {'p','t'});
            obj.grad{2} = matlabFunction(dot(dtheta_part, phi_part), 'vars', {'p','t'});   
            
            if obj.coeffs.hess_chk
                obj.hessian{1,1} = matlabFunction(diff(dot(theta_part,dphi_part),p),'vars',{'p','t'});
                obj.hessian{1,2} = matlabFunction(diff(dot(theta_part,dphi_part),t)+(cos(t)/sin(t))*dot(theta_part,dphi_part),'vars',{'p','t'});
                obj.hessian{2,2} = matlabFunction(diff(dot(dtheta_part,phi_part)+(sin(2*t)/2)*dot(dtheta_part,phi_part),t),'vars',{'p','t'});
                obj.hessian{2,1} = obj.hessian{1,2};
            end
        end
        function obj = get_values(obj, grid) %gets function values on the grid
            obj.values = obj.func(grid{1},grid{2});
        end
        function obj = get_nodal(obj, grid) 
                N = find_zeros(obj.values); %get nodal lines
                obj.nodal = [grid{1}(N), grid{2}(N)];
                obj.nodal_doms = HoshenKop(double(N)); %counts nodal domains
        end
        function obj = get_extrema(obj, coo)    
            
            res = coo.resolution;
            angRes = coo.angularResolution;
            V = obj.func(coo.grid{1},coo.grid{2});
            dydt = diff(V,1,1); %theta variation of function
            dydt = [dydt; dydt(end,:)];
            dydp = diff(V,1,2); %phi variation of function
            dydp = [dydp, dydp(:,1)];

            %gets the intersection of the nodal lines of the gradient components
            susp_points = find_zeros(dydt)&find_zeros(dydp);    
            if all(sign(dydt(1,:))==sign(dydt(1,1)))
                susp_points(1,:) = ones(1,2*res);
            else
                susp_points(1,:) = zeros(1,2*res);
            end
            if all(sign(dydt(end,:))==sign(dydt(end,1)))
                susp_points(end-1,:) = ones(1,2*res);
                susp_points(end,:) = ones(1,2*res);
            else
                susp_points(end-1,:) = zeros(1,2*res);
                susp_points(end,:) = zeros(1,2*res);
            end         
            %seperates the clusters of zeros located around each
            %extrema/saddle into a list of coordinates for the suspected
            %points called "susp_points"
            [~,colors, susp_points]=HoshenKop(double(susp_points)); 
            %refines the locations of suspected points.
            %this seems to me as the most important part of the code in
            %preventing and diminishing bugs. Here one should implement
            %Newton-Raphson to accurate the location of the point. A trial
            %for such implementation is located in the working directory.
            susp_points = obj.refine_susp(coo, colors, susp_points);            
            
            if isempty(susp_points)
                return;
            end
            
            %sort suspected points into minima/maxima/saddles
            surr2 = 5*(pi/res)*[cos(linspace(0,2*pi,angRes))',sin(linspace(0,2*pi,angRes))'];
            
            for i=1:length(susp_points(:,1))
                p = susp_points(i,:)+(pi/(2*res))*ones(1,2);
                if sin(p(2))
                    surr = ([1/sin(p(2))^2, 1; 0, 1]*surr2')';
                else
                    surr = surr2;
                end
                if ~isempty([obj.minima;obj.maxima;obj.saddles])                 
                    if any(pdist2(p,[obj.minima;obj.maxima;obj.saddles])<=(50*abs(sin(p(2)))*pi/(res*sqrt(obj.coeffs.l))))
                        continue;
                    end
                end   
                convexity = obj.func(p(1)+surr(:,1),p(2)+surr(:,2)) + obj.func(p(1)-surr(:,1),p(2)-surr(:,2)) - 2*obj.func(p(1),p(2));
                if all(convexity>=0) %minimum found
                    obj.minima = [obj.minima; p];
                elseif all(convexity<=0) %maximum found
                    obj.maxima = [obj.maxima; p];
                else %saddle found
                    obj.saddles = [obj.saddles; p];
                    %now find the starting points for the 4 neumann lines
                    %lot of bugs come from here as the gradient descent
                    %might be unstable in some cases, like when close to
                    %the poles.
                    surr_samp = p + surr;
                    [~,M] = max(obj.func(surr_samp(:,1),surr_samp(:,2)));
                    [~,m] = min(obj.func(surr_samp(:,1),surr_samp(:,2)));
                    avg = (surr(m,:)+surr(M,:))/2;
                    if det([avg;surr(m,:)])>0
                        U =(1/sqrt(2))*[1,0;0,sin(p(2))]*[1,-1;1,1];
                        Uinv = (1/sqrt(2))*[1,0;0,sin(p(2))]*[1,1;-1,1];
                        R = [U*avg(:),Uinv*avg(:)]';
                    else
                        U = (1/sqrt(2))*[1,0;0,sin(p(2))]*[1,-1;1,1];
                        Uinv = (1/sqrt(2))*[1,0;0,sin(p(2))]*[1,1;-1,1];                       
                        R = [Uinv*avg(:),U*avg(:)]';
                    end
                    for km=1:2 %creates 4 neumann lines object with now only the starting point for the lines
                        for kM=1:2                            
                            obj.neumann = [obj.neumann;NeumannDomain(res,p,p+((-1)^km)*R(1,:),p+((-1)^kM)*R(2,:))];
                        end
                    end
                end
            end
        end
        function obj = get_neumann_lines(obj, wb)                  
            grad_ = @(v) [obj.grad{1}(v(1),v(2)),obj.grad{2}(v(1),v(2))];
            if obj.coeffs.hess_chk
                hessian_ = @(v) [obj.hessian{1,1}(v(1),v(2)),obj.hessian{1,2}(v(1),v(2));obj.hessian{2,1}(v(1),v(2)),obj.hessian{2,2}(v(1),v(2))];
            else
                hessian_ = 0;
            end
            l = length(obj.neumann);
            for i=1:l
                if wb~=0
                    waitbar(0.6+i*0.35/l,wb,sprintf('Getting Neumann lines... %d/%d',ceil(i/2),l/2));
                end
                obj.neumann(i) = get_line(obj.neumann(i), grad_, hessian_); %creates the neuman line with gradient descent
            end
        end
        function obj = merge_neumann_lines(obj)
            %since now we have that each domain had only 2 lines, there is
            %a multiplicity of neumann domains, and we have to merge each
            %domain with the counter domain that is associated with the
            %same extrema. Hence we first assign each domain its extrema,
            %and then merge any two domains with the same extrema.
            for i=1:length(obj.neumann)
                obj.neumann(i) = obj.neumann(i).get_associated_extrima(obj.maxima,obj.minima);
            end
            i = 1;
            new_neumann = [];
            while(i<length(obj.neumann))
                for j=(i+1):length(obj.neumann)
                    if isequal(obj.neumann(i).minimum,obj.neumann(j).minimum)&&isequal(obj.neumann(i).maximum,obj.neumann(j).maximum)
                        new_neumann = [new_neumann; obj.neumann(i).devour(obj.neumann(j),obj.coeffs.l)];
                        break;
                    end
                end
                i = i + 1;
            end
            obj.neumann = new_neumann;
        end
        function obj = get_extrema_deg(obj) %count the number of domains around each extrema
                maxDeg = zeros(length(obj.maxima(:,1)),1);
                minDeg = zeros(length(obj.minima(:,1)),1);
                for k=1:length(obj.neumann)
                        for i=1:length(maxDeg)
                                if isequal(obj.neumann(k).maximum, obj.maxima(i,:))
                                        maxDeg(i) = maxDeg(i)+1;
                                        break;
                                end
                        end
                        for j=1:length(minDeg)
                                if isequal(obj.neumann(k).minimum, obj.minima(j,:))
                                        minDeg(j) = minDeg(j)+1;
                                        break;
                                end
                        end
                end
                obj.maxDeg = maxDeg; %#ok<*PROP>
                obj.minDeg = minDeg;
        end
        function X = refine_susp(obj, coo, colors, M) %takes the center of mass of each cluster
            X_ = [];
            for i=1:length(colors)           
                cluster = [coo.grid{1}(M==colors(i)), coo.grid{2}(M==colors(i))];
                if size(cluster,1)<coo.resolution
                    X_ = [X_; mean(cluster,1)];
                elseif cluster(1,2)<(2*pi/coo.resolution)||cluster(1,2)>(pi*(1-2/coo.resolution))
                    X_ = [X_; mean(cluster,1)];
                end
            end
            X = X_;
        end
    end        
    methods %plot methods
        function graph = print_function(obj, grid)
            [x, y, z] = sph2car(grid{1}, grid{2}, 1);
            graph = surf(x, y, z, obj.values, 'EdgeColor', 'none');
        end
        function graph = print_nodal(obj, pt_size, pt_style)
            [x, y, z] = sph2car(obj.nodal(:,1), obj.nodal(:,2), 1);
            graph = scatter3(x, y, z, pt_size, pt_style);
        end
        function [graph_max, graph_min] = print_extrema(obj, rad, pt_size, max_style, min_style)
            if ~isempty(obj.maxima)
                [x, y, z] = sph2car(obj.maxima(:,1), obj.maxima(:,2), rad);
                graph_max = scatter3(x, y, z, pt_size, max_style, 'filled');                
            end            
           
            if ~isempty(obj.minima)
                [x, y, z] = sph2car(obj.minima(:,1), obj.minima(:,2), rad);
                graph_min = scatter3(x, y, z, pt_size, min_style, 'filled');                
            end  
        end
        function graph = print_saddles(obj, rad, pt_size, color)
            if isempty(obj.saddles)
                graph = [];
                return;
            end
                [x, y, z] = sph2car(obj.saddles(:,1), obj.saddles(:,2), rad);
                graph = scatter3(x, y, z, pt_size, color, 'filled');                        
        end
        function [graphs_, graph_prob] = print_neumann(obj, color, line_wdth)
            if isempty(obj.neumann)
                graphs_ = [];
                graph_prob = [];
                return;
            end         
            graphs = [];      
            problems = [];
           
           for neu = obj.neumann(:)'
                   problems = [problems; neu.problems];                 
                   
                   graphs = [graphs; neu.print_dom(color, line_wdth)];                    
           end

           if ~isempty(problems)
               [x, y, z] = sph2car(problems(:,1), problems(:,2), 1);
               graph_prob = scatter3(x, y, z, 'ro');
           else
               graph_prob = [];
           end
           
           graphs_ = graphs;           
        end
        function [graph, marks, marked] = print_diamonds(obj, color, diam_size, rad)
            if isempty(obj.neumann)
                graph = [];    
                marked = [];
                marks = [];
                return;
            end  
            domains = [];
            rhoes = [];
            for neu = obj.neumann(:)'
                rhoes = [rhoes; neu.rho];
                [x, y, z] = sph2car(neu.COM(:,1), neu.COM(:,2), 1);
                domains = [domains; x, y, z];
            end
            rho_colors = (rhoes-min(rhoes))*color/max(rhoes);
           graph = scatter3(domains(:,1), domains(:,2), domains(:,3), diam_size, rho_colors, 'd', 'filled');
           marks = rad*domains;
           marked = [];
        end
        function graphs_ = print_hess(obj, color, line_wdth)    
            if isempty(obj.neumann)
                graphs_ = [];                   
                return;
            end  
            graphs = [];
            for neu = obj.neumann(:)'
                for j=1:size(neu.hess_eig,3)
                    [x, y, z] = sph2car(neu.hess_eig(:,1,j), neu.hess_eig(:,2,j), 1);
                    for s=1:2
                        graphs = [graphs; plot3(x([s,s+2]), y([s,s+2]), z([s,s+2]), 'color', color, 'linewidth', line_wdth)];
                    end
                end
            end
            graphs_ = graphs;
        end
        function [graph, ax_lim_] = print_rhoes(obj, ax, vis)
            if isempty(obj.neumann)
                graph = [];                   
                return;
            end  
            rhoes = [];
            for neu = obj.neumann(:)'
                rhoes = [rhoes, neu.rho];
            end
            [y, x] = histcounts(rhoes);
            graph = plot(ax, x(1:end-1), y, 'Visible', vis);
            ax_lim =[x(1), x(end), 0, 1.5*max(y)];
            if strcmp(vis,'on')
                axis(ax, ax_lim);
            end
            ax_lim_ = ax_lim;
        end
        function [graph, ax_lim_] = print_extdeg(obj, ax, vis)
             M = max([obj.maxDeg; obj.minDeg]);
             if M>0
                 graph = [histogram(ax, obj.maxDeg, (0:1:M), 'FaceColor', 'r', 'Visible',vis)...
                     ;histogram(ax, obj.minDeg, (0:1:M), 'FaceColor', 'b', 'Visible', vis)];
                  ax_lim =[0, M, 0, max(length(obj.maxima), length(obj.minima))] ;
                  if strcmp(vis,'on')
                      axis(ax, ax_lim);
                  end
                  ax_lim_ = ax_lim;
             end
        end
    end
end