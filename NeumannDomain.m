classdef NeumannDomain
    properties        
        maximum
        minimum
        saddles = cell(2,1)
        lines = cell(4,1)        
        resolution
        area
        circumference
        rho
        problems = [] %in case that there are problems in the gradient descent
        hess_eig = []
        COM %center of mass. Used for ploting the domain's diamond.
    end
    methods
        function obj = NeumannDomain(res, sad, low_stpt, high_stpt)           
            obj.saddles{1} = sad;           
            obj.resolution = res;
            if (length(low_stpt(:))>2)||(length(high_stpt(:))>2)
                fprintf('error');
            end
            obj.lines{1} = [sad;low_stpt];
            obj.lines{2} = [sad;high_stpt];
        end     
        function obj = get_line(obj, grad, hessian) %uses gradient descent to create the neumann lines         
            dl_ = 0.1*pi/obj.resolution;
            C = 0;
            S = 0;
            for n=1:2
                p = obj.lines{n}(2,:);
                G = p-obj.lines{n}(1,:);
                G = G/norm(G);
                while(1)   
                    %step                    
                    if sin(p(2))
                        dl = dl_*[1/abs(sin(p(2))), 1].*G;
                    else
                        dl = dl_*G;
                    end
%                     dl = dl_*[1/max(abs(sin(p(2))),0.001), 1].*G;
                    p = p + dl;                                       
                    obj.lines{n} = [obj.lines{n}; p];
                    C = C + norm(dl); %get length
                    S = S + ((-1)^n)*cos(p(2))*(dl(1)); %get area using green
                     p(1) = mod(p(1),2*pi);
%                     if abs(S)>3
%                         dl
%                     end                   
                    
                    %this is just checking my hypotesis that the hessian
                    %eigenvectores are tangent to the neumann lines. But
                    %this hypotesis is false.
                    if hessian~=0&&~mod(length(obj.lines{n}),100)
                        [HV,~] = eig(hessian(p));
                        HV = 0.03*[abs(sin(p(2))),0;0,1]*HV';
                        HV = p + [HV; -HV];
                        obj.hess_eig = cat(3, obj.hess_eig, HV);
                    end
                    
                    G_ =  ((-1)^n)*grad(p);
                    G_ = (2/pi)*G_/norm(G_);
                    if dot(G,G_)<0
                        if length(obj.lines{n})<10                            
                            fprintf('problem: gradient return at early stage. point: (%f,%f)\n',p(1)/(2*pi),p(2)/pi);
                            obj.problems = [obj.problems; p];
                            G = G_;                       
                        else
                            break;
                        end
                    else
                        G = G_;
                    end                   
                end
            end
            obj.circumference = C;
            obj.area = S;
        end
        function obj = get_associated_extrima(obj, maxima, minima) %finds the neumann domain's extrema from the lists of extrema           
            [~,i] = min(pdist2(obj.spheric_mod(obj.lines{1}(end,:)),minima));
            obj.minimum = minima(i,:);   
            [~,i] = min(pdist2(obj.spheric_mod(obj.lines{2}(end,:)),maxima));
            obj.maximum = maxima(i,:);
            if isempty(obj.minimum)||isempty(obj.maximum)
                fprintf('min or max was not found. point: (%f,%f)',obj.saddles{1}(1),obj.saddles{1}(2));
            end
        end
        function obj = devour(obj, obj2, L) %the story when one neumann domain devours its evil twin
            obj.saddles{2} = obj2.saddles{1};
            obj.lines{3} = obj2.lines{1};
            obj.lines{4} = obj2.lines{2};
            
            obj.circumference = obj.circumference + obj2.circumference;
            obj.area = abs(obj.area-obj2.area);
            obj.rho = sqrt(L)*obj.area/obj.circumference;        
            
            [x1, y1, z1] = sph2car(obj.saddles{1}(:,1), obj.saddles{1}(:,2), 1);
            [x2, y2, z2] = sph2car(obj.saddles{2}(:,1), obj.saddles{2}(:,2), 1);
            m = mean([x1, y1, z1; x2, y2, z2], 1);
            [ph, th] = cart2sph(m(1), m(2), m(3));
            obj.COM =  [mod(ph, 2*pi), pi/2-th];
            
            if obj.area>5
                obj.area
            end
            
        end          
        function graph_ = print_dom(obj, color, line_wdth) %ploting method for the neuman lines
            graph = [];
            for line = obj.lines(:)'             
                [x, y, z] = sph2car(line{:}(:,1), line{:}(:,2), 1);
                graph = [graph; plot3(x, y, z,'color', color, 'linewidth', line_wdth)];
            end
            graph_ = graph;
        end
    end
    methods (Static)
        function p_ = spheric_mod(p)
            if p(2)<0
                p(2) = -p(2);
                p(1) = p(1) + pi;
            elseif p(2)>pi
                p(2) = p(2) - pi;
                p(1) = p(1) + pi;
            end
            p(1) = mod(p(1),2*pi);
            p_ = p;
        end
    end
end