function graphs_ = plot_function(coo, SH, ax, rho_ax, wb)
    if wb~=0
        waitbar(0.95,wb,'All function data acquired. Plotting everything...');
    end
    axes(ax);
    colormap([linspace(0,0.7,coo.colorResolution)',0.1*ones(coo.colorResolution,1),linspace(0.7,0,coo.colorResolution)']);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    colorbar;
    daspect([1 1 1]); 
    ED = 1.01; 
    legs = {'Nodal lines'};
    S = coo.pointSize;
    
    %plot function
    graphs.SH = SH.print_function(coo.grid);    
    graphs.SH.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
 
    %plot nodal dots
    graphs.nodal = SH.print_nodal(10, 'k.');
    
    
    %plot maxima
    if ~isempty(SH.maxima)&&~isempty(SH.minima)
        [graphs.maxima, graphs.minima]= SH.print_extrema(ED, S, 'ro', 'bo');
        legs = [legs, 'Maxima', 'Minima'];
    end 
    
    
    %plot saddles
    if ~isempty(SH.saddles)        
        graphs.saddles = SH.print_saddles(ED, S, [0.9,0.9,0.9]);
        legs = [legs, 'Saddles'];
    end
    
    
    if ~isempty(SH.neumann)
        %plot neumann lines
        [graphs.neumann_lines, graphs.problems] = SH.print_neumann([0.5,0.5,0.5], 1.5);
        
        %plot diamonds
        [graphs.neu_domains_marks, graphs.neu_domains, graphs.neu_marked] = SH.print_diamonds([0,1,0], 50, ED);
        
        %hessian eigenvectors
        if SH.coeffs.hess_chk
            graphs.hessian = SH.print_hess([0,1,0], 1);
        end
        
        legs = [legs, 'Neumann lines'];    
        
        %rho distribution
        [graphs.rho, graphs.rho_M] = SH.print_rhoes(rho_ax, 'on');
        
        %extreme degrees
        [graphs.extDeg, graphs.extDeg_M] = SH.print_extdeg(rho_ax, 'off');
        
        axes(ax);
    end
    

    %get z and x axes
    plot3([0,0],[0,0],[-1.5,1.5],'color',[0,0.5,0.5],'linewidth',1);
    plot3([-1.5,1.5],[0,0],[0,0],'color',[0,0.5,0.5],'linewidth',1);
    
    
    %get legends
    legend(legs,'location','northeast');
    
    
    %some stats for table
    switch SH.coeffs.style
        case 'uni'
            tag = 'Uniform';
        case 'gauss'
            tag = 'Gaussian';
        case 'manual'
            tag = 'Manual';
    end
    graphs.data = {SH.coeffs.l, tag, SH.nodal_doms, length(SH.neumann), length(SH.maxima), length(SH.minima), length(SH.saddles)};
    
    graphs_ = graphs;
