function SH_ = get_SH(coeffs, coo, wb)
    if wb~=0
        waitbar(0.05,wb,'Coordinate system acquired. Getting function...');
    end
    
    %creates the object of the trial function
    SH = SphericalHarmonic(coeffs);
    
    %gets function handle, gradient function and hessian
    SH = SH.get_function(wb);
    
    %gets function values on the sphere
    if wb~=0
        waitbar(0.5, wb, 'Function acquired. Getting numerical data...');
    end
    SH = SH.get_values(coo.grid);
    
    %get nodal lines
    if wb~=0
        waitbar(0.6, wb, 'Getting nodal lines...');
    end
    SH = SH.get_nodal(coo.grid);
    
    %get extrema and saddles
    if wb~=0
        waitbar(0.6, wb, 'Nodal lines acquired. Getting extrema and saddles...');
    end
    SH = SH.get_extrema(coo);
    
    %get neumann lines
    if wb~=0
        waitbar(0.6, wb, 'Extrema and saddles acquired. Getting Neumann lines...');
    end
    SH = SH.get_neumann_lines(wb);
    SH = SH.merge_neumann_lines; 
    
    if ~isempty(SH.maxima)||~isempty(SH.minima)||~isempty(SH.saddles)
        SH = SH.get_extrema_deg;
    end
    
    SH_ = SH;
end