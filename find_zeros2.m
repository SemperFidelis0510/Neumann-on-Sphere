function Z = find_zeros2(M)     
    move_horizontal = M(1:end-1,2:end).*M(1:end-1,1:end-1);    
    move_vertical = M(2:end,1:end-1).*M(1:end-1,1:end-1);
    Z = (move_vertical<=0)|(move_horizontal<=0);