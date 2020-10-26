function Z = find_zeros(M)     
    S = size(M);
    move_horizontal = circshift(M,1,2).*M;    
    move_vertical = [M(2:end,:).*M(1:end-1,:); ones(1,S(2))];
    Z = (move_vertical<=0)|(move_horizontal<=0);