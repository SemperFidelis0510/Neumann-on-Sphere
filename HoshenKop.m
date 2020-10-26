function [colors_n, colors_, M] = HoshenKop(N)
        %Only fit for sphere. Periodic on the columns (phi axis).       
        S = size(N);
        N = [zeros(1,S(2)+3); zeros(S(1),1), N, N(:,1), zeros(S(1),1)];       
        [cols, rows] = find(N');
        colors = [];
        
        for i=1:length(rows)
                n = rows(i);
                m = cols(i);
                c = 0;
                flag = 0;
                
                if N(n-1,m)
                        c = N(n-1,m);                        
                elseif N(n-1,m+1)
                        c = N(n-1,m+1);
                        flag = 1;
                elseif N(n-1,m-1)
                        c = N(n-1,m-1);
                elseif N(n,m-1)
                        c = N(n,m-1);
                end
                
                if flag
                        if N(n-1,m-1)
                                colors = merge_colors(colors, colors(N(n-1,m-1)), colors(c));
                        elseif N(n,m-1)
                                colors = merge_colors(colors, colors(N(n,m-1)), colors(c));
                        end
                end              
                
                if ~c
                        c = length(colors)+1;
                        colors = [colors, c];
                end
                
                N(n,m) = c;
        end
        
        for i=2:S(1)
                if N(i,2)
                        colors = merge_colors(colors, colors(N(i,end-1)), colors(N(i,2)));
                end
        end
        
        for i=1:length(colors)
                N(N==i) = colors(i);
        end
        
        colors_n=length(unique(colors));
        colors_=unique(colors);
        M = N(2:end,2:end-2);
        
        
end

function c_ = find_color(colors, c)
        while colors(c)~=c
                c = colors(c);
        end
        c_ = c;
end

function V_ = merge_colors(V, old_col, new_col)         
        V(V==old_col) = new_col;
        V_ = V;
end