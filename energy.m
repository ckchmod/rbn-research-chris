function [N_S_W_E] = energy(N_S_W_E)

    for i = 1:length(N_S_W_E)
        
        if (N_S_W_E(i) == 0)
            N_S_W_E(i) = -1;
        end
          
    end

end

