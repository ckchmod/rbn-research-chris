function [energy] = energy(s_i , N_S_W_E, J, h_interac, h_f)

    N_S_W_E_temp = N_S_W_E;
    s_i_temp = s_i;
    h_f_temp = h_f;

    for i = 1:length(N_S_W_E_temp)
        
        if (N_S_W_E_temp(i) == 0)
            N_S_W_E_temp(i) = -1;
        end
          
    end
    
    if (s_i_temp == 0)
        s_i_temp = -1;
    end
    
    if (h_f_temp == 0)
        h_f_temp = -1;
    end
    
    energy = -1*J * s_i_temp* (N_S_W_E_temp(1) + N_S_W_E_temp(2) + N_S_W_E_temp(3) + N_S_W_E_temp(4)) ...
        - h_interac*s_i_temp*h_f_temp; 
    
end

