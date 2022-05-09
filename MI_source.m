function [Hz_tot, Hrho_tot]...
            = MI_source(i_SRC1, j_SRC1, H_Tx, Hz_tot, Hrho_tot, n,...
                i_SRC3, j_SRC3, i_SRC4, j_SRC4, use_two_Tx_coils)
    
    % SRC1 (right) and SRC2 (due to cylindrical symmetry)
    Hz_tot(i_SRC1, j_SRC1) = H_Tx(n);
    Hz_tot(i_SRC1-1, j_SRC1) = -H_Tx(n);
    Hrho_tot(i_SRC1, j_SRC1) = -H_Tx(n);
    Hrho_tot(i_SRC1, j_SRC1-1) = H_Tx(n);
    
    if (use_two_Tx_coils)
        
        % SRC3 (top)
        Hz_tot(i_SRC3, j_SRC3) = -H_Tx(n);
        Hrho_tot(i_SRC3, j_SRC3) = H_Tx(n);
        Hrho_tot(i_SRC3, j_SRC3-1) = -H_Tx(n);
        
        % SRC4 (bottom)
        Hz_tot(i_SRC4, j_SRC4) = H_Tx(n);
        Hrho_tot(i_SRC4, j_SRC4) = -H_Tx(n);
        Hrho_tot(i_SRC4, j_SRC4-1) = H_Tx(n);
        
    end
    
end