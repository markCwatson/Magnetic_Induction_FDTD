function [sample_Rx_Hz, sample_Rx_Hrho, sample_Rx_Ephi,...
        sample_Tx_Hz, sample_Tx_Hrho, sample_Tx_Ephi,...
        sample_Rx_Hz_flux, sample_Rx_Hrho_flux, sample_Rx_opt_flux]...
            = MI_sample(Hz_tot, Hrho_tot, Ephi_tot, n,...
                i_Rx, j_Rx, i_Tx, j_Tx, radius_Rx,...
                sample_Rx_Hz, sample_Rx_Hrho, sample_Rx_Ephi,...
                sample_Tx_Hz, sample_Tx_Hrho, sample_Tx_Ephi,...
                sample_Rx_Hz_flux, sample_Rx_Hrho_flux, drho, ...
                sample_Rx_opt_flux, number_of_Rx, i_Rx_opt, j_Rx_opt,...
                slope, dx)

    %sampling H/E at Rx location for FFT analysis
    sample_Rx_Hz(n) = Hz_tot(i_Rx_opt, j_Rx_opt);
    sample_Rx_Hrho(n) = Hrho_tot(i_Rx_opt, j_Rx_opt);
    sample_Rx_Ephi(n) = Ephi_tot(i_Rx_opt, j_Rx_opt);
    
    %sampling H/E at Tx location for FFT analysis
    sample_Tx_Hz(n) = Hz_tot(i_Tx, j_Tx);
    sample_Tx_Hrho(n) = Hrho_tot(i_Tx, j_Tx);
    sample_Tx_Ephi(n) = Ephi_tot(i_Tx, j_Tx);
    
    %sampling H at Rx location for calculation of flux for Vo
    for i = 1:number_of_Rx

        Rx_i = (i_Rx(i)-ceil((radius_Rx-dx(i))/drho)):(i_Rx(i)+ceil((radius_Rx-dx(i))/drho));
        
        for j = 1:1:length(Rx_i)
            
            y = slope(i)*(Rx_i(j) - i_Rx(i)) + j_Rx(i);
            sample_Rx_Hz_flux(n, j, i) = Hz_tot(Rx_i(j), ceil(y));
            sample_Rx_Hrho_flux(n, j, i) = Hrho_tot(Rx_i(j), ceil(y));
            
        end
        
    end
    
    % for comparison with Rx coil optimally aligned with Tx coil
    if (length(i_Rx_opt-floor(radius_Rx/drho):i_Rx_opt+ceil(radius_Rx/drho))== length(sample_Rx_opt_flux(n, :)))
        
        sample_Rx_opt_flux(n, :) = Hz_tot(i_Rx_opt-floor(radius_Rx/drho):i_Rx_opt+ceil(radius_Rx/drho), j_Rx_opt);
        
    else
        
        sample_Rx_opt_flux(n, :) = Hz_tot(i_Rx_opt-floor(radius_Rx/drho):i_Rx_opt+floor(radius_Rx/drho), j_Rx_opt);
    
    end
    
end