function [sample_Rx_Hz, sample_Rx_Hrho, sample_Rx_Ephi, sample_Tx_Hz,...
            sample_Tx_Hrho, sample_Tx_Ephi, sample_Rx_Hz_flux, sample_Rx_Hrho_flux,...
            Ephi_tot, Hz_tot, Hrho_tot, IEphi, Curl_H_phi, curl_E_rho, ICErho, curl_E_z, ICEz,...
            IHrho, IHz, ICHphi, diff, sample_Rx_opt_flux, vidObj]...
                = MI_initialize(T_STEPS, Nrho, Nz, do_annimation, radius_Rx, i_Rx, i_Rx_opt,...
                    number_of_Rx, drho, dx)

    % Initialize for analysis
    sample_Rx_Hz = zeros(1,T_STEPS);
    sample_Rx_Hrho = zeros(1,T_STEPS);
    sample_Rx_Ephi = zeros(1,T_STEPS);
    sample_Tx_Hz = zeros(1,T_STEPS);
    sample_Tx_Hrho = zeros(1,T_STEPS);
    sample_Tx_Ephi = zeros(1,T_STEPS);
    
    % initialize for convergence analysis
    diff = 0;
    
    % Initialize fields
    Ephi_tot = zeros(Nrho, Nz);
    Hz_tot = zeros(Nrho, Nz);
    Hrho_tot = zeros(Nrho, Nz);
    
    % for FDTD update equations
    IEphi = zeros(Nrho, Nz);
    Curl_H_phi = zeros(Nrho, Nz);
    curl_E_rho = zeros(Nrho, Nz);
    ICErho = zeros(Nrho, Nz);
    curl_E_z = zeros(Nrho, Nz);
    ICEz = zeros(Nrho, Nz);
    IHrho = zeros(Nrho, Nz);
    IHz = zeros(Nrho, Nz);
    ICHphi = zeros(Nrho, Nz);
    
    % for finding flux through Rx coil
    max_len = 0;
    
    for i = 1:number_of_Rx   

        Rx_i = (i_Rx(i)-ceil((radius_Rx-dx(i))/drho)):(i_Rx(i)+ceil((radius_Rx-dx(i))/drho));
        max_len = max(length(Rx_i), max_len);
        
    end
    
    sample_Rx_Hz_flux = zeros(T_STEPS, max_len, number_of_Rx);
    sample_Rx_Hrho_flux = zeros(T_STEPS, max_len, number_of_Rx);
    
    sample_Rx_opt_flux = zeros(T_STEPS, length(i_Rx_opt-radius_Rx/drho:i_Rx_opt+radius_Rx/drho));
    
    % initialize video
    if (do_annimation == true)
        
        vidObj = VideoWriter('FDTDsim_fullSim_change_name', 'MPEG-4');
        open(vidObj);
    
    else 
        
        vidObj = 0;
        
    end
    
end