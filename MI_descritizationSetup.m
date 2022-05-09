function  [dz, drho, dt, T_STEPS, lambda_res, lambda_H2O, lambda_air, skin_depth] =...
    MI_descritizationSetup(freq, u0, c0, Er_H2O, Eps0, sigma_H2O, lambda_res,...
    smallest_dim, dim_res, depthInH2O_Tx, heightFromH2O_Rx, num_periods, input_type, FFT_resolution)
    
    % gamma := complex propagation constant
    % in general, gamma = alpha + i*beta = sqrt(i*2*pi*f*U0*(sigma + i*2*pi*f*epsilon))
    % alpha = 2*pi*f*sqrt(U0*epsilon/2)*sqrt(sqrt(1 + (sigma/(2*pi*f*epsilon))^2) - 1)
    % beta = 2*pi*f*sqrt(U0*epsilon/2)*sqrt(sqrt(1 + (sigma/(2*pi*f*epsilon))^2) + 1)
    % C_H2O = f*lambda_H2O (speed of wave in water)
    beta = 2*pi*freq*sqrt(u0*Er_H2O*Eps0/2)*sqrt(sqrt(1+(sigma_H2O/(2*pi*freq*Er_H2O*Eps0))^2)+1);
    lambda_H2O = 2*pi/beta;
    lambda_air = c0/freq;
    
    alpha = 2*pi*freq*sqrt(u0*Er_H2O*Eps0/2)*sqrt(sqrt(1+(sigma_H2O/(2*pi*freq*Er_H2O*Eps0))^2)-1);
    skin_depth = 1/alpha; % attenuates to 1/e
    
    vp_H2O = 2*pi*freq/beta;
    T_H2O = depthInH2O_Tx/vp_H2O;
    T_air = heightFromH2O_Rx/c0;
    T_tot = T_air + T_H2O; % time for wave to travel from Tx to Rx
    
    drho = min([lambda_H2O/lambda_res, smallest_dim/dim_res]);
    dz = drho;
    
    dt = 0.9/(c0*sqrt((1/drho)^2+(1/dz)^2)); % Courant stability criteria
    fs = 1/dt; % sampling freq
    
    T_period = 1/freq;
    time_steps_per_period = T_period/dt;
    
    if (input_type == 't')
        
        T_STEPS = ceil(num_periods*time_steps_per_period);
        
    elseif (input_type == 'p')
        
        % note: The frequency resolution is equal to the sampling...
        %... frequency divided by # bins
        T_STEPS = ceil(fs/FFT_resolution);
       
    end

end