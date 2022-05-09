function [number_of_Rx, angle_degrees_start, angle_degrees_diff, c0, Eps0,...
            u0, freq, sigma_H2O, sigma_air, Er_H2O, Er_air, lambda_res,...
            convergence_iteration, max_convergence_iterations, dim_res, t0, tau,...
            do_convergence_iterations, threshold, input_mode, I_max, do_PML,...
            Vin_Amp, N_Tx, radius_Tx, b_Tx, f0_desired, R_Tx, C_Tx, sigma_wire, num_periods,...
            height_Tx, height_H2O, height_Rx, height_max, center_Rx, center_Tx, width, increase,...
            do_annimation, fmax, fmin, plot_fields, plot_fft, input_type,...
            N_Rx, radius_Rx, b_Rx, C_Rx, R_Rx, converg_test_type, scale_high, scale_low,...
            fmax_pulse, FFT_resolution, iterations_until_getframe, do_initial_validity_test,...
            do_initial_validity_test_annimation, do_initial_validity_test_video, last_output_max,...
            use_two_Tx_coils, height_Rx_opt]...
                = MI_setup()

    % signal parameters
    freq = 1000*10^3;            % frequency of signal [Hz]
    c0 = 299792458;             % [m/s]   
    fmax = 20*10^6;             % max freq interested in [Hz]
    fmin = 10*10^3;             % min freq interested in [Hz]
    fmax_pulse = 20*10^6;       % max freq in pulse [Hz]
    tau = 2/(pi*fmax_pulse);    %Gaussian pulse width = 2*Tau (seconds)          
    t0 = 50*tau;                %Gaussian pulse delay (seconds)
    input_mode = 'I';           % enter 'V' for voltage or 'I' for current
    I_max = 1;                  % [A]
    Vin_Amp = 10;               % [V]
    input_type = 't';           % enter 't' for a single tone defined by freq or...
                                % 'p' for a broadband pulse or...
                                % 'm' for a modulated signal (using PAM with a carrier defined by freq)
                                % note: 'm' is not functional at this time!
    
    % geometry
    height_max = 4;             % top of domain
    width = 4;                  % width of domain [m]
    height_H2O = 2;             % height of H2O line [m]
    
    height_Tx = 1.5;              % height of Tx coil [m]
    center_Tx = 2;              % center of Tx coil [m]
    
    number_of_Rx = 2;           % placed on a circle (not including the optimally placed one)
    distance_Rx_from_Tx = 1;    % radius of circle
    angle_degrees_start = 15;    % angle of first Rx coil WRT horizon
    angle_degrees_diff = 15;
    height_Rx_opt = height_Tx + distance_Rx_from_Tx; % center always aligned with Tx
    
    % generate Rx coils
    height_Rx = zeros(1, number_of_Rx); % height of Rx coil [m]
    center_Rx = zeros(1, number_of_Rx); % center of Rx coil [m]
    angle = angle_degrees_start;
    
    for i = 1:number_of_Rx
       
        height_Rx(i) = height_Tx+distance_Rx_from_Tx*sin((90-angle)*pi/180);
        center_Rx(i) = center_Tx-distance_Rx_from_Tx*cos((90-angle)*pi/180);
        angle = angle + angle_degrees_diff;
    
    end
    
    % material parameters
    Eps0 = (1e-9)/(36*pi);      % [Farad/m]
    u0 = (1e-7)*4*pi;           % [Henry/m]
    sigma_H2O = 4.0;            % conductivity in [S/m]
    sigma_air = 0;
    Er_H2O = 81;                % relative permitivity    
    Er_air = 1;  
    do_PML = true;
    
    % runtime parameters
    convergence_iteration = 1;          % initialized to 1 (don't change)
    last_output_max = 0;                % for convergence testing
    max_convergence_iterations = 1;
    do_convergence_iterations = true;   % keep 'true' to enable output (bug)
    do_initial_validity_test = false;    % injects pulse and compares with theoretical freq modes
    do_initial_validity_test_annimation = false;    % show annimation
    do_initial_validity_test_video = false;         % create video file
    lambda_res = 10;
    threshold = 1;                      % threshold for convergence testing [percent]
    dim_res = 5;                % resolution of the smallest physical dimension (radius of loop)
    increase = 1;               % for convergence testing where dim_res (new) = dim_res (old) + increase
    converg_test_type = 'V';    % 'V' or 'I' for voltage or current based convergence testing
    FFT_resolution = 100*10^3;  % desired resolutiong of FFT (effects T_STEPS if input_type = 'p')
    num_periods = 3;            % number of periods of input signal to show (increases run time!)
    use_two_Tx_coils = false;   % NOT IMPLEMENTED
    
    % annimation
    do_annimation = false;          % to show animation 
    scale_high = 5;                 % set colorbar scale for animation [A/m]
    scale_low = 250;                % is set as caxis [mA/m]
    iterations_until_getframe = 20; % number of iterations before a frame is captured
    
    % Tx coil parameters
    N_Tx = 5;                       % [turns]
    radius_Tx = 0.1;                % [m]
    b_Tx = 0.0002;                  % [m] realistic wire radius for calculation of L
                                    % 18 AWG has diameter ~ 1mm (rated to 10 - 16 A)
                                    % 20 AWG has diameter ~ 0.8mm (rated to 5 - 11 A)
    f0_desired = freq;              % [Hz]
    R_Tx = 4;                       % [Ohm]
    C_Tx = 9.4*10^-9;               % [F]
    sigma_wire = 58.14*10^6;        % [S/m] copper -> ~58.14*10^6 S/m
    
    % Rx coil parameters
    N_Rx = 5;                       % [turns]
    radius_Rx = 0.1;                % [m]
    b_Rx = b_Tx;                    % [m]
    C_Rx = 9.4*10^-9;               % [F]
    R_Rx = 50;                      % [Ohm]
    
    %plot fields
    plot_fields = false;              % enter true/false for plot H/E fields option
    
    %plot FFT
    plot_fft = false;                 % enter 'true/false for FFT analysis option
    
end