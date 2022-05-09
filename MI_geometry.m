function [RX2_rho, RX1_rho, RX2_z, RX1_z, SRC2_rho, SRC1_rho, SRC2_z, SRC1_z,...
        lambda_res, dt, T_STEPS, Nrho, Nz, i_SRC1, j_SRC1, i_SRC2, skin_depth,...
        j_SRC2, depthInH2O_Tx, i_Rx, j_Rx, i_Tx, j_Tx, RHO, Z, PML_thickness, lambda_H2O, lambda_air,...
        drho, dz, j_H2O_height, heightFromH2O_Rx, dmin_lambda, dmin_dim,...
        i_SRC3, j_SRC3, i_SRC4, j_SRC4, i_Rx_opt, j_Rx_opt, angle_degrees, slope, dx]...
            = MI_geometry(c0, freq, u0, Eps0, sigma_H2O, Er_H2O, lambda_res, dim_res,...
            radius_Tx, num_periods, height_Tx, height_H2O, height_Rx,...
            number_of_Rx, angle_degrees_start, angle_degrees_diff,...
            height_max, center_Rx, center_Tx, width, radius_Rx, input_type,...
            FFT_resolution, use_two_Tx_coils, height_Rx_opt)

%%
% The problem geometry looks like this:
%{

              |----------------d6----------------|
  -           |-------d5--------|----------------|
  |      
  |                                                    
  |                                                     
  |  -                          Rx                    
  |  |                                                 
  |  |                                                
  |  |                                                 
  |  |  -      ------------------------------------      
  d4 |  |      ------------------------------------    
  |  |  |      ------------------------------------      
  |  d3 |      ------------------------------------      
  |  |  |      ------------------------------------    
  |  |  |  -   -----------------Tx-----------------      
  |  | d2  |   ------------------------------------       
  |  |  |  |   ------------------------------------       
  |  |  |  |   ------------------------------------       
  |  |  |  d1  ------------------------------------       
  |  |  |  |   ------------------------------------       
  |  |  |  |   ------------------------------------       
  -  -  -  -   ------------------------------------             

    Tx := Source Injection Point (transmitter)
    Rx := Receiver

    Note2:  1) Dirichlet boundary conditions from nxE=0 (cross-product) 
                and n*B=0 (dot-product) (n is the normal vector)
                ->  along i = 1/Nrho: Ephi = Hrho = 0
                    along j = 1/Nz: Ephi = Hz = 0
            2) Using surf() maps matrix to grid (i.e. Ephi(1,1) |--> origin); 
            therefore, one can think of Ephi matrix as Ephi(rho,z) @ t
            3) Simulating TM-mode (Ephi != 0, Hphi = 0)
            4) Origin is lower left corner

%}
    angle_radians = zeros(1, number_of_Rx);
    slope = zeros(1, number_of_Rx);
    dx = zeros(1, number_of_Rx);
    angle_degrees = zeros(1, number_of_Rx);
    
    for i = 1:number_of_Rx
        
        angle_degrees(i) = angle_degrees_start + (i - 1)*angle_degrees_diff;
        angle_radians(i) = angle_degrees(i)*pi/180;
        slope(i) = tan(angle_radians(i));
        dx(i) = radius_Rx*(1-cos(angle_radians(i)));
       
    end
    
    depthInH2O_Tx = height_H2O - height_Tx;
    heightFromH2O_Rx = height_Rx - height_H2O;
    
    smallest_dim = min([height_Tx, height_H2O, height_Rx, height_max, width, radius_Tx, radius_Rx*cos((angle_degrees_start+(number_of_Rx - 1)*angle_degrees_diff)*pi/180)]);
    
    [dz, drho, dt, T_STEPS, lambda_res, lambda_H2O, lambda_air, skin_depth]...
        = MI_descritizationSetup(freq, u0, c0, Er_H2O, Eps0, sigma_H2O,...
        lambda_res, smallest_dim, dim_res, depthInH2O_Tx, heightFromH2O_Rx,...
        num_periods, input_type, FFT_resolution);
    
    dmin_lambda = lambda_H2O/lambda_res;
    dmin_dim = smallest_dim/dim_res;
    
    % determine # of cells (ensures Nrho is an even #)
    if (mod(ceil(width/drho), 2) == 0)       
        Nrho = ceil(width/drho);
    else
        width = width + drho;
        Nrho = ceil((width)/drho);
    end
    Nz = ceil(height_max/dz);
    drho = width/Nrho; 
    
    % PML (in # of grid Yee cell)
    PML_thickness = floor(Nrho/3);
    
    if (PML_thickness > 30)
        
            PML_thickness = 30;
            
    end
    
    % for plotting Tx line in annimation
    SRC1_rho = radius_Tx;     % location of source [m]
    SRC1_z = height_Tx;
    SRC2_rho = -radius_Tx;      % location of source [m]
    SRC2_z = height_Tx;
    
    i_SRC1 = ceil(((width/2)+SRC1_rho)/drho);           % i/j corresponding to source
    j_SRC1 = floor(SRC1_z/dz);
    i_SRC2 = floor(((width/2)+SRC2_rho)/drho);         
    j_SRC2 = floor(SRC2_z/dz);
    
    i_Tx = ceil(center_Tx/drho);                        % i/j corresponding to transmitter
    j_Tx = ceil(height_Tx/dz);
        
    i_Rx = zeros(1, number_of_Rx);
    j_Rx = zeros(1, number_of_Rx);
    
    for i = 1:number_of_Rx
        
        i_Rx(i) = ceil(center_Rx(i)/drho);              % i/j corresponding to receiver
        j_Rx(i) = ceil(height_Rx(i)/dz);
    
    end
    
    j_Rx_opt = ceil(height_Rx_opt/dz);                  % i/j corresponding to optimally aligned receiver
    
    j_H2O_height = ceil(height_H2O/dz);
    
    % for use with surf()
    [RHO,Z] = meshgrid(drho*(-(Nrho/2)+1):drho:drho*(Nrho/2), dz:dz:dz*(Nz));
    
    grid_setup = zeros(Nrho, Nz) + 0.29;
    grid_setup(:, j_H2O_height:end) = 0.36;  
    grid_setup(i_Tx-radius_Tx/drho:i_Tx+radius_Tx/drho, j_Tx) = 0.99;
    
    RX1_rho = zeros(1, number_of_Rx);
    RX2_rho = zeros(1, number_of_Rx);
    RX1_z = zeros(1, number_of_Rx);
    RX2_z = zeros(1, number_of_Rx);
    
    for i = 1:number_of_Rx
        
        Rx_i = (i_Rx(i)-ceil((radius_Rx-dx(i))/drho)):(i_Rx(i)+ceil((radius_Rx-dx(i))/drho));
        for j = 1:1:length(Rx_i)
            y = slope(i)*(Rx_i(j) - i_Rx(i)) + j_Rx(i);
            grid_setup(Rx_i(j), ceil(y)) = 0.99;
        end
        
        RX1_rho(i) = center_Rx(i)  + radius_Rx*cos(angle_radians(i)) - width/2;     % location of receiver [m]
        RX2_rho(i) = center_Rx(i) - radius_Rx*cos(angle_radians(i)) - width/2;      % location of receiver [m]
        RX1_z(i) = height_Rx(i) + radius_Rx*sin(angle_radians(i));
        RX2_z(i) = height_Rx(i) - radius_Rx*sin(angle_radians(i));
    
    end
    
    % For optimally aligned Rx coil
    i_Rx_opt = i_Tx;
    grid_setup(i_Rx_opt-radius_Rx/drho:i_Rx_opt+radius_Rx/drho, j_Rx_opt) = 0.99;
    
    % setup second transmitter coil
    i_SRC3 = 0; j_SRC3 = 0; i_SRC4 = 0; j_SRC4 = 0;
    if (use_two_Tx_coils == true)
       
        i_SRC3 = ceil((width/2)/drho);         % i/j corresponding to second source
        j_SRC3 = floor((SRC1_z/dz)+SRC1_rho);
        i_SRC4 = floor((width/2)/drho);         
        j_SRC4 = floor((SRC2_z/dz)+SRC2_rho);
        
        grid_setup(i_Tx, j_Tx-radius_Tx/drho:j_Tx+radius_Tx/drho) = 0.99;
        
    end
    
    fig1 = figure(1);
    set(fig1, 'Name', 'Debug - Grid and Coil Setup', 'NumberTitle','off');
    movegui('west');
    axis equal
    axis tight
    cla
    hold on
    surf(RHO, Z, grid_setup', 'FaceAlpha', 0.7);
    view(2);
    colormap jet;
    caxis([0 1]);
    title('The Grid Setup (Rx coil on top, Tx on bottom)')
    xlabel('\rho-position [m]')
    ylabel('z-position [m]')
    drawnow
    
end