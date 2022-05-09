function [H_row, H_z] = MI_myPQintegral(I, N_Tx, radius_Tx, freq, sigma, mu, depth, height, offset)
    
    % loop antennas are horizontal and coaxial

    % PQintegral performs Sommerfeld integral
    % mode = 1 returns P for uplink direction (row component)
    % mode = 2 returns P for downlink direction (row component)
    % mode = 3 returns Q for either direction (z component)

    % Z is normalized height of observation point Z = height/depth
    % D is normalized offset of observation point D = offset/depth
    % T is normalized depth of transmitter T = sqrt(2)*depth/skin_depth
    % .. where skin_depth = sqrt(1/(pi*mu*freq*sigma))
    
    Z = height/depth;
    D = offset/depth;
    skin_depth = sqrt(1/(pi*mu*freq*sigma));
    T = sqrt(2)*depth/skin_depth;
    
    TOL = 1e-5; % defult tolerance is 1e-6 (made bigger to avoid signularity warning)
    TRACE = 0;
    XLIMIT = 30;
    
    for mode = 1:2:3 % only interested in uplink
        
        % quadl numerically evaluates Sommerfled integral using...
        % ...adaptive Gauss/Lobatto quadrature
        a = quadl(@MI_sommerfeld, 0, XLIMIT, TOL, TRACE, D, T, Z, mode);
        %a = integral(@sommerfeld, 0, XLIMIT, TOL, TRACE, D, T, Z, mode);
        %a = integral(@sommerfeld, 0, XLIMIT, 'RelTol', TOL, 'AbsTol', 1e-13, D, T, Z, mode);
            % for oscillatory f'n -> integral(fun,0,Inf,'RelTol',1e-8,'AbsTol',1e-13)
        
        A = pi*radius_Tx^2;
        mag_dipole = N_Tx*I*A; % [Am^2]

        if mode == 1
            H_row = a*mag_dipole/(2*pi*depth^3);
        else
            H_z = a*mag_dipole/(2*pi*depth^3);
        end
    
    end
    
end