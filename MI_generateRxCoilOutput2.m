function [Vo, Vo_max, Vo_opt_max, I_Rx_max]...
            = MI_generateRxCoilOutput2(T_STEPS, dt, N_Rx,...
                u0, radius_Rx, b_Rx, f0_desired, sample_Rx_Hz_flux, ...
                sample_Rx_Hrho_flux, R_Rx, I_Tx, Vin, input_mode, ...
                input_type, sample_Rx_opt_flux, angle_degrees, number_of_Rx)
            
    Havg = zeros(number_of_Rx, T_STEPS);
    flux = zeros(number_of_Rx, T_STEPS);
    Vemf = zeros(number_of_Rx, T_STEPS);
    Havg_opt = zeros(1, T_STEPS);
    flux_opt = zeros(1, T_STEPS);
    Vemf_opt = zeros(1, T_STEPS);
        
    % find flux of H through Rx coil as function of time
    for j = 1:number_of_Rx
        
        len = length(sample_Rx_Hz_flux(1, :, 1));
        
        for nn = 1:T_STEPS

            % average H over 2D line representing Rx coil at each timestep

            for i = 1:len
                
                Havg(j, nn) = Havg(j, nn) + sample_Rx_Hz_flux(nn, i, j)*cos(angle_degrees(j)*pi/180) + sample_Rx_Hrho_flux(nn, i, j)*(-sin(angle_degrees(j)*pi/180));
            
            end

            Havg(j, nn) = Havg(j, nn)/len;
            % approximate integral over circle with H ~ const w.r.t. space
            flux(j, nn) = pi*radius_Rx*radius_Rx*u0*Havg(j, nn);

        end    

        % find induced voltage (Vemf) using Vemf = -N*d(flux)/dt
        for nn = 2:T_STEPS-1

            Vemf(j, nn) = -N_Rx*(flux(j, nn+1) - flux(j, nn-1))/(2*dt); % using central diff (O(h^2))

        end
    
    end
    
    for nn = 1:T_STEPS
    
        % For optimal Rx coil alignment
        % average H over 2D line representing Rx coil at each timestep
        Havg_opt(nn) = sum(sample_Rx_opt_flux(nn, :))/length(sample_Rx_opt_flux(nn, :));
        % approximate integral over circle with H ~ const w.r.t. space
        flux_opt(nn) = pi*radius_Rx*radius_Rx*u0*Havg_opt(nn);
                
    end
    
    % find induced voltage (Vemf) using Vemf = -N*d(flux)/dt
    for nn = 2:T_STEPS-1

        Vemf_opt(nn) = -N_Rx*(flux_opt(nn+1) - flux_opt(nn-1))/(2*dt); % using central diff (O(h^2))

    end
    
    % input Vemf into RLC circuit
    L_Rx = N_Rx*N_Rx*u0*radius_Rx*(log(8*radius_Rx/b_Rx)-2); %Self inductance of loop from Balanis
    C_required = (1/L_Rx)*(1/(f0_desired*2*pi))^2;
    C_Rx = C_required; 

    a1 = 1/(dt*dt) + R_Rx/(dt*2*L_Rx);
    a2 = 1/(2*L_Rx*dt);
    a3 = 2/(dt*dt) - 1/(L_Rx*C_Rx);
    a4 = R_Rx/(dt*2*L_Rx) - 1/(dt*dt);

    I_Rx = zeros(number_of_Rx, T_STEPS);
    I_Rx_opt = zeros(1,T_STEPS);
    
    I_Rx_opt(1) = (a2/a1)*Vemf_opt(1); %n=0
    I_Rx_opt(2) = (a2/a1)*Vemf_opt(2) + (a3/a1)*(I_Rx_opt(1)); %n=1
    
    for j = 1:number_of_Rx
        
        I_Rx(j, 1) = (a2/a1)*Vemf(j, 1); %n=0
        I_Rx(j, 2) = (a2/a1)*Vemf(j, 2) + (a3/a1)*(I_Rx(j, 1)); %n=1

        for nn = 2:T_STEPS-1

            I_Rx(j, nn+1) = (a2/a1)*(Vemf(j, nn+1)-Vemf(j, nn-1)) + (a3/a1)*I_Rx(j, nn) + (a4/a1)*I_Rx(j, nn-1);

        end    
    
    end
    
    for nn = 2:T_STEPS-1

        I_Rx_opt(nn+1) = (a2/a1)*(Vemf_opt(nn+1)-Vemf_opt(nn-1)) + (a3/a1)*I_Rx_opt(nn) + (a4/a1)*I_Rx_opt(nn-1);

    end
    
    f0 = 1/(2*pi*sqrt(L_Rx*C_Rx));
    Q_Rx = 2*pi*f0*L_Rx/R_Rx;
    
    % find the output voltage Vo = I*R
    Vo = I_Rx*R_Rx;
    Vo_opt = I_Rx_opt*R_Rx;
    
    % generate max values of Vo, IRx, and H
    Vo_max = zeros(1, number_of_Rx);
    I_Rx_max = zeros(1, number_of_Rx);
    for j = 1:number_of_Rx
        
        Vo_max(j) = max(abs(Vo(j, ceil(length(Vo)/2):end)));0
        I_Rx_max(j) = max(abs(I_Rx(j, ceil(length(I_Rx)/2):end)));
    
    end
    
    Vo_opt_max = max(abs(Vo_opt(ceil(length(Vo_opt)/2):end)));
    
    if (input_type == 'p') % perform an FFT if input is a pulse
       
        if (input_mode == 'V')
            
            MI_FFTanalysis(Vin, Vo_opt, dt, 776, 'V')
            
        elseif (input_mode == 'I')
            
            MI_FFTanalysis(I_Tx, I_Rx_opt, dt, 776, 'I')
            
        end
        
    end
    
    nn = 1:T_STEPS;
    
    figure(2)
    subplot(2,1,2)
    hold off
    yyaxis left
    plot(nn*dt*10^6, Vemf_opt*1000, 'k-',...
        nn*dt*10^6, Vo_opt*1000, 'b--')
    yyaxis right
    plot(nn*dt*10^6, I_Rx_opt*10^6, 'r.')
    txt1 = sprintf('Generating Rx Coil Current, V_e_m_f, and V_o\r\n');
    txt2 = sprintf('C_R_x = %1.1f nF, L_R_x = %1.1f \\muH, R_R_x = %1.0f \\Omega, N_R_x = %d, Q_R_x = %1.1f', C_Rx*10^9, L_Rx*10^6, R_Rx, N_Rx, Q_Rx);
    title([txt1, txt2])
    legend('V_e_m_f', 'V_o', 'I_L_R_x', 'Location', 'NorthWest');
    txt3 = sprintf('time [\\musec]');
    xlabel(txt3)
    xlim([0 dt*T_STEPS*10^6])
    yyaxis left
    txt4 = sprintf('V_e_m_f and V_o [mV]');
    ylabel(txt4)
    ylim([-1.1*max(abs(Vemf_opt*1000)) 1.1*max(abs(Vemf_opt*1000))])
    yyaxis right
    txt5 = sprintf('I_L_R_x [\\muA]');
    ylabel(txt5)
    ylim([-1.1*max(abs(I_Rx_opt*10^6)) 1.1*max(abs(I_Rx_opt*10^6))])
    grid minor
    drawnow   
    
end