function [H_Tx, Vin, I_Tx]...
            = MI_generateTxCoilOutput(input_type, input_mode,...
                I_max, Vin_Amp, T_STEPS, dt, freq, N_Tx, u0, radius_Tx, b_Tx,...
                f0_desired, R_Tx, sigma_wire, t0, tau, fmax_pulse, drho)

    L = N_Tx*N_Tx*u0*radius_Tx*(log(8*radius_Tx/b_Tx)-2); %Self inductance of loop from Balanis
    %note: if L is too big, the system response is very slow

    C_required = (1/L)*(1/(f0_desired*2*pi))^2;
    C_Tx = C_required;
    
    f0 = 1/(2*pi*sqrt(L*C_Tx));
    Q = 2*pi*f0*L/R_Tx;

    nn = 1:T_STEPS;
    Vin = 0;
    
    if (input_type == 't') % single tone
    
        if (input_mode == 'V') % input voltage

            %sine input
            Vin = Vin_Amp*sin((2*pi*freq*nn*dt));
            %sine with envelope
            %tau = 2*(1/freq);
            %Vin = vin_Amp*sin((2*pi*freq*nn*dt)).*(1-exp(-(nn-1)*dt/tau));
         
            % I might replace this with Runge Kutta (RK4)
            a1 = 1/(dt*dt) + R_Tx/(dt*2*L);
            a2 = 1/(2*L*dt);
            a3 = 2/(dt*dt) - 1/(L*C_Tx);
            a4 = R_Tx/(dt*2*L) - 1/(dt*dt);

            I_Tx = zeros(1,T_STEPS);
            I_Tx(1) = (a2/a1)*Vin(1); %n=0
            I_Tx(2) = (a2/a1)*Vin(2) + (a3/a1)*(I_Tx(1)); %n=1
            
            for n = 2:T_STEPS-1

                I_Tx(n+1) = (a2/a1)*(Vin(n+1)-Vin(n-1)) + (a3/a1)*I_Tx(n) + (a4/a1)*I_Tx(n-1);

            end

        elseif (input_mode == 'I') % input current

            I_Tx = I_max.*sin((2*pi*freq*nn*dt));

        end
        
    elseif (input_type == 'p') % pulse input

        if (input_mode == 'V') % input voltage

            %voltage pulse
            Vin = Vin_Amp*exp(-((nn*dt - t0)/tau).^2);
            
            a1 = 1/(dt*dt) + R_Tx/(dt*2*L);
            a2 = 1/(2*L*dt);
            a3 = 2/(dt*dt) - 1/(L*C_Tx);
            a4 = R_Tx/(dt*2*L) - 1/(dt*dt);

            I_Tx = zeros(1,T_STEPS);
            I_Tx(1) = (a2/a1)*Vin(1); %n=0
            I_Tx(2) = (a2/a1)*Vin(2) + (a3/a1)*(I_Tx(1)); %n=1
            
            for n = 2:T_STEPS-1

                I_Tx(n+1) = (a2/a1)*(Vin(n+1)-Vin(n-1)) + (a3/a1)*I_Tx(n) + (a4/a1)*I_Tx(n-1);

            end

        elseif (input_mode == 'I') % input current
            
            % current pulse
            I_Tx = I_max.*exp(-((nn*dt - t0)/tau).^2);
            
        end
        
    elseif (input_type == 'm') % modulated using PAM with carrier (does not work right now)
        
        N = num_symbols; % number of symbols
        M = 4; % Alphabet size
        x = randi([0 M-1],N,1); % make symbols using 2 bit blocks (0 to 3)|dec
        y = real(pammod(x,M)); % map symbols to PAM signal: 3 -> M-1=3, 2 -> 1, 1 -> -1, 0 -> -3
        up = 25;
        L = 5; % in terms of symbols 
        a = 1;
        % h = rcosfir(a,L,up,up,'sqrt');
        h = rcosfir(a,L,up,up,'normal');

        w = upsample(y,up); % inserts up-1 zeros between points
        
        if (input_mode == 'V') % input voltage
            
            Vin = Vin_Amp*conv(h,w);
            
            a1 = 1/(dt*dt) + R_Tx/(dt*2*L);
            a2 = 1/(2*L*dt);
            a3 = 2/(dt*dt) - 1/(L*C_Tx);
            a4 = R_Tx/(dt*2*L) - 1/(dt*dt);

            I_Tx = zeros(1,T_STEPS);
            I_Tx(1) = (a2/a1)*Vin(1); %n=0
            I_Tx(2) = (a2/a1)*Vin(2) + (a3/a1)*(I_Tx(1)); %n=1
            
            for n = 2:T_STEPS-1

                I_Tx(n+1) = (a2/a1)*(Vin(n+1)-Vin(n-1)) + (a3/a1)*I_Tx(n) + (a4/a1)*I_Tx(n-1);

            end
            
        elseif (input_mode == 'I') % input current
            
            I_Tx = I_max*conv(h,w);
            
        end        
        
        figure (904)
        plot(Vin)
        title('Modulated Input Signal (Voltage)')
 
        nn = 1:length(Vin);
        T_STEPS = length(Vin);
        
    end
    
    fig2 = figure(2);
    set(fig2, 'Name', 'Output - Tx and Rx Coil Voltages/Currents', 'NumberTitle','off');
    subplot(2,1,1)
    
    if (input_mode == 'V')
        
        yyaxis left
        plot(nn*dt*10^6, Vin, 'b.')
        ylabel('V_i_n [V]')
        ylim([-1.1*max(abs(Vin)) 1.1*max(abs(Vin))])
        yyaxis right
        plot(nn*dt*10^6, I_Tx, 'r-')
        ylabel('I_L_T_x [A]')
        ylim([-1.1*max(abs(I_Tx)) 1.1*max(abs(I_Tx))])
        legend('V_i_n', 'I_L_T_x', 'Location', 'SouthWest');
        
    elseif (input_mode == 'I')
        
        plot(nn*dt*10^6, I_Tx, 'r.')
        ylabel('I_L_T_x [A]')
        ylim([-1.1*max(abs(I_Tx)) 1.1*max(abs(I_Tx))])
        
    end
    
    txt1 = sprintf('Generating Tx Coil Current\r\n');
    
    if (input_type == 'p')
    
        txt2 = sprintf('f_m_a_x ~ %1.1f MHz, L_T_x = %1.1f \\muH, N_T_x = %d, Q_T_x = %1.1f', fmax_pulse/10^6, L*10^6, N_Tx, Q);
    
    else
        
        txt2 = sprintf('f = %1.1f kHz, L_T_x = %1.1f \\muH, N_T_x = %d, Q_T_x = %1.1f', freq/10^3, L*10^6, N_Tx, Q);
    
    end
    
    title([txt1, txt2])
    txt3 = sprintf('time [\\musec]');
    xlabel(txt3)
    xlim([0 dt*T_STEPS*10^6])
    grid minor
    drawnow
    
    % since A_effective < A < Yee cell
    H_Tx = N_Tx*I_Tx/(2*pi*drho/2);
    
%     % Generate E using Ohm's law J = sigma*E, where J = I/A
%     skinDepth = sqrt(1/(sigma_wire*pi*freq*u0));
%     if (skinDepth <= b_Tx)
%         A_effective = pi*(b_Tx)^2 - pi*(b_Tx-skinDepth)^2;
%     else
%         A_effective = pi*(b_Tx)^2;
%     end
%     
%     Jphi = I_Tx/(A_effective);          % [A/m^2 per turn]
%     Ephi_Tx = N_Tx*Jphi/sigma_wire;     % [A/m^2 / (A/V)/m = V/m]
    
end