function MI_printParams(freq, dt, T_STEPS, drho, sigma_H2O, lambda_H2O, skin_depth, lambda_air,...
            I_max, N_Tx, radius_Tx, depthInH2O_Tx, N_Rx, radius_Rx, heightFromH2O_Rx)
    
    fig100 = figure(100);
    set(fig100, 'Name', 'Debug - Important Parameters', 'NumberTitle','off');
    ah100 = axes(fig100,'position',[.2,.2,.6,.6]); 
    movegui('west');
    % timing/discretization
    message1 = sprintf('f = %1.1f kHz \r\n\\Deltat = %1.2f ps \r\nSTEPS = %d \r\n\\Delta\\rho = \\Deltaz = %1.2f mm\r\n',...
        freq/1000, dt*10^(12), T_STEPS, drho*1000);
    % EM effects
    message2 = sprintf('\\sigma_H_2_O = %1.1f S/m \r\n\\lambda_H_2_O = %1.1f m \r\n\\delta = %1.2f m \r\n\\lambda_0 = %1.1f m \r\n',...
        sigma_H2O, lambda_H2O, skin_depth, lambda_air);
    % Tx coil
    message3 = sprintf('I_m_a_x = %1.1f A \r\nN_T_x = %d turn(s) \r\nr_T_x = %1.1f cm \r\nd_T_x = %1.1f m \r\n',...
        I_max, N_Tx, 100*radius_Tx, depthInH2O_Tx);
    % Rx coil
    message4 = sprintf('N_R_x = %d turn(s) \r\nr_R_x = %1.1f cm \r\nh_R_x = %1.1f m \r\n',...
        N_Rx, 100*radius_Rx, heightFromH2O_Rx);
    text(ah100,.5,.5,[message1 message2 message3 message4],'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
    ah100.Visible = 'off';
     
end

%fin