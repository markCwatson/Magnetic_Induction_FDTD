function MI_plotHcomponents(sample_Tx_Hz, sample_Tx_Hrho, sample_Rx_Hz, sample_Rx_Hrho, dt, T_STEPS)
    
    fig7 = figure(7);
    set(fig7, 'Name', 'Output - Tx and Rx Coil H Fields', 'NumberTitle','off');
    subplot(2,1,1)
    cla
    hold on
    t = 0:dt*10^6:dt*(T_STEPS-1)*10^6;
    plot(t, sample_Tx_Hz, 'r-.', t, sample_Tx_Hrho, 'b--')
    title('H_z(t) and H_\rho(t) at Center of Tx Coil')
    legend('H_z(t)', 'H_\rho(t)', 'Location', 'northwest')
    xlabel('time [\musec]')
    ylabel('H-Fields [A/m]')
    xlim([0 max(t)])
    grid minor
    
    hold on
    subplot(2,1,2)
    plot(t, sample_Rx_Hz*1000, 'r-.', t, sample_Rx_Hrho*1000, 'b--')
    title('H_z(t) and H_\rho(t) at Center of Rx Coil')
    legend('H_z(t)', 'H_\rho(t)', 'Location', 'northwest')
    xlabel('time [\musec]')
    ylabel('H-Fields [mA/m]')
    xlim([0 max(t)])
    grid minor
        
end