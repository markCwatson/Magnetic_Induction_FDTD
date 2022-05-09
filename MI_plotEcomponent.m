function MI_plotEcomponent(sample_Tx_Ephi, sample_Rx_Ephi, dt, T_STEPS)
    
    fig8 = figure(8);
    set(fig8, 'Name', 'Output - Tx and Rx Coil H Fields', 'NumberTitle','off');
    subplot(2,1,1)
    cla
    hold on
    t = 0:dt*10^6:dt*(T_STEPS-1)*10^6;
    plot(t, sample_Tx_Ephi, 'r-.')
    title('E_\phi(t) at Center of Tx Coil')
    xlabel('time [\musec]')
    ylabel('E_\phi [V/m]')
    xlim([0 max(t)])
    grid minor
    
    hold on
    subplot(2,1,2)
    plot(t, sample_Rx_Ephi*1000, 'r-.')
    title('E_\phi(t) at Center of Rx Coil')
    xlabel('time [\musec]')
    ylabel('E_\phi [mV/m]')
    xlim([0 max(t)])
    grid minor
        
end