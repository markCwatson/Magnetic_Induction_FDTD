function MI_Validation(I_max, N_Tx, N_Rx, Vo_max, Vo_opt_max, radius_Tx, u0, sigma_H2O,...
            depthInH2O_Tx, heightFromH2O_Rx, freq, fmax, fmin)

    %from Gibson
    [Vo_Gibson, Hz_Gibson, f_Gibson] = MI_Gibson(I_max, N_Tx, N_Rx, radius_Tx, u0, depthInH2O_Tx, heightFromH2O_Rx, sigma_H2O);

    figure('Name','Output - Validation','NumberTitle','off');     
    semilogx(f_Gibson(2:end), 20*log10(Vo_Gibson(2:end)), 'k-')
    hold on
    semilogx(freq, 20*log10(Vo_opt_max), 'rd', freq, 20*log10(Vo_max), 'ro')
    txt = sprintf('Using Wait''s Sommerfeld Integral');
    title(txt)
    xlabel('frequency [Hz]')
    ylabel('20log(|Vo|) [dBV]')
    xlim([fmin/10 10*fmax])
    grid on
    
    if (max(20*log10(Vo_Gibson)) >= max(20*log10(Vo_max)))
        
        ylim([1.2*min(20*log10(Vo_Gibson)) 0.8*max(20*log10(Vo_Gibson))])
        
    else
        
        ylim([1.2*min(20*log10(Vo_Gibson)) 0.8*max(20*log10(Vo_max))])
        
    end
    
    legend('Gibson (Aligned)', 'FDTD (Aligned)', 'FDTD (Misaligned)','Location', 'NorthEast')
    
end