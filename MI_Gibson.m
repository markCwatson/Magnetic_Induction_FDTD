function [Vo_Gibson, Hz_Gibson, f_Gibson]  = MI_Gibson(I_max, N_Tx, N_Rx, radius, u0, depth, height, sigma)

    f_Gibson = [1000:1000:90000 100000:100000:900000 1000000:1000000:10000000];
    depth = [depth/2 depth depth*2 depth*5];
    line_styles = ["k-" "k--" "k:" "k-."];
    mag_Hz = zeros(length(depth), length(f_Gibson));
    mag_V = zeros(length(depth), length(f_Gibson));
    
    fig33 = figure(33);
    set(fig33, 'Name', 'Output - Wait''s Sommerfeld Integral: V_e_m_f', 'NumberTitle','off');
    fig34 = figure(34);
    set(fig34, 'Name', 'Output - Wait''s Sommerfeld Integral: H', 'NumberTitle','off');
    
    for j = 1:length(depth)

        for i = 1:length(f_Gibson)

            % myPQintegral(radius_Tx, freq, sigma, mu, depth, height, offset)
            [H_rho, H_z] = MI_myPQintegral(I_max, N_Tx, radius, f_Gibson(i), sigma, (1e-7)*4*pi, depth(j), height, 0);
            mag_Hz(j, i) = abs(H_z);

        end

        mag_V(j,:) = N_Rx*2*pi*pi*radius*radius*u0*f_Gibson.*mag_Hz(j,:);
        
        figure(33)
        loglog(f_Gibson, 20*log10(mag_V(j, :)), line_styles(j), 'LineWidth', 2);
        hold on
        figure(34)
        loglog(f_Gibson, 20*log10(mag_Hz(j, :)), line_styles(j), 'LineWidth', 2);
        hold on

    end
    
    % the output in row 2 is always at the desired depth
    Vo_Gibson = mag_V(2,:);
    Hz_Gibson = mag_Hz(2,:);
    
    leg1 = num2str(depth(1));
    leg1 = join([leg1 "m"]);
    leg2 = num2str(depth(2));
    leg2 = join([leg2 "m"]);
    leg3 = num2str(depth(3));
    leg3 = join([leg3 "m"]);
    leg4 = num2str(depth(4));
    leg4 = join([leg4 "m"]);
    
    figure(33)
    grid on
    legend([ leg1 leg2 leg3 leg4], 'Location', 'SouthWest')
    xlim([min(f_Gibson) max(f_Gibson)])
    ylim([10*min(20*log10(Vo_Gibson)) 0.5*max(20*log10(Vo_Gibson))])
    xlabel('frequency [Hz]')
    ylabel('20log(|V_e_m_f|) [dBV]')
    txt1 = sprintf('Induced Voltage at Varying Depths Using Wait');
    txt2 = sprintf('(I_T_x = %1.1f A, \\sigma_H_2_O = %1.1f S/m, h_R_x = %1.1f m, N_T_x = %d, N_R_x = %d, r = %1.1f cm)\r\n', I_max, sigma, height, N_Tx, N_Rx, 100*radius);
    title({txt1, txt2});
    
    figure(34)
    grid on
    legend([ leg1 leg2 leg3 leg4], 'Location', 'SouthWest')
    xlim([min(f_Gibson) max(f_Gibson)])
    ylim([10*min(20*log10(Hz_Gibson)) 0.5*max(20*log10(Hz_Gibson))])
    xlabel('frequency [Hz]')
    ylabel('20log(|H|) [dBA/m]')
    txt3 = sprintf('Low-Frequency Window at Varying Depths Using Wait');
    txt4 = sprintf('(I_T_x = %1.1f A, \\sigma_H_2_O = %1.1f S/m, h_R_x = %1.1f m, N_T_x = %d, N_R_x = %d, r = %1.1f cm)\r\n', I_max, sigma, height, N_Tx, N_Rx, 100*radius);
    title({txt3, txt4});

end

% fin