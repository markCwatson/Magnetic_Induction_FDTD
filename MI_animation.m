function MI_animation(n, T_STEPS, RHO, Z, Hz_tot, Hrho_tot, dt,...
            drho, Nrho, dz, Nz, height_H2O, sigma_H2O,...
            vidObj, freq, angle_degrees_start, angle_degrees_diff,...
            scale_high, scale_low, PML_thickness,...
            d_SRC2_rho, d_SRC1_rho, d_SRC2_z, d_SRC1_z, height_Rx,...
            d_RX2_rho, d_RX1_rho, d_RX1_z, d_RX2_z, do_annimation,...
            sigma_air, height_Rx_opt, number_of_Rx)
        
        
    fontsize = 15;
    
    % Plot results
    fig221 = figure(221);
    set(fig221, 'Name', 'Output - Animation', 'NumberTitle','off');  
    cla
    hold on
  
    % plot entire domain on left
    x = 0.1; y = 0.1; w = 0.4; h = 0.8;
    subplot('Position', [x y w h])
    surf(RHO, Z, Hz_tot', 'EdgeColor', 'None', 'facecolor',...
        'interp', 'FaceAlpha', 0.7)
    
    line([d_SRC2_rho+drho d_SRC1_rho],[d_SRC2_z d_SRC1_z], [100 100], 'LineWidth', 5)
    line([d_SRC2_rho+drho d_SRC1_rho],[height_Rx_opt height_Rx_opt], [100 100], 'LineWidth', 3, 'LineStyle', ':')
    
    for i = 1:number_of_Rx
        
        line([d_RX2_rho(i) d_RX1_rho(i)],[d_RX2_z(i) d_RX1_z(i)], [100 100], 'LineWidth', 5)
        
    end

    view(2)
    colormap jet;  
    colorbar;
    h = colorbar;
    ylabel(h, 'H_z [A/m]', 'fontsize', fontsize)

    % set color bar limits for best visuals
    caxis([-scale_high scale_high]);
    
    % insert water line
    if ((sigma_H2O ~= 0) && (sigma_air == 0))
    
        yline(height_H2O);
        
    end

    txt1 = sprintf('H_z: n = %d of %d, \\Deltat = %1.2f ps, t = %1.2f \\mus, f = %1.2f MHz, \\theta_{Rx} = %1.1f^o, \\Delta\\theta_{Rx} = %1.1f^o, \\Delta\\rho = \\Deltaz = %1.2f mm',...
        n, T_STEPS, dt*10^(12), n*dt*10^(6), freq*10^(-6), angle_degrees_start, angle_degrees_diff, drho*1000);
       
    title({txt1})
    ax = gca;
    ax.FontSize = fontsize;
    
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    daspect([1 1 1])
    xlabel('\rho-position [m]', 'fontsize', fontsize)
    ylabel('z-position [m]', 'fontsize', fontsize)
    xlim([drho*((-Nrho/2)+PML_thickness) drho*(((Nrho/2)+1)-PML_thickness)])
    %xticks(drho*(-Nrho/2):1:drho*((Nrho/2)+1))
    ylim([dz*PML_thickness dz*((Nz+1)-PML_thickness)])
    yticks(0:0.25:dz*(Nz+1))
    
    if ((sigma_H2O ~= 0) && (sigma_air == 0))
        
        txt = 'H_2O';
        text(0.8*drho*((-Nrho/2)+PML_thickness), 0.95*height_H2O,txt, 'HorizontalAlignment', 'center', 'fontsize', fontsize)
        txt = 'Air';
        text(0.8*drho*((-Nrho/2)+PML_thickness), 1.05*height_H2O,txt, 'HorizontalAlignment', 'center', 'fontsize', fontsize)
    
    end
    
    txt = 'Tx Coil \rightarrow';
    text(0.25*drho*((-Nrho/2)+PML_thickness), d_SRC1_z, txt, 'HorizontalAlignment', 'right', 'fontsize', fontsize)
    txt = '\leftarrow Optimal Rx Coil';
    text(0.15*drho*(((Nrho/2)+1)-PML_thickness), height_Rx_opt, txt, 'HorizontalAlignment', 'left', 'fontsize', fontsize)
    txt = 'Misaligned Rx Coil';
    text(d_RX2_rho, 1.1*height_Rx, txt, 'HorizontalAlignment', 'left', 'fontsize', fontsize)
    txt = '\downarrow';
    text(d_RX2_rho, 1.07*height_Rx, txt, 'HorizontalAlignment', 'left', 'fontsize', fontsize)
   
    % plot close-up of interface on right top
    x = 0.6; y = 0.5; w = 0.3; h = 0.4;
    subplot('Position', [x y w h])
    surf(RHO, Z, 1000*Hz_tot', 'EdgeColor', 'None', 'facecolor',...
        'interp', 'FaceAlpha', 0.7)

    for i = 1:number_of_Rx
        
        line([d_RX2_rho(i) d_RX1_rho(i)],[d_RX2_z(i) d_RX1_z(i)], [100 100], 'LineWidth', 5)
        
    end
    
    line([d_SRC2_rho+drho d_SRC1_rho],[height_Rx_opt height_Rx_opt], [100 100], 'LineWidth', 3, 'LineStyle', ':')

    view(2)
    colormap jet;  
    colorbar;
    h = colorbar;
    ylabel(h, 'H_z [mA/m]', 'fontsize', fontsize)

    % set color bar limits for best visuals
    caxis([-scale_low scale_low]);

    title('Scale for Field-Strength Adjusted Near Rx Coil.')

    ax = gca;
    ax.FontSize = fontsize;
    
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    daspect([1 1 1])
    ylabel('z-position [m]', 'fontsize', fontsize)
    xlim([drho*((-Nrho/2)+PML_thickness) drho*(((Nrho/2)+1)-PML_thickness)])
    ylim([height_H2O dz*((Nz+1)-PML_thickness)])
    set(gca,'xticklabel',[])
    yticks(0:0.25:dz*(Nz+1))  
    
    if ((sigma_H2O ~= 0) && (sigma_air == 0))
        
        txt = 'Air';
        text(0.8*drho*((-Nrho/2)+PML_thickness), 1.05*height_H2O,txt, 'HorizontalAlignment', 'center', 'fontsize', fontsize)
    
    end
    
    txt = '\leftarrow Optimal Rx Coil';
    text(0.15*drho*(((Nrho/2)+1)-PML_thickness), height_Rx_opt, txt, 'HorizontalAlignment', 'left', 'fontsize', fontsize)
    txt = 'Misaligned Rx Coil';
    text(d_RX2_rho, 1.1*height_Rx, txt, 'HorizontalAlignment', 'left', 'fontsize', fontsize)
    txt = '\downarrow';
    text(d_RX2_rho, 1.07*height_Rx, txt, 'HorizontalAlignment', 'left', 'fontsize', fontsize)
   
	ax = gca;
    ax.FontSize = fontsize;
    
    % plot close-up of interface on right bottom
    x = 0.6; y = 0.1; w = 0.3; h = 0.4;
    subplot('Position', [x y w h])
    surf(RHO, Z, Hz_tot', 'EdgeColor', 'None', 'facecolor',...
        'interp', 'FaceAlpha', 0.7)
    
    line([d_SRC2_rho+drho d_SRC1_rho],[d_SRC2_z d_SRC1_z], [100 100], 'LineWidth', 5)

    view(2)
    colormap jet;  
    colorbar;
    h = colorbar;
    ylabel(h, 'H_z [A/m]', 'fontsize', fontsize)

    % set color bar limits for best visuals
    caxis([-scale_high scale_high]);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    daspect([1 1 1])
    xlabel('\rho-position [m]', 'fontsize', fontsize)
    ylabel('z-position [m]', 'fontsize', fontsize)
    xlim([drho*((-Nrho/2)+PML_thickness) drho*(((Nrho/2)+1)-PML_thickness)])
    ylim([dz*PML_thickness height_H2O])
    yticks(0:0.25:dz*(Nz+1)) 
    
    if ((sigma_H2O ~= 0) && (sigma_air == 0))
        
        txt = 'H_2O';
        text(0.8*drho*((-Nrho/2)+PML_thickness), 0.95*height_H2O,txt,'HorizontalAlignment','center', 'fontsize', fontsize)
    
    end
    
    txt = '\leftarrow Tx Coil';
    text(0.15*drho*(((Nrho/2)+1)-PML_thickness), d_SRC1_z, txt, 'HorizontalAlignment', 'left', 'fontsize', fontsize)
    
    ax = gca;
    ax.FontSize = fontsize;
    
    drawnow
    currFrame = getframe(gcf);
    
    if (do_annimation == true)
        writeVideo(vidObj, currFrame);
    end
    
end