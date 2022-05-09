function [Eps, mySigma]...
            = MI_materialGrid(Nrho, Nz, Eps0, j_H2O_height, sigma_H2O,...
                sigma_air, Er_H2O, Er_air, drho, dz, RHO, Z)

    %material  
    E_H2O = Er_H2O*Eps0;      
    E_air = Er_air*Eps0;
    
    mySigma = zeros(Nrho, Nz);
    Eps = zeros(Nrho, Nz);
    
    % define permitivity
    for i = 1:Nrho
        
        for j = 1:Nz
            
            if (j == j_H2O_height)
            
                Eps(i,j) = (E_H2O + E_air)/2;
                
            elseif (j < j_H2O_height)
                
                Eps(i,j) = E_H2O;
                
            else
                
                Eps(i,j) = E_air;
                
            end
            
        end
        
    end
    
    fig3 = figure(3);
    set(fig3, 'Name', 'Debug - Permitivities', 'NumberTitle','off');
    movegui('west');
    axis equal
    axis tight
    cla
    hold on
    surf(RHO, Z, Eps'/Eps0, 'FaceAlpha', 0.7);
    view(2);
    colormap(flipud(winter));
    colorbar
    h = colorbar;
    txt = sprintf('\\epsilon_r [rad]');
    ylabel(h, txt)
    txt = sprintf('Relative Permitivity matrix');
    title(txt)
    xlabel('\rho-position [m]')
    ylabel('z-position [m]')
    drawnow
    
    % define conductivity
    for i = 1:Nrho
        
        for j = 1:Nz
            
            if (j == j_H2O_height) 
                
                mySigma(i,j) = (sigma_H2O + sigma_air)/2;
                
            elseif (j < j_H2O_height)
                
                mySigma(i,j) = sigma_H2O;
                
            else
                
                mySigma(i,j) = sigma_air; % in air
                
            end
            
        end
        
    end    
    
    fig4 = figure(4);
    set(fig4, 'Name', 'Debug - Conductivities', 'NumberTitle','off');
    movegui('west');
    axis equal
    axis tight
    cla
    hold on
    surf(RHO, Z, mySigma', 'FaceAlpha', 0.7);
    view(2);
    colormap(flipud(parula));
    colorbar
    h = colorbar;
    txt = sprintf('\\sigma [S/m]');
    ylabel(h, txt)
    txt = sprintf('Conductivity matrix');
    title(txt)
    xlabel('\rho-position [m]')
    ylabel('z-position [m]')
    drawnow
    
end