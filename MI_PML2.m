function [A, B, C, D, E, F, G, H, I, J, K, L]...
            = MI_PML2(Nrho, Nz, dt, Eps, u0, mySigma, PML_thickness, do_PML)             
    %%
    %initialize conductivities of PMLs
    sigma_rho = zeros(Nrho, Nz);
    sigma_z = zeros(Nrho, Nz);
    sigma_phi = zeros(Nrho, Nz);
    
    % for Hrho
    mHrho0 = zeros(Nrho, Nz);
    mHrho1 = zeros(Nrho, Nz);
    mHrho2 = zeros(Nrho, Nz);
    mHrho3 = zeros(Nrho, Nz);
    mHrho4 = zeros(Nrho, Nz);
    A = zeros(Nrho, Nz);
    B = zeros(Nrho, Nz);
    C = zeros(Nrho, Nz); 
    D = zeros(Nrho, Nz); 
    
    % for Hz
    mHz0 = zeros(Nrho, Nz);
    mHz1 = zeros(Nrho, Nz);
    mHz2 = zeros(Nrho, Nz);
    mHz3 = zeros(Nrho, Nz);
    mHz4 = zeros(Nrho, Nz);
    E = zeros(Nrho, Nz);
    F = zeros(Nrho, Nz);
    G = zeros(Nrho, Nz);  
    H = zeros(Nrho, Nz);  
    
    % for Ephi
    mEphi0 = zeros(Nrho, Nz);
    mEphi1 = zeros(Nrho, Nz);    
    mEphi2 = zeros(Nrho, Nz);
    mEphi3 = zeros(Nrho, Nz);
    mEphi4 = zeros(Nrho, Nz);
    I = zeros(Nrho, Nz);    
    J = zeros(Nrho, Nz);
    K = zeros(Nrho, Nz);
    L = zeros(Nrho, Nz);
    
    %%
    % build PML
    if (do_PML == true)
        
        for (j = 1:1:Nz) 

            for (i = 1:1:PML_thickness)

                sigma_rho(PML_thickness-i+1,j) = (Eps(PML_thickness-i+1,j)/(2*dt))*(i/PML_thickness)^3; 
                sigma_rho(Nrho-PML_thickness+i,j) = (Eps(Nrho-PML_thickness+i,j)/(2*dt))*(i/PML_thickness)^3;
                
                sigma_phi(PML_thickness-i+1,j) = (Eps(PML_thickness-i+1,j)/(2*dt))*(1/PML_thickness)^3;
                sigma_phi(Nrho-PML_thickness+i,j) = (Eps(Nrho-PML_thickness+i,j)/(2*dt))*(1/PML_thickness)^3;
                
            end

        end

        for (i = 1:1:Nrho)

            for (j = 1:1:PML_thickness)

                sigma_z(i,PML_thickness-j+1) = (Eps(i,PML_thickness-j+1)/(2*dt))*(j/PML_thickness)^3; 
                sigma_z(i,Nz-PML_thickness+j) = (Eps(i,Nz-PML_thickness+j)/(2*dt))*(j/PML_thickness)^3; 
                
                sigma_phi(i,PML_thickness-j+1) = (Eps(i,PML_thickness-j+1)/(2*dt))*(1/PML_thickness)^3; 
                sigma_phi(i,Nz-PML_thickness+j) = (Eps(i,Nz-PML_thickness+j)/(2*dt))*(1/PML_thickness)^3;
                
            end

        end
    
    else % (do_PML == false) no PML
        
        % do nothing (sigmas all zero)      
        
    end
    
    %%
    % generate constants for FDTD algorithm
    %for Ezphi
    for (j = 1:Nz)
        
        for (i = 1:Nrho)
            
           mEphi0(i,j) = (Eps(i,j)/dt)+((sigma_rho(i,j)+sigma_z(i,j)+mySigma(i,j))/2)+((sigma_rho(i,j)*sigma_z(i,j)*dt)/(4*Eps(i,j)));
           mEphi1(i,j) = (Eps(i,j)/dt)-((sigma_rho(i,j)+sigma_z(i,j)+mySigma(i,j))/2)-((sigma_rho(i,j)*sigma_z(i,j)*dt)/(4*Eps(i,j)));
           mEphi2(i,j) = 1;
           mEphi3(i,j) = -dt*sigma_phi(i,j)/Eps(i,j);
           mEphi4(i,j) = -(dt*sigma_rho(i,j)*sigma_z(i,j)/Eps(i,j));
           I(i,j) = mEphi1(i,j)/mEphi0(i,j); % for previous time step
           J(i,j) = mEphi2(i,j)/mEphi0(i,j); % for curl term
           K(i,j) = mEphi3(i,j)/mEphi0(i,j); % for integral term involving curl
           L(i,j) = mEphi4(i,j)/mEphi0(i,j); % for integral term involving E
           
        end
        
    end     
    
    %for Hrho
    for (j = 1:Nz)
        
        for (i = 1:Nrho)
            
           mHrho0(i,j) = (1/dt) + ((sigma_z(i,j)+sigma_phi(i,j))/(2*Eps(i,j))) + (dt*sigma_z(i,j)*sigma_phi(i,j))/(4*Eps(i,j)*Eps(i,j));
           mHrho1(i,j) = (1/dt) - ((sigma_z(i,j)+sigma_phi(i,j))/(2*Eps(i,j))) - (dt*sigma_z(i,j)*sigma_phi(i,j))/(4*Eps(i,j)*Eps(i,j));
           mHrho2(i,j) = -(1/u0);
           mHrho3(i,j) = -(dt*sigma_rho(i,j))/(u0*Eps(i,j));
           mHrho4(i,j) = -(dt*sigma_z(i,j)*sigma_phi(i,j))/(Eps(i,j)*Eps(i,j));
           A(i,j) = mHrho1(i,j)/mHrho0(i,j); % for previous time step
           B(i,j) = mHrho2(i,j)/mHrho0(i,j); % for curl term
           C(i,j) = mHrho3(i,j)/mHrho0(i,j); % for integral term involving curl
           D(i,j) = mHrho4(i,j)/mHrho0(i,j); % for integral term involving H
           
        end
        
    end 
    
    %for Hz
    for (j = 1:Nz)
        
        for (i = 1:Nrho)
            
           mHz0(i,j) = (1/dt) + ((sigma_rho(i,j)+sigma_phi(i,j))/(2*Eps(i,j))) + (dt*sigma_rho(i,j)*sigma_phi(i,j))/(4*Eps(i,j)*Eps(i,j));
           mHz1(i,j) = (1/dt) - ((sigma_rho(i,j)+sigma_phi(i,j))/(2*Eps(i,j))) - (dt*sigma_rho(i,j)*sigma_phi(i,j))/(4*Eps(i,j)*Eps(i,j));
           mHz2(i,j) = (1/u0);
           mHz3(i,j) = (dt*sigma_z(i,j))/(u0*Eps(i,j));
           mHz4(i,j) = (dt*sigma_rho(i,j)*sigma_phi(i,j))/(Eps(i,j)*Eps(i,j));
           E(i,j) = mHz1(i,j)/mHz0(i,j); % for previous time step
           F(i,j) = mHz2(i,j)/mHz0(i,j); % for curl term
           G(i,j) = mHz3(i,j)/mHz0(i,j); % for integral term involving curl
           H(i,j) = mHz4(i,j)/mHz0(i,j); % for integral term involving H
           
        end
        
    end 
    
    % for debugging PML
    fig5 = figure(5);
    set(fig5, 'Name', 'Debug - Perfectly Matched Layer', 'NumberTitle','off');
    movegui('west');
    axis equal
    axis tight
    cla
    hold on
    surf(sigma_rho' + sigma_z', 'EdgeColor', 'None', 'FaceAlpha', 0.7); 
    view(2);
    xlim([0 Nrho+1])
    ylim([0 Nz+1])
    colormap cool;
    caxis([0 max(sigma_rho(1,:))]);
    colorbar;
    h = colorbar;
    txt = sprintf('\\sigma [S/m]');
    ylabel(h, txt)
    xlabel('i')
    ylabel('j')
    title('The Perfectly Matched Layer')
    drawnow
    
    % free up memory
    clear mHz0 mHz1 mHz2 mHz3 mHz4
    clear mHrho0 mHrho1 mHrho2 mHrho3 mHrho4
    clear mEphi0 mEphi1 mEphi2 mEphi3 mEphi4
    
end