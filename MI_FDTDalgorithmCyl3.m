function [Ephi_tot, Hz_tot, Hrho_tot, IEphi, ICErho, ICEz, IHrho, IHz, ICHphi]...
            = MI_FDTDalgorithmCyl3(Nrho, Nz, drho, dz, Ephi_tot, Hz_tot, Hrho_tot,...
                A, B, C, D, E, F, G, H, I, J, K, L,...
                IEphi, Curl_H_phi, curl_E_rho, ICErho, IHrho, curl_E_z, ICEz, IHz, ICHphi)
            
    % Main FDTD loop
    % solves 2-D scattering problem with PML and Dirichlet BC
    % in cylindrical coordinates
    for j = 2:Nz
        
        for i = (Nrho/2)+1:Nrho
            
            % Integrating Ephi, at each point, over time
            IEphi(i,j) = IEphi(i,j) + Ephi_tot(i,j);
            % phi component of the curl of H
            Curl_H_phi(i,j) = ((Hz_tot(i,j) - Hz_tot(i-1,j))/drho - (Hrho_tot(i,j) - Hrho_tot(i,j-1))/dz);
            % Integrating curl_H_phi, at each point, over time
            ICHphi(i,j) = ICHphi(i,j) + Curl_H_phi(i,j);
            %update equation for Ephi
            Ephi_tot(i,j) = I(i,j)*Ephi_tot(i,j) + J(i,j)*Curl_H_phi(i,j) + K(i,j)*ICHphi(i,j) + L(i,j)*IEphi(i,j); 
            
        end
        
    end
    
    for j = 1:Nz-1
        
        % Update Hrho from E
        for i = (Nrho/2)+1:Nrho
            
            % rho component of the curl of E
            curl_E_rho(i,j) = (Ephi_tot(i,j+1) - Ephi_tot(i,j))/dz;
            % Integrating curl_E_rho, at each point, over time
            ICErho(i,j) = ICErho(i,j) + curl_E_rho(i,j);
            % Integrating Hrho, at each point, over time
            IHrho(i,j) = IHrho(i,j) + Hrho_tot(i,j);
            % update equation for Hrho
            Hrho_tot(i,j) = A(i,j)*Hrho_tot(i,j) + B(i,j)*curl_E_rho(i,j) + C(i,j)*ICErho(i,j) + D(i,j)*IHrho(i,j);

        end
        
    end
    
    % BC at j = Nz
    % rho component of the curl of E
    curl_E_rho(:,Nz) = (0 - Ephi_tot(:,Nz))/dz;
    % Integrating curl_E_rho, at each point, over time
    ICErho(:,Nz) = ICErho(:,Nz) + curl_E_rho(:,Nz);
    % Integrating Hrho, at each point, over time
    IHrho(:,Nz) = IHrho(:,Nz) + Hrho_tot(:,Nz);
    % update equation for Hrho
    Hrho_tot(:,Nz) = A(:,Nz).*Hrho_tot(:,Nz) + B(:,Nz).*curl_E_rho(:,Nz) + C(:,Nz).*ICErho(:,Nz) + D(:,Nz).*IHrho(:,Nz);
    
    for j = 2:Nz
        
        % Update Hz from E
        for i = (Nrho/2)+1:Nrho-1

            % z component of the curl of E
            curl_E_z(i,j) = (1/(((i-((Nrho/2)+1))*drho)+drho/2))*Ephi_tot(i,j) + (Ephi_tot(i+1,j) - Ephi_tot(i,j))/drho; % nonconservitive form
            % Integrating curl_E_z, at each point, over time
            ICEz(i,j) = ICEz(i,j) + curl_E_z(i,j);
            % Integrating Hz, at each point, over time
            IHz(i,j) = IHz(i,j) + Hz_tot(i,j);
            % update equation for Hz
            Hz_tot(i,j) = E(i,j)*Hz_tot(i,j) + F(i,j)*curl_E_z(i,j) + G(i,j)*ICEz(i,j) + H(i,j)*IHz(i,j);    
            
        end
        
    end
    
    % BC at i = Nrho
    % z component of the curl of E
    curl_E_z(Nrho,:) = (1/(Nrho*drho))*Ephi_tot(Nrho,:) + (0 - Ephi_tot(Nrho,:))/drho; % nonconservitive form
    % Integrating curl_E_z, at each point, over time
    ICEz(Nrho,:) = ICEz(Nrho,:) + curl_E_z(Nrho,:);
    % Integrating Hz, at each point, over time
    IHz(Nrho,:) = IHz(Nrho,:) + Hz_tot(Nrho,:);
    % update equation for Hy
    Hz_tot(Nrho,:) = E(Nrho,:).*Hz_tot(Nrho,:) + F(Nrho,:).*curl_E_z(Nrho,:) + G(Nrho,:).*ICEz(Nrho,:) + H(Nrho,:).*IHz(Nrho,:);
    
    
    % apply cylindrical symmetry
    Hrho_tot(Nrho/2:-1:1, :) = Hrho_tot((Nrho/2)+1:Nrho, :);
    Hz_tot(Nrho/2:-1:1, :) = Hz_tot((Nrho/2)+1:Nrho, :);
    Ephi_tot(Nrho/2:-1:1, :) = Ephi_tot((Nrho/2)+1:Nrho, :);
       
end