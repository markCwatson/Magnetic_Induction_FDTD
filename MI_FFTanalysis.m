function MI_FFTanalysis(sample_Tx, sample_Rx, dt, figNum, type)
    
    % for windowing set window_start != 1
    window_start = 1;
    sample = [sample_Tx(window_start:end); sample_Rx(window_start:end)];    
    Fs = 1/dt;
    T_STEPS = length(sample_Tx(window_start:end));
    f = Fs*(0:(T_STEPS/2))/T_STEPS;
    % note: The frequency resolution is equal to the sampling...
    %... frequency divided by FFT size
        
    for (i = 1:1:2)
        
        % Fourier analysis on signal
        Y = fft(sample(i,:));
        P2_mag = abs(Y/T_STEPS); % magnitude of Y in [A/m] or [V/m]
        P1_mag = P2_mag(1:floor(T_STEPS/2+1));
        P1_mag(2:end-1) = 2*P1_mag(2:end-1);
        
        P2_phase = 180*angle(Y)/pi; % phase of Y in degrees
        P1_phase = P2_phase(1:floor(T_STEPS/2+1));
        
        fig = figure(figNum);
        set(fig, 'Name', 'FFT Analysis', 'NumberTitle','off');
        %subplot(2,1,1) % uncomment this and below if you want phase
        hold on
        
        if (i == 1) % Tx
            
            semilogx(f, 20*log10(P1_mag), 'b-o') % in dB
            
        elseif (i == 2) % Rx
            
            semilogx(f, 20*log10(P1_mag), 'r-o') % in dB
        
        end
        
        %xlim([0 0.1*max(f)*10^(-6)])
        xlabel('f (MHz)')
        
        if (i == 1)
            
            if (type == 'H')
            
                ylabel('20*log_1_0(|H_y(f)|) [dB]')
                title('FFT of H_y at Tx and Rx Locations')
            
            elseif (type == 'E')
            
                ylabel('20*log_1_0(|E_z(f)|) [dB]')
                title('FFT of E_z at Tx and Rx Locations')
            
            elseif (type == 'V')
            
                ylabel('20*log_1_0(|V(f)|) [dBV]')
                title('FFT of Voltages in Tx and Rx')
            
            elseif (type == 'I')
            
                ylabel('20*log_1_0(|I(f)|) [dBA]')
                title('FFT of Currents in Tx and Rx')
            
            end
            
        elseif (i == 2)
            
           legend('Tx','Rx', 'Location', 'NorthEast') 
           
        end    
        
        % for phase
%         fig = figure(figNum);
%         set(fig, 'Name', 'FFT Analysis', 'NumberTitle','off');
%         subplot(2,1,2)
%         hold on
%         
%         plot(f.*10^(-6),P1_phase, '-o')
%         %xlim([0 2*freq*10^(-6)])
%         
%         % for zoomed in mag
% %         plot(f.*10^(-6),20*log10(P1_mag), '-o') % in dB
% %         xlim([0 1000]) %in MHz       
%         
%         if (i == 1)
%             
%             if (type == 'H')
%                 
%                 ylabel('phase(H_y(f)) [deg]')
%                 
%             elseif (type == 'E')
%                 
%                 ylabel('phase(E_z(f)) [deg]')
%                 
%             elseif (type == 'V')
%             
%                 ylabel('phase(V(f)) [deg]')
%             
%             elseif (type == 'I')
%             
%                 ylabel('phase(I(f)) [deg]')
%             
%             end
%             
%         elseif (i == 2)
%             
%             legend('Tx','Rx', 'Location', 'East')
%             
%         end
%         
%         xlabel('f (kHz)')
        
    end

end