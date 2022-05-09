function a = sommerfeld(x, D, T, Z, mode)

    u = sqrt(x.*x + (T.^2)*1i);
    
    if T == 0 % special case
        a = x.*exp(-u).*exp(-x.*Z)/2;
    else
        a = (x.^2).*exp(-u).*exp(-x.*Z)/(x+u);
    end
    
    if mode == 1
        a = x.*besselj(1, x.*D).*a; % uplink P field (used here)
    elseif mode == 2
        a = u.*besselj(1, x.*D).*a; % downlink P field
    elseif mode == 3
        a = x.*besselj(0, x.*D).*a; % Q field (used here)
    else
        a = 0;
    end

end