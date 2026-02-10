function G = getG(M,N,chanParams,padLen,padType)
    % Form time domain channel matrix from detected DD paths
    if strcmp(padType,'ZP') || strcmp(padType,'CP')
        Meff = M + padLen;  % account for subsymbol pad length in forming channel
        lmax = padLen;      % max delay
    else
        Meff = M;
        lmax = max(chanParams.pathDelays);  % max delay
    end
    MN = Meff*N;
    P = length(chanParams.pathDelays);  % number of paths
    
    % Form an array of channel responses for each path
    g = zeros(lmax+1,MN);
    for p = 1:P
        gp = chanParams.pathGains(p);
        lp = chanParams.pathDelays(p);
        vp = chanParams.pathDopplers(p); 

        % For each DD path, compute the channel response.
        % Each path is a complex sinusoid at the Doppler frequency (kp)
        % shifted by a delay (lp) and scaled by the path gain (gp)
        g(lp+1,:) = g(lp+1,:) + gp*exp(1i*2*pi/MN * vp*((0:MN-1)-lp));
    end    

    % Form the MN-by-MN channel matrix G
    G = zeros(MN,MN);
    % Each DD path is a diagonal in G offset by its path delay l
    for l = unique(chanParams.pathDelays).'
        G = G + diag(g(l+1,l+1:end),-l);
    end
end