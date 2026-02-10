function Es = computeSymbolEnergy(mod_order)
%COMPUTESYMBOOLENERGY Computes average symbol energy for a given modulation order.
%   Es = computeSymbolEnergy(mod_order) returns the average symbol energy Es 
%   for the modulation scheme defined by mod_order.
%   Supports BPSK (2), QPSK (4), 16-QAM (16), and 64-QAM (64).

    switch mod_order
        case 2  % BPSK
            symbols = pskmod(0:mod_order-1, mod_order, pi);
        case 4  % QPSK
            symbols = pskmod(0:mod_order-1, mod_order, pi/4);
        case {16, 64}  % QAM
            symbols = qammod(0:mod_order-1, mod_order);
        otherwise
            error('Unsupported modulation order: %d', mod_order);
    end

    % Compute average symbol energy
    Es = mean(abs(symbols).^2);
end
