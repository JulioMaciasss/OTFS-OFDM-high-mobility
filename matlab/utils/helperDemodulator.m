function bits = helperDemodulator(simb_mod, mod_order)
    switch mod_order
        case 2
            bits = pskdemod(simb_mod, mod_order, pi, 'OutputType', 'bit', 'OutputDataType', 'logical');
        case 4
            bits = pskdemod(simb_mod, mod_order, pi/4, 'OutputType', 'bit', 'OutputDataType', 'logical');
        case {16, 64}
            bits = qamdemod(simb_mod, mod_order, 'OutputType', 'bit', 'OutputDataType', 'logical');
        otherwise
            error('Modulation not supported');
    end
end