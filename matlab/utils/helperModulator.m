function simb_mod = helperModulator(data, mod_order)
    switch mod_order
        case 2
            simb_mod = pskmod(data, mod_order, pi, 'InputType', 'bit');
        case 4
            simb_mod = pskmod(data, mod_order, pi/4, 'InputType', 'bit');
        case {16, 64}
            simb_mod = qammod(data, mod_order, 'InputType', 'bit');
        otherwise
            error('Modulación no soportada');
    end
end