% SIMULATION SETUP
M = 64;          % Number of subcarriers
N = 30;          % Number of subsymbols per frame
df = 15e3;       % LTE subcarrier spacing (15 kHz)
fc = 5e9;        % Carrier frequency (5 GHz)
padLen = 10;     % Number of padding samples, should be greater than channel dispersion
padType = 'ZP';  % Zero Padding (ZP) to mitigate Inter-Symbol Interference (ISI)
SNRdB_values = linspace(0, 20, 12);
numb_slots = 100; % Number of slots, each slot contains 30 symbols

channel_type = 'awgn';

% Pilot generation and grid population
pilotBin = floor(N/2)+1;
Pdd = zeros(M,N);
Pdd(1,pilotBin) = exp(1i*pi/4); % populate just one bin to se

% Configure paths
chanParams.pathDelays      = [0 5 8]; % num muestras de retraso 
%al principio 0 estamos sincronizados con base, luego 5 muestras -3, se
%aleja, luego 8 muestra se acerca 5 en doppler
chanParams.pathGains       = [1  0.7 0.5]; % ganancia del camino 
chanParams.pathDopplers    = [0  0  1 ]; % indices doppler(multiplos fsamp/MN)

fsamp = M*df;            % freq muestreo M*df
Meff = M + padLen;       % numero efect muestras por simb M+padlen(padding)
numSamps = Meff * N;     % Num muestras por OTFS sub-simb incluye padding
T = ((M+padLen)/(M*df)); % duracion del simbolo, Meff/numSamps

% Para calcular la freq doppler en hz
chanParams.pathDopplerFreqs = chanParams.pathDopplers * 1/(N*T); % Hz

% Transmit pilots over all subcarriers and symbols to sound the channel
txOut0 = ofdmmod(exp(1i*pi/4)*ones(M,N),M,padLen);       % transmit pilots over the entire grid

% Add white Gaussian noise
Es = mean(abs(pskmod(0:3,4,pi/4).^ 2));

berOFDM_values = zeros(size(SNRdB_values));
blerOFDM_values = zeros(size(SNRdB_values));
num_tx_bits = 2*M*N;

for i_snr=1:length(SNRdB_values)
    SNRdB = SNRdB_values(i_snr);
    n_errors = 0;
    block_error = 0;
    n_blk_errors = 0;
    n0 = Es/(10^(SNRdB/10));
    if ~strcmp(channel_type, 'awgn')
        dopplerOut1 = dopplerChannel(txOut0,fsamp,chanParams);% send through channel
        chOut = awgn(dopplerOut1,SNRdB,'measured');              % add noise
        Yofdm = ofdmdemod(chOut(1:(M+padLen)*N),M,padLen);      % demodulate
        Hofdm = Yofdm * conj(Pdd(1,pilotBin)) / (abs(Pdd(1,pilotBin))^2 + n0); % LMMSE channel estimate
    end
        for i_slot = 1:numb_slots
            % Data generation
            Xgrid = zeros(M,N);
            Xdata = randi([0,1],2*M,N);
            Xgrid(1:M,:) = pskmod(Xdata,4,pi/4,InputType="bit");
        
            % Transmit data over the same channel and use channel estimates to equalize
            txOut = ofdmmod(Xgrid,M,padLen);                        % transmit data grid
            
            if ~strcmp(channel_type, 'awgn')
                dopplerOut2 = dopplerChannel(txOut,fsamp,chanParams);    % send through channel
            else
                dopplerOut2 = txOut;
            end
            chOut = awgn(dopplerOut2,SNRdB,'measured');              % add noise
            
            rxWindow = chOut(1:(M+padLen)*N);
            
            if ~strcmp(channel_type, 'awgn')
                Yofdm = ofdmdemod(rxWindow,M,padLen);                   % demodulate
                Xhat_ofdm = conj(Hofdm) .* Yofdm ./ (abs(Hofdm).^2+n0); % equalize with LMMSE
            else
                Xhat_ofdm= ofdmdemod(rxWindow,M,padLen);
            end
            XhatDataOFDM = pskdemod(Xhat_ofdm,4,pi/4, ...
                OutputType="bit",OutputDataType="logical");         % decode
            bit_errors = xor(Xdata,XhatDataOFDM);
            errores_ahora = sum(bit_errors(:));
            % Cuentas errores y acumulas
            n_errors = errores_ahora + n_errors;
    
            if(n_errors>0)
                block_error = 1 ;% Si es erróneo
            else
                block_error = 0;
            end
    
            n_blk_errors = n_blk_errors + block_error;
        
        end 
    % Imprimir la BER después de completar todos los slots de una SNR
    %fprintf('SNR = %2d dB -> BER = %.6f | BLER = %.6f | Errores: %d de %d bits\n', ...
     %   SNRdB, berOTFS_values(i_snr), blerOTFS_values(i_snr), errores_ahora, num_tx_bits);
    berOFDM_values(i_snr) = n_errors/(num_tx_bits*numb_slots);
    blerOFDM_values(i_snr) = n_blk_errors/numb_slots;

end

% Reemplazar valores 0 por un número pequeño para evitar problemas con log
berOFDM_values(berOFDM_values == 0) = 1e-10;
blerOFDM_values(blerOFDM_values == 0) = 1e-10;

% Graficar la BER contra los valores de SNR
figure;
semilogy(SNRdB_values, berOFDM_values, '-o', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs. SNR para OFDM');
grid on;
ylim([1e-4, 1]); % Ajusta los límites del eje Y (ejemplo)

% Graficar la BLER contra los valores de SNR
figure;
semilogy(SNRdB_values, blerOFDM_values, '-o', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('BLER');
title('BLER vs. SNR para OFDM');
grid on;
ylim([1e-1, 1]); % Ajusta los límites del eje Y (ejemplo)

