%SIMULATION SET UP
M = 64;          % num subportadoras
N = 30;          % number of subsymb por trama
df = 15e3;       % espaciado de freq LTE
fc = 5e9;        % freq portadora
padLen = 10;     % num muestras relleno, mayor q dispersion del canal
padType = 'ZP';  % zero padding pa mitigar ISI
SNRdB_values = linspace(0, 20, 12);
numb_slots = 100; %cada slot 30 simb

channel_type = '';

% Pilot generation and grid population
pilotBin = floor(N/2)+1;
Pdd = zeros(M,N);
Pdd(1,pilotBin) = exp(1i*pi/4); % populate just one bin to se

% OTFS modulation
txOut0 = helperOTFSmod(Pdd,padLen,padType);

chanParams.pathDelays      = [0  5   8  ]; % number of samples that path is delayed
chanParams.pathGains       = [1  0.7 0.5]; % complex path gain
chanParams.pathDopplers    = [0 -3 5]; % Doppler index as a multiple of fsamp/MN


fsamp = M*df;            % freq muestreo M*df
Meff = M + padLen;       % numero efect muestras por simb M+padlen(padding)
numSamps = Meff * N;     % Num muestras por OTFS sub-simb incluye padding
T = ((M+padLen)/(M*df)); % duracion del simbolo, Meff/numSamps

% Para calcular la freq doppler en hz
chanParams.pathDopplerFreqs = chanParams.pathDopplers * 1/(N*T); % Hz

% Add white Gaussian noise
Es = mean(abs(pskmod(0:3,4,pi/4).^ 2));

berOTFS_values = zeros(size(SNRdB_values));
blerOTFS_values = zeros(size(SNRdB_values));
num_tx_bits = 2*M*N;

for i_snr=1:length(SNRdB_values)
    SNRdB=SNRdB_values(i_snr);
    n_errors = 0;
    block_error = 0;
    n_blk_errors = 0;
    n0 = Es/(10^(SNRdB/10));
    
    %HAGO LA ESTIMACION DEL CANAL
        if ~strcmp(channel_type, 'awgn')

            dopplerOut1 = dopplerChannel(txOut0, fsamp, chanParams);
            chOut = awgn(dopplerOut1,SNRdB,'measured');
            % Get a sample window
            rxIn = chOut(1:numSamps);
            %txOut es el canal, fsamp la freq de muestreo y chanParams los parametros
            % OTFS demodulation
            Ydd = helperOTFSdemod(rxIn,M,padLen,0,padType);
            
            % LMMSE channel estimate in the delay-Doppler 
            Hdd = Ydd * conj(Pdd(1,pilotBin)) / (abs(Pdd(1,pilotBin))^2 + n0);
            
            %threshold = max(0.2, 0.1 * max(abs(Hdd(:)))); % Evita detectar ruido débil
            [lp, vp] = find(abs(Hdd) >= 0.05);
        
            chanEst.pathGains = diag(Hdd(lp,vp));   % get path gains
            chanEst.pathDelays = lp - 1;            % get delay indices
            chanEst.pathDopplers = vp - pilotBin;   % get Doppler indices
        end

    for i_slot = 1:numb_slots
        % Data generation -> GENERO BITS PARA TRANSMITIR
        Xgrid = zeros(M,N);
        Xdata = randi([0,1],2*M,N);
        Xgrid(1:M,:) = pskmod(Xdata,4,pi/4,InputType="bit");
        
        % OTFS modulation -> SEÑAL A TRANSMITIR 
        txOut = helperOTFSmod(Xgrid,padLen,padType);
        
        % Add channel and noise -> METO CANAL Y AÑADO RUIDO
        if ~strcmp(channel_type, 'awgn')
            dopplerOut2 = dopplerChannel(txOut,fsamp,chanParams);
        else
            dopplerOut2 = txOut;
        end

        chOut = awgn(dopplerOut2,SNRdB,'measured');%dopplerout es la señal a la q añado
            %ruido con la realcion señal ruido, y measured indica que debe calcular
            %automente la potencia de la señal antes q el ruid
        
        if ~strcmp(channel_type, 'awgn')
            % Form G matrix using channel estimates
            G = getG(M,N,chanEst,padLen,padType);
        
            rxWindow = chOut(1:numSamps);
            
            %ESTIMO CON LMMSE 
            y_otfs = ((G'*G)+n0*eye(Meff*N)) \ (G'*rxWindow); % LMMSE
        else
            y_otfs = chOut(1:numSamps);
        end

        %DEMODULO Y DETECTO BITS 
        Xhat_otfs = helperOTFSdemod(y_otfs,M,padLen,0,padType); % OTFS demodulation
        XhatDataOTFS = pskdemod(Xhat_otfs,4,pi/4,OutputType="bit",OutputDataType="logical");
        bit_errors = xor(Xdata,XhatDataOTFS);
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
    berOTFS_values(i_snr) = n_errors/(num_tx_bits*numb_slots);
    blerOTFS_values(i_snr) = n_blk_errors/numb_slots;
       
end
% Reemplazar valores 0 por un número pequeño para evitar problemas con log
berOTFS_values(berOTFS_values == 0) = 1e-10;
blerOTFS_values(blerOTFS_values == 0) = 1e-10;

% Graficar la BER contra los valores de SNR
figure;
semilogy(SNRdB_values, berOTFS_values, '-o', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs. SNR para OTFS');
grid on;
ylim([1e-4, 1]); % Ajusta los límites del eje Y (ejemplo)

% Graficar la BLER contra los valores de SNR
figure;
semilogy(SNRdB_values, blerOTFS_values, '-o', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('BLER');
title('BLER vs. SNR para OTFS');
grid on;
ylim([1e-1, 1]); % Ajusta los límites del eje Y (ejemplo)
