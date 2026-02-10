%SIMULATION SET UP
M = 64;          % num subportadoras
N = 30;          % number of subsymb por trama
df = 15e3;       % espaciado de freq LTE
fc = 3e9;        % freq portadora
padLen = 30;     % num muestras relleno, mayor q dispersion del canal
padType = 'ZP';  % zero padding pa mitigar ISI
SNRdB_values = linspace(0, 20, 12);
numb_slots = 100; %cada slot 30 simb
idealCE =  false;
mod_order = 4; %2->BPSK, 4-> QPSK, 16-> 16-QAM, 64-> 64-QAM
images_turn = false;
channel_response = true;

% create the name depend on the modulation
switch mod_order
    case 2
        modName = 'BPSK';
    case 4
        modName = 'QPSK';
    case 16
        modName = '16QAM';
    case 64
        modName = '64QAM';
    otherwise
        error('Unsupported modulation order: %d', mod_order);
end

% define the basepath
basePath = fullfile('C:\Users\jm-y2\Desktop\TFG\SIMULACIONES', ['OTFS_3GHz_cdlA/', modName]);

% Pilot generation and grid population
pilotBin = floor(N/2)+1;
Pdd = zeros(M,N);
Pdd(1,pilotBin) = exp(1i*pi/4); % populate just one bin to se

% OTFS modulation
txOut0 = helperOTFSmod(Pdd,padLen,padType);

fsamp = M*df;            % freq muestreo M*df
Meff = M + padLen;       % numero efect muestras por simb M+padlen(padding)
numSamps = Meff * N;     % Num muestras por OTFS sub-simb incluye padding
T = ((M+padLen)/(M*df)); % duracion del simbolo, Meff/numSamps

%% Modelo TDL
channel = nrCDLChannel;
% Configuración común para el canal
channel.SampleRate = fsamp;

%Params to change
channel.DelayProfile = 'CDL-A'; 
channel.MaximumDopplerShift = 3;
channel.DelaySpread = 30e-9; 
cdl.TransmitAntennaArray.Size = [1 1 1 1 1];
cdl.ReceiveAntennaArray.Size = [1 1 1 1 1];

threshold = 0.0305;

% Obtener valores numéricos para construir nombre
doppler = channel.MaximumDopplerShift;
delay_ns = round(channel.DelaySpread * 1e9); % convertir a nanosegundos

% Crear nombre de carpeta con formato limpio
dop = sprintf('OTFS_tdlA_Doppler%dHz_Delay_hdd', doppler, delay_ns);

if idealCE
    folderName = strcat(dop, '_idealCE');  % Para Windows

else 
    folderName =dop;  % Para Windows
end

% Ruta completa (usa barra correcta según el SO)
fullFolderPath = fullfile(basePath, folderName);

% Crear la carpeta si no existe
if ~exist(fullFolderPath, 'dir')  
    mkdir(fullFolderPath);
end

Es = computeSymbolEnergy(mod_order);

num_tx_bits = log2(mod_order)*M*N;

berOTFS_values = zeros(size(SNRdB_values));
blerOTFS_values = zeros(size(SNRdB_values));
SINR_W = 0;

for i_snr=1:length(SNRdB_values)
    SNRdB=SNRdB_values(i_snr);
    n_errors = 0;
    block_error = 0;
    n_blk_errors = 0;
    SINR_Slot = 0; %reset SNR for a slot
    n0 = Es/(10^(SNRdB/10));

    for i_slot = 1:numb_slots
        release(channel);
        channel.Seed = i_slot;
        reset(channel);
        
        % Aplico canal 
        dopplerOut1 = channel(txOut0);  % Filtrar la señal transmitida a través del canal
        if idealCE
            chOut = dopplerOut1;  % Si se usa canal perfecto, no se añade ruido
        else
            % Si no se usa canal perfecto, añadir ruido con la relación SNR especificada
            chOut = awgn(dopplerOut1, SNRdB, 'measured');
        end
        
        % Obtener una ventana de muestra de la señal recibida
        rxIn = chOut(1:numSamps, 1);  % Extraer las muestras relevantes de la señal recibida
        
        % Demodulación OTFS de la señal recibida
        Ydd = helperOTFSdemod(rxIn, M, padLen, 0, padType);  % Demodular utilizando OTFS
        
        % Estimación del canal LMMSE en el dominio Delay-Doppler
        if idealCE
            % Si se utiliza estimación ideal de canal, calcular la estimación usando los pilotos
            Hdd = Ydd * conj(Pdd(1, pilotBin)) / (abs(Pdd(1, pilotBin))^2);
        else
            % Si no se utiliza estimación ideal de canal, añadir ruido (n0)
            Hdd = Ydd * conj(Pdd(1, pilotBin)) / (abs(Pdd(1, pilotBin))^2 + n0);
        end
        
        
        %threshold =0.1*max(abs(Hdd(:)));
        %threshold = 0.1 * max(abs(Hdd(:))); % Evita detectar ruido débil
        
        % Filtrar los valores de Hdd que sean suficientemente grandes (umbral 0.05)
        [lp, vp] = find(abs(Hdd) >= threshold);  % Encontrar índices donde la magnitud es suficientemente grande
        
        % Guardar la estimación del canal en la estructura 'chanEst'
        chanEst.pathGains = diag(Hdd(lp, vp));   % Obtener las ganancias de los diferentes caminos
        chanEst.pathDelays = lp - 1;             % Obtener los índices de retardo (delay) de los caminos
        chanEst.pathDopplers = vp - pilotBin;    % Obtener los índices de Doppler de los caminos
        
        reset(channel); 

        % Data generation -> GENERO BITS PARA TRANSMITIR
        Xgrid = zeros(M,N);
        Xdata = randi([0,1],log2(mod_order)*M,N);
        Xgrid(1:M,:) = helperModulator(Xdata, mod_order);

        % OTFS modulation -> SEÑAL A TRANSMITIR 
        txOut = helperOTFSmod(Xgrid,padLen,padType);
        
        % Add channel and noise -> METO CANAL Y AÑADO RUIDO
        dopplerOut2 = channel(txOut);
        chOut = awgn(dopplerOut2,SNRdB,'measured');
        
        % Form G matrix using channel estimates
        G = getG(M,N,chanEst,padLen,padType);
    
        rxWindow = chOut(1:numSamps,1);
        
        %ESTIMO CON LMMSE 
        y_otfs = ((G'*G)+n0*eye(Meff*N)) \ (G'*rxWindow); % LMMSE
        
        %DEMODULO Y DETECTO BITS 
        Xhat_otfs = helperOTFSdemod(y_otfs,M,padLen,0,padType); % OTFS demodulation
        %demodulo
        XhatDataOTFS = helperDemodulator(Xhat_otfs, mod_order);

        SINR_now_W = 1/mean(abs(Xhat_otfs(:)-Xgrid(:)).^2);%SNRI
        SINR_Slot = SINR_Slot + SINR_now_W;

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
    %Print channel response
    if channel_response && i_snr == 1
        % Graficar la Respuesta al Canal (Hdd)
        figure;
        xa = -N/2 : N/2 - 1;  % Eje Doppler centrado
        ya = 0:1:M-1;         % Eje de muestras de tiempo (o delay)
        
        mesh(xa, ya, abs(Hdd));  % Muestra la magnitud de Hdd en un gráfico 3D
        view([-9.441 62.412]);   % Ángulo de vista para el gráfico 3D
        
        title('Delay-Doppler Channel Response H_{dd} from Channel Sounding');
        xlabel('Normalized Doppler');
        ylabel('Normalized Delay');
        zlabel('Magnitude');
        
        %guardar grafica 
        savefig(gcf, fullfile(fullFolderPath, 'hdd_channel_response.fig'));  % Guarda en la carpeta creada

        % Asegurarse de que solo se grafique una vez
        channel_response = false;
    end
    SINR_slot_dB =  pow2db(SINR_Slot/numb_slots);
    SINR_W = SINR_Slot+SINR_W;
    %print constelations
    if images_turn && mod(i_snr,3)==0
        % Llamada a la función 
        plotConstellationSubplot(Xgrid, Xhat_otfs, SNRdB,fullFolderPath,SINR_slot_dB);
    end

    berOTFS_values(i_snr) = n_errors/(num_tx_bits*numb_slots);
    blerOTFS_values(i_snr) = n_blk_errors/numb_slots;
       
end
%parte eficiencia espectral 
% Encuentra el índice del valor más cercano a 12 dB
%[~, idx_snr] = min(abs(SNRdB_values - 16));  % Calcula la diferencia y encuentra el índice mínimo

% Ahora puedes obtener el valor de BER correspondiente a ese índice
% ber_fijo = berOTFS_values(idx_snr);
% fprintf('ber_fijo %.2f \n',ber_fijo);
% effSpectral = (log2(mod_order) * M) / (M + padLen); % bps/Hz
% fprintf('effSpectral %.2f \n',effSpectral);

%Calculate the final value of SNRI
SINR_dB = pow2db(SINR_W/(numb_slots*length(SNRdB_values)));
fprintf('SINR_dB %.2f \n',SINR_dB);

archivo = fullfile(fullFolderPath, 'SINR_dB.mat');
save(archivo, 'SINR_dB');  % Guardar la variable 'valor' en el archivo .mat

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
savefig(gcf, fullfile(fullFolderPath, 'BER_vs_SNR.fig'));  % Guarda en la carpeta creada

% Graficar la BLER contra los valores de SNR
figure;
semilogy(SNRdB_values, blerOTFS_values, '-o', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('BLER');
title('BLER vs. SNR para OTFS');
grid on;
ylim([1e-1, 1]); % Ajusta los límites del eje Y (ejemplo)
savefig(gcf, fullfile(fullFolderPath, 'BLER_vs_SNR.fig'));  % Guarda en la carpeta creada

% Guardar todas las variables en un archivo .mat dentro de la carpeta creada
save(fullfile(fullFolderPath, 'resultados.mat'));

% Mensaje de confirmación
fprintf('Resultados guardados en: %s\n', fullFolderPath);