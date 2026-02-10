% Definir el path base donde quieres guardar los resultados
basePath = 'C:\Users\jm-y2\Desktop\TFG\SIMULACIONES\OFDM_3GHz_tdl';  % Cambia esta ruta según tu preferencia

% SIMULATION SETUP
M = 64;          % Number of subcarriers
N = 30;          % Number of subsymbols per frame
df = 15e3;       % LTE subcarrier spacing (15 kHz)
fc = 3e9;        % Carrier frequency (5 GHz)
padLen = 10;     % Number of padding samples, should be greater than channel dispersion
padType = 'ZP';  % Zero Padding (ZP) to mitigate Inter-Symbol Interference (ISI)
SNRdB_values = linspace(0, 20, 12);
numb_slots = 100; % Number of slots, each slot contains 30 symbols
idealCE = false;
mod_order = 4; %2->BPSK, 4-> QPSK, 16-> 16-QAM, 64-> 64-QAM
images_turn = false;
channel_response = true;

% Pilot generation and grid population
pilotBin = floor(N/2)+1;
Pdd = zeros(M,N);
Pdd(1,pilotBin) = exp(1i*pi/4); % populate just one bin to se

% Transmit pilots over all subcarriers and symbols to sound the channel
txOut0 = ofdmmod(exp(1i*pi/4)*ones(M,N),M,padLen);       % transmit pilots over the entire grid

fsamp = M*df;            % freq muestreo M*df
Meff = M + padLen;       % numero efect muestras por simb M+padlen(padding)
numSamps = Meff * N;     % Num muestras por OTFS sub-simb incluye padding
T = ((M+padLen)/(M*df)); % duracion del simbolo, Meff/numSamps

% Modelo TDL
channel = nrTDLChannel;

% Configuración común para el canal
channel.SampleRate = fsamp;

mobilityLevel = input(['Seleccione el nivel de movilidad:\n 1: TDL-A\n 2: TDL-B\n 3: TDL-C\n 4: TDL-D\n 5: TDL-E\n Seleccione (1-5): ']);
%'6: TDL-F\n' ...
switch mobilityLevel
    case 1
        channel.DelayProfile = 'TDL-A'; 
        channel.MaximumDopplerShift = 10;    % Usuario caminando
    case 2
        channel.DelayProfile = 'TDL-B'; 
        channel.MaximumDopplerShift = 50;    % Vehículo en ciudad
    case 3
        channel.DelayProfile = 'TDL-C'; 
        channel.MaximumDopplerShift = 100;   % Vehículo rápido
    case 4
        channel.DelayProfile = 'TDL-D'; 
        channel.MaximumDopplerShift = 150;   % Alta dispersión, rural
    case 5
        channel.DelayProfile = 'TDL-E'; 
        channel.MaximumDopplerShift = 200;   % Entorno severo, tren o vehículo muy rápido
    % case 6
    %     channel.DelayProfile = 'TDL-F';
    %     channel.MaximumDopplerShift = 300;   % Tren a alta velocidad, entorno extremo
end

fprintf('Movilidad seleccionada: %d\n', mobilityLevel);

% Convertir los valores del parámetro en una cadena de texto
mobilityLevelStr = strjoin(string(mobilityLevel), '_');  % Crear string separado por "_"
% Concatenar la ruta manualmente con una barra invertida (Windows) o slash (Linux/Mac)
folderName = strcat(basePath, '\OFDM_', mobilityLevelStr);  % Para Windows
% Convertir a char para evitar problemas con exist()
folderName = char(folderName);  
% Crear la carpeta si no existe
if ~exist(folderName, 'dir')  
    mkdir(folderName);
end

Es = computeSymbolEnergy(mod_order);

num_tx_bits = log2(mod_order)*M*N;

berOFDM_values = zeros(size(SNRdB_values));
blerOFDM_values = zeros(size(SNRdB_values));

for i_snr=1:length(SNRdB_values)
    SNRdB = SNRdB_values(i_snr);
    n_errors = 0;
    block_error = 0;
    n_blk_errors = 0;
    n0 = Es/(10^(SNRdB/10));
    
    for i_slot = 1:numb_slots
        release(channel);
        channel.Seed = i_slot;
        reset(channel);
        dopplerOut1 = channel(txOut0);
    
        if idealCE
            chOut = dopplerOut1;
        else
            chOut = awgn(dopplerOut1,SNRdB,'measured');
        end
                      % add noise
        Yofdm = ofdmdemod(chOut(1:(M+padLen)*N,1),M,padLen);      % demodulate
        P = exp(1i*pi/4);  % el valor transmitido en toda la grilla
        if idealCE
            Hofdm = Yofdm * conj(P) / (abs(P)^2);
        else
            Hofdm = Yofdm * conj(P) / (abs(P)^2 + n0);
        end
        
        reset(channel); 

        % Generación de datos
        Xgrid = zeros(M, N);  % Inicialización del grid
        Xdata = randi([0,1], log2(mod_order)*M, N);  % Bits aleatorios
        Xgrid(1:M,:) = helperModulator(Xdata, mod_order);
    
        % Transmisión por canal y ecualización usando estimaciones de canal
        txOut = ofdmmod(Xgrid, M, padLen);  % OFDM modulación
    
        dopplerOut2 = channel(txOut);
        chOut = awgn(dopplerOut2, SNRdB, 'measured');             % Añadir ruido AWGN
    
        % Recepción y demodulación
        rxWindow = chOut(1:(M+padLen)*N,1);         % Tomar muestras útiles
        Yofdm = ofdmdemod(rxWindow, M, padLen);   % OFDM demodulación
    
        % Ecualización LMMSE
        if idealCE
            Xhat_ofdm = conj(Hofdm) .* Yofdm ./ (abs(Hofdm).^2);  % Canal ideal
        else
            Xhat_ofdm = conj(Hofdm) .* Yofdm ./ (abs(Hofdm).^2 + n0);  % Canal ruidoso
        end
    
        %demodulo
        XhatDataOFDM = helperDemodulator(Xhat_ofdm, mod_order);
            
    
        % Contar errores
        bit_errors = xor(Xdata, XhatDataOFDM);
        errores_ahora = sum(bit_errors(:));
        n_errors = n_errors + errores_ahora;
    
        if errores_ahora > 0
            block_error = 1;  % Hay error en el bloque
        else
            block_error = 0;
        end
    
        n_blk_errors = n_blk_errors + block_error;
    end  
    
    % Visualización de la constelación solo cada 3 pasos
    if images_turn && mod(i_snr, 3) == 0
        % Crear constelación de referencia según el orden de modulación
        if mod_order == 2
            refConst = pskmod(0:mod_order-1, mod_order, pi);
        elseif mod_order == 4
            refConst = pskmod(0:mod_order-1, mod_order, pi/4);
        elseif mod_order == 16 || mod_order == 64
            refConst = qammod(0:mod_order-1, mod_order);
        else
            error('Modulación no soportada para la constelación');
        end
    
        % Título de la constelación
        title_str = sprintf('OFDM, SNR = %.2f dB', SNRdB);
    
        % Mostrar diagrama de constelación
        constDiagOFDM = comm.ConstellationDiagram( ...
            'ReferenceConstellation', refConst, ...
            'XLimits', [-2 2], ...
            'YLimits', [-2 2], ...
            'Title', title_str);
    
        constDiagOFDM(Xhat_ofdm(:));
    
        % Puedes descomentar para guardar la imagen si lo necesitas
        % saveas(gcf, fullfile(folderName, strcat(title_str, '.png')));
    end
       
    % Visualizar la respuesta del canal solo para el último SNR (bajo ruido)
    if channel_response && i_snr == length(SNRdB_values)
        x = 1:M;  % Índices de subportadoras (frecuencia Doppler)
        y = 1:N;  % Índices de símbolos (retardo)
    
        % Visualización de la magnitud de H en el dominio tiempo-frecuencia (Hofdm)
        figure;
        surf(y, x, abs(Hofdm));  % Ejes intercambiados para representar H(t,f)
    
        % Mejorar apariencia del gráfico
        shading interp;           % Suavizar los colores
        colormap('jet');          % Mapa de colores más visual
        colorbar;                 % Mostrar barra de magnitudes
    
        % Etiquetas y título
        xlabel('Time Index (t)');         % Eje X = Tiempo
        ylabel('Frequency Index (f)');    % Eje Y = Frecuencia
        zlabel('Magnitude');              % Eje Z = Magnitud del canal
        title('Channel Estimation H(t,f) - Direct Representation');
        savefig(gcf, fullfile(folderName, 'channel_response.fig'));  % Guarda en la carpeta creada
    
    end

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
savefig(gcf, fullfile(folderName, 'BER_vs_SNR.fig'));  % Guarda en la carpeta creada

% Graficar la BLER contra los valores de SNR
figure;
semilogy(SNRdB_values, blerOFDM_values, '-o', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('BLER');
title('BLER vs. SNR para OFDM');
grid on;
ylim([1e-1, 1]); % Ajusta los límites del eje Y (ejemplo)
savefig(gcf, fullfile(folderName, 'BLER_vs_SNR.fig'));  % Guarda en la carpeta creada

% Guardar todas las variables en un archivo .mat dentro de la carpeta creada
save(fullfile(folderName, 'resultados.mat'));

% Mensaje de confirmación
fprintf('Resultados guardados en: %s\n', folderName);
