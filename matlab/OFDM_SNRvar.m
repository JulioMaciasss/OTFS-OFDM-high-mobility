% ---- Load central config if available ----
if exist("CFG","var")
    M = CFG.M; N = CFG.N; df = CFG.df; fc = CFG.fc;
    SNRdBvalues = CFG.SNRdBvalues; numbslots = CFG.numbslots;
    modorder = CFG.modorder; idealCE = CFG.idealCE;
    imagesturn = CFG.imagesturn; channelresponse = CFG.channelresponse;
    padType = char(CFG.padType);
end
% ------------------------------------------


% SIMULATION SETUP
M = 64;          % Number of subcarriers
N = 30;          % Number of subsymbols per frame
df = 15e3;       % LTE subcarrier spacing (15 kHz)
fc = 3e9;        % Carrier frequency (5 GHz)
padLen = 10;     % Number of padding samples, should be greater than channel dispersion
padType = 'ZP';  % Zero Padding (ZP) to mitigate Inter-Symbol Interference (ISI)
SNRdB_values = linspace(0, 20, 12);
numb_slots = 100; % Number of slots, each slot contains 30 symbols
idealCE = true;
mod_order = 4; %2->BPSK, 4-> QPSK, 16-> 16-QAM, 64-> 64-QAM
images_turn = false;
channel_response = false;

% Crear sufijo de nombre según la modulación
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

% Definir el path base dinámicamente según la modulación
basePath = CFG.out.results_dir;

% Pilot generation and grid population
pilotBin = floor(N/2)+1;
Pdd = zeros(M,N);
Pdd(1,pilotBin) = exp(1i*pi/4); % populate just one bin to se

% Configure paths
chanParams.pathDelays      = 0; % number of samples that path is delayed
chanParams.pathGains       = 1; % complex path gain
chanParams.pathDopplers    = -2; % Doppler index as a multiple of fsamp/MN

% Convertir los valores del parámetro en una cadena de texto
dopplerStr = strjoin(string(chanParams.pathDopplers), '_');  % Crear string separado por "_"
% Concatenar la ruta manualmente con una barra invertida (Windows) o slash (Linux/Mac)
folderName = strcat(basePath, '\OFDM_', dopplerStr);  % Para Windows
% Convertir a char para evitar problemas con exist()
folderName = char(folderName);  
% Crear la carpeta si no existe
if ~exist(folderName, 'dir')  
    mkdir(folderName);
end

fsamp = M*df;            % freq muestreo M*df
Meff = M + padLen;       % numero efect muestras por simb M+padlen(padding)
numSamps = Meff * N;     % Num muestras por OTFS sub-simb incluye padding
T = ((M+padLen)/(M*df)); % duracion del simbolo, Meff/numSamps

% Para calcular la freq doppler en hz
chanParams.pathDopplerFreqs = chanParams.pathDopplers * 1/(N*T); % Hz

% Transmit pilots over all subcarriers and symbols to sound the channel
txOut0 = ofdmmod(exp(1i*pi/4)*ones(M,N),M,padLen);       % transmit pilots over the entire grid

%Calculate the symbol energy
Es = computeSymbolEnergy(mod_order);

num_tx_bits = log2(mod_order)*M*N; %numb of tx bits

berOFDM_values = zeros(size(SNRdB_values));
blerOFDM_values = zeros(size(SNRdB_values));

for i_snr=1:length(SNRdB_values)
    SNRdB = SNRdB_values(i_snr);
    n_errors = 0;
    block_error = 0;
    n_blk_errors = 0;
    n0 = Es/(10^(SNRdB/10));
    
    dopplerOut1 = dopplerChannel(txOut0,fsamp,chanParams);% send through channel
     
    if idealCE
        chOut = dopplerOut1;
    else
        chOut = awgn(dopplerOut1,SNRdB,'measured');
    end

    Yofdm = ofdmdemod(chOut(1:(M+padLen)*N),M,padLen);      % demodulate
    P = exp(1i*pi/4);  % el valor transmitido en toda la grilla
    
    if idealCE
        Hofdm = Yofdm * conj(P) / (abs(P)^2);
    else
        Hofdm = Yofdm * conj(P) / (abs(P)^2 + n0);
    end
    
    % % Para graficar la respuesta del canal, solo grafico el ultimo con
    % bajo ruido
    if channel_response && i_snr == length(SNRdB_values)
        x = 1:M;  % Índices para las subportadoras (frecuencia Doppler)
        y = 1:N;  % Índices para los símbolos (retardo)

        % Visualización de la magnitud de Hofdm con ejes intercambiados
        figure;
        surf(y, x, abs(Hofdm));  % Se intercambian los ejes para representar H(t,f)

        % Mejorar la visualización
        shading interp;  % Interpolación para suavizar el color
        colormap('jet');  % Usar un mapa de colores para mejorar el contraste
        colorbar;  % Mostrar la barra de colores para la magnitud

        % Etiquetas y título corregidos
        xlabel('Time Index (t)');  % Representación del tiempo en el eje X
        ylabel('Frequency Index (f)');  % Representación de la frecuencia en el eje Y
        zlabel('Magnitude');
        title('Channel Estimation H(t,f) - Direct Representation');
    end
    
        for i_slot = 1:numb_slots
            % Data generation
            Xgrid = zeros(M,N);
            Xdata = randi([0,1],log2(mod_order)*M,N);
            %modulo
            Xgrid(1:M,:) = helperModulator(Xdata, mod_order);
        
            % Transmit data over the same channel and use channel estimates to equalize
            txOut = ofdmmod(Xgrid,M,padLen);                        % transmit data grid
            
            dopplerOut2 = dopplerChannel(txOut,fsamp,chanParams);    % send through channel
            chOut = awgn(dopplerOut2,SNRdB,'measured');              % add noise
            
            rxWindow = chOut(1:(M+padLen)*N);
            Yofdm = ofdmdemod(rxWindow,M,padLen);                   % demodulate
            if idealCE
                Xhat_ofdm = conj(Hofdm) .* Yofdm ./ (abs(Hofdm).^2); % equalize with LMMSE
            else
                Xhat_ofdm = conj(Hofdm) .* Yofdm ./ (abs(Hofdm).^2+n0); % equalize with LMMSE
            end
            %demodulo
            XhatDataOFDM = helperDemodulator(Xhat_ofdm, mod_order);
            
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

        %Hago print del ultimo slot 
        if images_turn && mod(i_snr,3)==0
            % Crear constelación de referencia según mod_order
            if mod_order == 2
                refConst = pskmod(0:mod_order-1, mod_order, pi);
            elseif mod_order == 4
                refConst = pskmod(0:mod_order-1, mod_order, pi/4);
            elseif mod_order == 16 || mod_order == 64
                refConst = qammod(0:mod_order-1, mod_order);
            else
                error('Modulación no soportada para la constelación');
            end
            title_str = sprintf('OFDM, SNR = %.2f dB', SNRdB);
            constDiagOFDM = comm.ConstellationDiagram( ...
                'ReferenceConstellation', refConst, ...
                'XLimits', [-2 2], ...
                'YLimits', [-2 2],...
                'Title',title_str);
            
            constDiagOFDM(Xhat_ofdm(:));
            %saveas(gcf, fullfile(folderName, strcat(title_str, '.png')));  % Guarda en la carpeta creada
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

