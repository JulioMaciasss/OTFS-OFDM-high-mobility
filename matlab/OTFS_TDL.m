% Definir el path base donde quieres guardar los resultados
basePath = 'C:\Users\jm-y2\Desktop\TFG\SIMULACIONES\OTFS_3GHz_tdl';  % Cambia esta ruta según tu preferencia

%SIMULATION SET UP
M = 64;          % num subportadoras
N = 30;          % number of subsymb por trama
df = 15e3;       % espaciado de freq LTE
fc = 3e9;        % freq portadora
padLen = 10;     % num muestras relleno, mayor q dispersion del canal
padType = 'ZP';  % zero padding pa mitigar ISI
SNRdB_values = linspace(0, 20, 12);
numb_slots = 100; %cada slot 30 simb
idealCE =  false;
mod_order = 4; %2->BPSK, 4-> QPSK, 16-> 16-QAM, 64-> 64-QAM
images_turn = false;
channel_response = false;

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

% Modelo TDL
channel = nrTDLChannel;

% Configuración común para el canal
channel.SampleRate = fsamp;
% Pedir al usuario el nivel de movilidad
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
folderName = strcat(basePath, '\OTFS_', mobilityLevelStr);  % Para Windows
% Convertir a char para evitar problemas con exist()
folderName = char(folderName);  
% Crear la carpeta si no existe
if ~exist(folderName, 'dir')  
    mkdir(folderName);
end

Es = computeSymbolEnergy(mod_order);

num_tx_bits = log2(mod_order)*M*N;

berOTFS_values = zeros(size(SNRdB_values));
blerOTFS_values = zeros(size(SNRdB_values));


for i_snr=1:length(SNRdB_values)
    SNRdB=SNRdB_values(i_snr);
    n_errors = 0;
    block_error = 0;
    n_blk_errors = 0;
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
        
        % Filtrar los valores de Hdd que sean suficientemente grandes (umbral 0.05)
        [lp, vp] = find(abs(Hdd) >= 0.05);  % Encontrar índices donde la magnitud es suficientemente grande
        
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
        chOut = awgn(dopplerOut2,SNRdB,'measured');%dopplerout es la señal a la q añado
            %ruido con la realcion señal ruido, y measured indica que debe calcular
            %automente la potencia de la señal antes q el ruid
        
        % Form G matrix using channel estimates
        G = getG(M,N,chanEst,padLen,padType);
    
        rxWindow = chOut(1:numSamps,1);
        
        %ESTIMO CON LMMSE 
        y_otfs = ((G'*G)+n0*eye(Meff*N)) \ (G'*rxWindow); % LMMSE
        
        %DEMODULO Y DETECTO BITS 
        Xhat_otfs = helperOTFSdemod(y_otfs,M,padLen,0,padType); % OTFS demodulation
        %demodulo
        XhatDataOTFS = helperDemodulator(Xhat_otfs, mod_order);

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
        
        % Asegurarse de que solo se grafique una vez
        channel_response = false;
    end

    %print constelations
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

        constDiagOTFS = comm.ConstellationDiagram( ...
            'ReferenceConstellation', refConst, ...
            'XLimits', [-2 2], ...
            'YLimits', [-2 2], ...
            'Title', 'OTFS');
        subplot(2, 3, i_snr);
        constDiagOTFS(Xhat_otfs(:));
        title(['Simulación ' num2str()]); % Opcional: título para cada gráfico
    end

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
savefig(gcf, fullfile(folderName, 'BER_vs_SNR.fig'));  % Guarda en la carpeta creada

% Graficar la BLER contra los valores de SNR
figure;
semilogy(SNRdB_values, blerOTFS_values, '-o', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('BLER');
title('BLER vs. SNR para OTFS');
grid on;
ylim([1e-1, 1]); % Ajusta los límites del eje Y (ejemplo)
savefig(gcf, fullfile(folderName, 'BLER_vs_SNR.fig'));  % Guarda en la carpeta creada

% Guardar todas las variables en un archivo .mat dentro de la carpeta creada
save(fullfile(folderName, 'resultados.mat'));

% Mensaje de confirmación
fprintf('Resultados guardados en: %s\n', folderName);