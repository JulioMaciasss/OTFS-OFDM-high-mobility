% Definición de la función Q usando normcdf
Q = @(x) 1 - normcdf(x);

% Valores de SNR en dB
SNRdB_values = linspace(0, 20, 100); 
SNR_linear = 10.^(SNRdB_values/10); % Pasamos a escala lineal

M = 4; % QPSK
EbN0 = SNR_linear / log2(M); 

% BER teórica para QPSK en AWGN
ber_theoretical = 2 * Q(sqrt(2*EbN0)); 

% Gráfica
figure;
semilogy(SNRdB_values, ber_theoretical, '--r', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('BER');
title('BER Teórica de QPSK en canal AWGN');
grid on;
ylim([1e-4, 1]);
