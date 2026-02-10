%%------------------------------------------------------------------------------------------------------------%%
% Eng.Mahmoud Aldababseh
% Research and Teaching assistant in Al-Quds University in Palestine
% eng.mahmoud_87@hotmail.comm
%%------------------------------------------------------------------------------------------------------------%%
% This code show BER vs SNR for OFDM under AWGN condition
%%------------------------------------------------------------------------------------------------------------%%
    clear all
    clc
% Initialization
    M  = 2048;            % Number of subchannels
    Guard_Interval = M/4;% Guard interval length
    q  = 4;              % Modulation level
    nloop = 100;          % Number of iteration  
    SNR= 0:20;           % signal to noise ratio vector in dB
    Rate0= zeros(1, length(SNR)); % initializing bit error rate
%%
for( snr = 1 : length(SNR))
    snr
    for(i = 1 : nloop) 
%% Transmitter
        %% Data generation
            Data_In = randi(q-1, M, 1);
        %% Modualtion
            Data_Mod = qammod(Data_In,q);
        %% IFFT
            Data_IFFT = ifft(Data_Mod,M);
        %% Guard interval insertion
            Data_Guard = [Data_IFFT(M- Guard_Interval + 1 : M);Data_IFFT];
        %% P/S conversion
        [a b]=size(Data_Guard);
        Data_Tx=reshape(Data_Guard,a*b,1);
%% Channel Effect
        %% AWGN 
        Data_Noise = awgn(Data_Tx,snr,'measured');
        % signal_in_dB=10*log10(std(Data_Tx)^2);
    %noise_in_dB=signal_in_dB-SNR(snr); 
    %noise=(10^(noise_in_dB/10)^(1/2))*randn(size( Data_Tx, 1),size( Data_Tx, 2)); 
    %Data_Noise=Data_Tx+noise;
%% Receiver
        %% S/P conversion
        Data_Rx=reshape(Data_Noise,a,b);
        %% Guard interval removal
        Data_Removeal_Gaurd  = Data_Rx(Guard_Interval+1:M+Guard_Interval);
        %% FFT
         Data_FFT   = fft(Data_Removeal_Gaurd,M);
        %% Demodualtion
        Data_Out   = qamdemod(Data_FFT,q);
%% BER
        % 
        [nErr rate0] = symerr(Data_Out,Data_In);
        Rate0(snr)= Rate0(snr) + rate0;
    end
        Rate0(snr)= Rate0(snr)/nloop;
end

hf = figure;
semilogy(SNR,Rate0,'r-o')
hold on
xlabel('Signal to Noise Ratio [SNR] in dB');
ylabel('Bit Error Rate [BER]')
grid on
