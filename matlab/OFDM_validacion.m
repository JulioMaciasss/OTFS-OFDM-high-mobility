% OFDM_validacion.m — OFDM validation (AWGN + optional Doppler channel)
% Reproducible version (uses CFG when available), fixed BLER logic

clearvars -except CFG;
close all; clc;

%% Paths
thisDir = fileparts(mfilename("fullpath"));
addpath(fullfile(thisDir,"utils"));

%% ---- Load central config if available ----
if exist("CFG","var")
    M = CFG.M; N = CFG.N; df = CFG.df; fc = CFG.fc; %#ok<NASGU>
    padLen = CFG.padLen_default;
    padType = char(CFG.padType);

    SNRdBvalues = CFG.SNRdBvalues;
    numbslots  = CFG.numbslots;

    % Validation is AWGN by default
    channel_type = "awgn";

    % If you want Doppler test inside validation, set this to true manually:
    useDoppler = false;

    basePath = CFG.out.results_dir;

    % Doppler params (only used if useDoppler==true)
    chanParams.pathDelays   = CFG.doppler.pathDelays;
    chanParams.pathGains    = CFG.doppler.pathGains;
    chanParams.pathDopplers = CFG.doppler.pathDopplers;
else
    M = 64; N = 30; df = 15e3; fc = 3e9; %#ok<NASGU>
    padLen = 10; padType = 'ZP';

    SNRdBvalues = linspace(0,20,12);
    numbslots  = 100;

    channel_type = "awgn";
    useDoppler = false;

    basePath = fullfile(pwd,"results");

    chanParams.pathDelays   = [0 5 8];
    chanParams.pathGains    = [1 0.7 0.5];
    chanParams.pathDopplers = [0 0 1];
end
% ------------------------------------------

%% Derived
fsamp = M*df;
Meff = M + padLen;
numSamps = Meff * N; %#ok<NASGU>
T = (M+padLen)/(M*df); %#ok<NASGU>

chanParams.pathDopplerFreqs = chanParams.pathDopplers * (1/(N*T));

%% Output folder
folderName = fullfile(basePath, sprintf("OFDM_validation_%s_padLen%d", channel_type, padLen));
if ~exist(folderName,'dir'); mkdir(folderName); end

%% Pilot transmission (full grid pilot)
P = exp(1i*pi/4);
txOut0 = ofdmmod(P*ones(M,N), M, padLen);

%% Energy / counters
Es = computeSymbolEnergy(4); % QPSK reference energy
numtxbits = 2*M*N;           % because QPSK => 2 bits/sym

berOFDMvalues  = zeros(size(SNRdBvalues));
blerOFDMvalues = zeros(size(SNRdBvalues));

%% SNR loop
for isnr = 1:length(SNRdBvalues)
    SNRdB = SNRdBvalues(isnr);

    nerrors = 0;
    nblkerrors = 0;

    n0 = Es/(10^(SNRdB/10));

    % Channel estimate (only if using Doppler channel)
    if useDoppler && channel_type ~= "awgn"
        dopplerOut1 = dopplerChannel(txOut0, fsamp, chanParams);
        chOutPilot  = awgn(dopplerOut1, SNRdB, 'measured');
        Ypilot = ofdmdemod(chOutPilot(1:(M+padLen)*N), M, padLen);
        Hofdm  = Ypilot * conj(P) / (abs(P)^2 + n0);
    end

    for islot = 1:numbslots
        % Data generation (QPSK bits)
        Xdata = randi([0,1], 2*M, N);
        Xgrid = pskmod(Xdata, 4, pi/4, InputType="bit");

        txOut = ofdmmod(Xgrid, M, padLen);

        if useDoppler && channel_type ~= "awgn"
            dopplerOut2 = dopplerChannel(txOut, fsamp, chanParams);
        else
            dopplerOut2 = txOut;
        end

        chOut = awgn(dopplerOut2, SNRdB, 'measured');
        rxWindow = chOut(1:(M+padLen)*N);

        if useDoppler && channel_type ~= "awgn"
            Yofdm = ofdmdemod(rxWindow, M, padLen);
            Xhat = conj(Hofdm).*Yofdm./(abs(Hofdm).^2 + n0);
        else
            Xhat = ofdmdemod(rxWindow, M, padLen);
        end

        XhatData = pskdemod(Xhat, 4, pi/4, OutputType="bit", OutputDataType="logical");

        biterrors = xor(Xdata, XhatData);
        erroresahora = sum(biterrors(:));

        nerrors = nerrors + erroresahora;
        nblkerrors = nblkerrors + (erroresahora > 0);
    end

    berOFDMvalues(isnr)  = nerrors/(numtxbits*numbslots);
    blerOFDMvalues(isnr) = nblkerrors/numbslots;
end

%% Avoid log(0)
berOFDMvalues(berOFDMvalues==0) = 1e-10;
blerOFDMvalues(blerOFDMvalues==0) = 1e-10;

%% Plots + save
fig1 = figure;
semilogy(SNRdBvalues, berOFDMvalues, '-o', 'LineWidth', 2);
xlabel('SNR (dB)'); ylabel('BER');
title('OFDM validation — BER vs SNR');
grid on; ylim([1e-4 1]);
savefig(fig1, fullfile(folderName, 'BER_vs_SNR.fig'));
exportgraphics(fig1, fullfile(folderName, 'BER_vs_SNR.png'));

fig2 = figure;
semilogy(SNRdBvalues, blerOFDMvalues, '-o', 'LineWidth', 2);
xlabel('SNR (dB)'); ylabel('BLER');
title('OFDM validation — BLER vs SNR');
grid on; ylim([1e-1 1]);
savefig(fig2, fullfile(folderName, 'BLER_vs_SNR.fig'));
exportgraphics(fig2, fullfile(folderName, 'BLER_vs_SNR.png'));

save(fullfile(folderName, 'resultados.mat'));
fprintf("Resultados guardados en: %s\n", folderName);
