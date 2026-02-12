% OFDM_SNRvar.m — OFDM over sparse delay–Doppler (custom dopplerChannel)
% Reproducible version (uses CFG when available)

clearvars -except CFG;
close all; clc;

%% Paths
thisDir = fileparts(mfilename("fullpath"));
addpath(fullfile(thisDir,"utils"));

%% ---- Load central config if available ----
if exist("CFG","var")
    basePath = CFG.out.results_dir;

    M = CFG.M; N = CFG.N; df = CFG.df; fc = CFG.fc; %#ok<NASGU>
    padLen = CFG.padLen_default;
    padType = char(CFG.padType);

    SNRdBvalues = CFG.SNRdBvalues;
    numbslots = CFG.numbslots;

    idealCE = CFG.idealCE;
    modorder = CFG.modorder;
    imagesturn = CFG.imagesturn;
    channelresponse = CFG.channelresponse;

    % Doppler channel params
    chanParams.pathDelays   = CFG.doppler.pathDelays;
    chanParams.pathGains    = CFG.doppler.pathGains;
    chanParams.pathDopplers = CFG.doppler.pathDopplers;
else
    basePath = fullfile(pwd,"results");

    M = 64; N = 30; df = 15e3; fc = 3e9; %#ok<NASGU>
    padLen = 10; padType = 'ZP';

    SNRdBvalues = linspace(0,20,12);
    numbslots = 100;

    idealCE = true;
    modorder = 4;

    imagesturn = false;
    channelresponse = false;

    chanParams.pathDelays   = 0;
    chanParams.pathGains    = 1;
    chanParams.pathDopplers = -2;
end
% ------------------------------------------

%% Derived
fsamp = M*df;
Meff = M + padLen;
numSamps = Meff * N; %#ok<NASGU>
T = (M+padLen)/(M*df); %#ok<NASGU>

% For reference (Hz)
chanParams.pathDopplerFreqs = chanParams.pathDopplers * (1/(N*T));

%% Output folder
dopplerStr = strjoin(string(chanParams.pathDopplers), "_");
folderName = fullfile(basePath, sprintf("OFDM_Doppler_k%s_padLen%d_mod%d", dopplerStr, padLen, modorder));
if ~exist(folderName,'dir'); mkdir(folderName); end

%% Pilot transmission (full grid pilot) + pilot symbol
P = exp(1i*pi/4);
txOut0 = ofdmmod(P*ones(M,N), M, padLen);

%% Energy / counters
Es = computeSymbolEnergy(modorder);
numtxbits = log2(modorder)*M*N;

berOFDMvalues  = zeros(size(SNRdBvalues));
blerOFDMvalues = zeros(size(SNRdBvalues));

%% SNR loop
for isnr = 1:length(SNRdBvalues)
    SNRdB = SNRdBvalues(isnr);

    nerrors = 0;
    nblkerrors = 0;

    n0 = Es/(10^(SNRdB/10));

    % --- channel estimation (same channel for whole SNR point) ---
    dopplerOut1 = dopplerChannel(txOut0, fsamp, chanParams);

    if idealCE
        chOutPilot = dopplerOut1;
    else
        chOutPilot = awgn(dopplerOut1, SNRdB, 'measured');
    end

    Ypilot = ofdmdemod(chOutPilot(1:(M+padLen)*N), M, padLen);

    if idealCE
        Hofdm = Ypilot * conj(P) / (abs(P)^2);
    else
        Hofdm = Ypilot * conj(P) / (abs(P)^2 + n0);
    end

    % Optional channel response plot (once)
    if channelresponse && (isnr == length(SNRdBvalues))
        figure;
        surf(1:N, 1:M, abs(Hofdm));
        shading interp; colormap('jet'); colorbar;
        xlabel('Time Index (t)'); ylabel('Frequency Index (f)'); zlabel('Magnitude');
        title('Channel Estimation H(t,f) - Direct Representation');
        savefig(gcf, fullfile(folderName, 'channel_response.fig'));
        exportgraphics(gcf, fullfile(folderName, 'channel_response.png'));
    end

    % --- Monte Carlo slots ---
    for islot = 1:numbslots
        Xdata = randi([0,1], log2(modorder)*M, N);
        Xgrid = helperModulator(Xdata, modorder);

        txOut = ofdmmod(Xgrid, M, padLen);

        dopplerOut2 = dopplerChannel(txOut, fsamp, chanParams);

        if idealCE
            chOut = dopplerOut2;
        else
            chOut = awgn(dopplerOut2, SNRdB, 'measured');
        end

        rxWindow = chOut(1:(M+padLen)*N);
        Yofdm = ofdmdemod(rxWindow, M, padLen);

        if idealCE
            Xhatofdm = conj(Hofdm).*Yofdm./(abs(Hofdm).^2);
        else
            Xhatofdm = conj(Hofdm).*Yofdm./(abs(Hofdm).^2 + n0);
        end

        XhatDataOFDM = helperDemodulator(Xhatofdm, modorder);

        biterrors = xor(Xdata, XhatDataOFDM);
        erroresahora = sum(biterrors(:));

        nerrors = nerrors + erroresahora;
        nblkerrors = nblkerrors + (erroresahora > 0);

        if imagesturn && mod(isnr,3)==0 && islot==numbslots
            figure; plot(real(Xhatofdm(:)), imag(Xhatofdm(:)), '.'); grid on;
            title(sprintf('OFDM constellation (SNR=%.1f dB)', SNRdB));
        end
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
title('OFDM — BER vs SNR (Sparse Doppler Channel)');
grid on; ylim([1e-4 1]);
savefig(fig1, fullfile(folderName, 'BER_vs_SNR.fig'));
exportgraphics(fig1, fullfile(folderName, 'BER_vs_SNR.png'));

fig2 = figure;
semilogy(SNRdBvalues, blerOFDMvalues, '-o', 'LineWidth', 2);
xlabel('SNR (dB)'); ylabel('BLER');
title('OFDM — BLER vs SNR (Sparse Doppler Channel)');
grid on; ylim([1e-1 1]);
savefig(fig2, fullfile(folderName, 'BLER_vs_SNR.fig'));
exportgraphics(fig2, fullfile(folderName, 'BLER_vs_SNR.png'));

save(fullfile(folderName, 'resultados.mat'));

fprintf("Resultados guardados en: %s\n", folderName);
