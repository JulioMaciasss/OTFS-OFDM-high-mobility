% OTFS_SNRvariable.m — OTFS over sparse delay–Doppler (custom dopplerChannel)
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
    padLen  = CFG.padLen_default;
    padType = char(CFG.padType);

    SNRdBvalues = CFG.SNRdBvalues;
    numbslots  = CFG.numbslots;

    idealCE = CFG.idealCE;
    modorder = CFG.modorder;

    imagesturn = CFG.imagesturn;
    channelresponse = CFG.channelresponse;

    threshold = CFG.beta;

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

    idealCE = false;
    modorder = 4;

    imagesturn = false;
    channelresponse = false;

    threshold = 0.05;

    chanParams.pathDelays   = 0;
    chanParams.pathGains    = 1;
    chanParams.pathDopplers = -2;
end
% ------------------------------------------

%% Derived
fsamp = M*df;
Meff = M + padLen;
numSamps = Meff * N;
T = (M+padLen)/(M*df); %#ok<NASGU>

% For reference (Hz)
chanParams.pathDopplerFreqs = chanParams.pathDopplers * (1/(N*T));

%% Output folder
dopplerStr = strjoin(string(chanParams.pathDopplers), "_");
folderName = fullfile(basePath, sprintf("OTFS_Doppler_k%s_padLen%d_mod%d_beta%.4g", dopplerStr, padLen, modorder, threshold));
if ~exist(folderName,'dir'); mkdir(folderName); end

%% Pilot in delay–Doppler (single pilot)
pilotBin = floor(N/2)+1;
Pdd = zeros(M,N);
Pdd(1,pilotBin) = exp(1i*pi/4);

txOut0 = helperOTFSmod(Pdd, padLen, padType);

%% Energy / counters
Es = computeSymbolEnergy(modorder);
numtxbits = log2(modorder)*M*N;

berOTFSvalues  = zeros(size(SNRdBvalues));
blerOTFSvalues = zeros(size(SNRdBvalues));

%% SNR loop
for isnr = 1:length(SNRdBvalues)
    SNRdB = SNRdBvalues(isnr);

    nerrors = 0;
    nblkerrors = 0;

    n0 = Es/(10^(SNRdB/10));

    % ---- Channel estimation (pilot) ----
    dopplerOut1 = dopplerChannel(txOut0, fsamp, chanParams);

    if idealCE
        chOutPilot = dopplerOut1;
    else
        chOutPilot = awgn(dopplerOut1, SNRdB, 'measured');
    end

    rxIn = chOutPilot(1:numSamps);
    Ydd = helperOTFSdemod(rxIn, M, padLen, 0, padType);

    if idealCE
        Hdd = Ydd * conj(Pdd(1,pilotBin)) / (abs(Pdd(1,pilotBin))^2);
    else
        Hdd = Ydd * conj(Pdd(1,pilotBin)) / (abs(Pdd(1,pilotBin))^2 + n0);
    end

    % Optional channel response plot (once)
    if channelresponse && (isnr == 1)
        figure;
        xa = -N/2 : N/2 - 1;
        ya = 0:M-1;
        mesh(xa, ya, abs(Hdd));
        view([-9.441 62.412]);
        title('Delay-Doppler Channel Response H_{dd} from Channel Sounding');
        xlabel('Normalized Doppler'); ylabel('Normalized Delay'); zlabel('Magnitude');
        savefig(gcf, fullfile(folderName, 'hdd_channel_response.fig'));
        exportgraphics(gcf, fullfile(folderName, 'hdd_channel_response.png'));
    end

    % Thresholding in delay–Doppler
    [lp, vp] = find(abs(Hdd) >= threshold);

    chanEst.pathGains    = diag(Hdd(lp,vp));
    chanEst.pathDelays   = lp - 1;
    chanEst.pathDopplers = vp - pilotBin;

    % ---- Monte Carlo slots ----
    for islot = 1:numbslots
        % Data generation
        Xdata = randi([0,1], log2(modorder)*M, N);
        Xgrid = helperModulator(Xdata, modorder);

        % OTFS modulation
        txOut = helperOTFSmod(Xgrid, padLen, padType);

        % Channel + noise
        dopplerOut2 = dopplerChannel(txOut, fsamp, chanParams);
        if idealCE
            chOut = dopplerOut2;
        else
            chOut = awgn(dopplerOut2, SNRdB, 'measured');
        end

        % Form G matrix from estimates
        G = getG(M, N, chanEst, padLen, padType);

        rxWindow = chOut(1:numSamps);

        % LMMSE equalization (time domain)
        y_otfs = ((G'*G) + n0*eye(Meff*N)) \ (G' * rxWindow);

        % OTFS demod + bit detect
        Xhat_otfs = helperOTFSdemod(y_otfs, M, padLen, 0, padType);
        XhatDataOTFS = helperDemodulator(Xhat_otfs, modorder);

        % Error counting
        biterrors = xor(Xdata, XhatDataOTFS);
        erroresahora = sum(biterrors(:));

        nerrors = nerrors + erroresahora;
        nblkerrors = nblkerrors + (erroresahora > 0);

        % Optional constellation
        if imagesturn && mod(isnr,3)==0 && islot==numbslots
            figure; plot(real(Xhat_otfs(:)), imag(Xhat_otfs(:)), '.'); grid on;
            title(sprintf('OTFS constellation (SNR=%.1f dB)', SNRdB));
        end
    end

    berOTFSvalues(isnr)  = nerrors/(numtxbits*numbslots);
    blerOTFSvalues(isnr) = nblkerrors/numbslots;
end

%% Avoid log(0)
berOTFSvalues(berOTFSvalues==0) = 1e-10;
blerOTFSvalues(blerOTFSvalues==0) = 1e-10;

%% Plots + save
fig1 = figure;
semilogy(SNRdBvalues, berOTFSvalues, '-o', 'LineWidth', 2);
xlabel('SNR (dB)'); ylabel('BER');
title('OTFS — BER vs SNR (Sparse Doppler Channel)');
grid on; ylim([1e-4 1]);
savefig(fig1, fullfile(folderName, 'BER_vs_SNR.fig'));
exportgraphics(fig1, fullfile(folderName, 'BER_vs_SNR.png'));

fig2 = figure;
semilogy(SNRdBvalues, blerOTFSvalues, '-o', 'LineWidth', 2);
xlabel('SNR (dB)'); ylabel('BLER');
title('OTFS — BLER vs SNR (Sparse Doppler Channel)');
grid on; ylim([1e-1 1]);
savefig(fig2, fullfile(folderName, 'BLER_vs_SNR.fig'));
exportgraphics(fig2, fullfile(folderName, 'BLER_vs_SNR.png'));

save(fullfile(folderName, 'resultados.mat'));

fprintf("Resultados guardados en: %s\n", folderName);
