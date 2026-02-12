% OTFS_CDL.m — OTFS over 3GPP CDL-A channel (reproducible)
% Uses CFG when available, portable paths, fixed BLER logic

clearvars -except CFG;
close all; clc;

%% Paths
thisDir = fileparts(mfilename("fullpath"));
addpath(fullfile(thisDir,"utils"));

%% ---- Load central config if available ----
if exist("CFG","var")
    basePath = CFG.out.results_dir;

    M = CFG.M; N = CFG.N; df = CFG.df; fc = CFG.fc; %#ok<NASGU>
    padLen = CFG.padLen_cdl;
    padType = char(CFG.padType);

    SNRdBvalues = CFG.SNRdBvalues;
    numbslots  = CFG.numbslots;

    idealCE = CFG.idealCE;
    modorder = CFG.modorder;

    imagesturn = CFG.imagesturn;
    channelresponse = CFG.channelresponse;

    threshold = CFG.beta;

    cdlProfile = char(CFG.cdl.DelayProfile);
    DF = CFG.cdl.MaximumDopplerShift;
    DS = CFG.cdl.DelaySpread;
else
    basePath = fullfile(pwd,"results");

    M = 64; N = 30; df = 15e3; fc = 3e9; %#ok<NASGU>
    padLen = 30; padType = 'ZP';

    SNRdBvalues = linspace(0,20,12);
    numbslots  = 100;

    idealCE = false;
    modorder = 4;

    imagesturn = false;
    channelresponse = true;

    threshold = 0.0305;

    cdlProfile = "CDL-A";
    DF = 3;
    DS = 30e-9;
end
% ------------------------------------------

%% Mod name (only for labeling)
switch modorder
    case 2,  modName = "BPSK";
    case 4,  modName = "QPSK";
    case 16, modName = "16QAM";
    case 64, modName = "64QAM";
    otherwise, error("Unsupported modulation order: %d", modorder);
end

%% Derived
fsamp = M*df;
Meff = M + padLen;
numSamps = Meff * N;
T = (M+padLen)/(M*df); %#ok<NASGU>

%% Output folder
DS_ns = round(DS*1e9);
tagIdeal = "";
if idealCE; tagIdeal = "_idealCE"; end
fullFolderPath = fullfile(basePath, sprintf("OTFS_%s_%s_DF%dHz_DS%dns_padLen%d_beta%.4g%s", ...
    cdlProfile, modName, round(DF), DS_ns, padLen, threshold, tagIdeal));
if ~exist(fullFolderPath,'dir'); mkdir(fullFolderPath); end

%% Pilot in delay–Doppler (single pilot)
pilotBin = floor(N/2)+1;
Pdd = zeros(M,N);
Pdd(1,pilotBin) = exp(1i*pi/4);

txOut0 = helperOTFSmod(Pdd, padLen, padType);

%% Channel (CDL)
channel = nrCDLChannel;
channel.SampleRate = fsamp;
channel.DelayProfile = cdlProfile;
channel.MaximumDopplerShift = DF;
channel.DelaySpread = DS;

% SISO arrays (leave default unless you really need to set them)
% NOTE: In your old code you used "cdl.TransmitAntennaArray", but the object is "channel".
channel.TransmitAntennaArray.Size = [1 1 1 1 1];
channel.ReceiveAntennaArray.Size  = [1 1 1 1 1];

%% Energy / counters
Es = computeSymbolEnergy(modorder);
numtxbits = log2(modorder)*M*N;

berOTFSvalues  = zeros(size(SNRdBvalues));
blerOTFSvalues = zeros(size(SNRdBvalues));

SINR_W = 0;

%% Main SNR loop
for isnr = 1:length(SNRdBvalues)
    SNRdB = SNRdBvalues(isnr);

    nerrors = 0;
    nblkerrors = 0;

    SINR_Slot = 0;
    n0 = Es/(10^(SNRdB/10));

    for islot = 1:numbslots
        release(channel);
        channel.Seed = islot;
        reset(channel);

        % ---- Channel sounding with pilot ----
        dopplerOut1 = channel(txOut0);
        if idealCE
            chOutPilot = dopplerOut1;
        else
            chOutPilot = awgn(dopplerOut1, SNRdB, 'measured');
        end

        rxIn = chOutPilot(1:numSamps, 1);
        Ydd = helperOTFSdemod(rxIn, M, padLen, 0, padType);

        if idealCE
            Hdd = Ydd * conj(Pdd(1, pilotBin)) / (abs(Pdd(1, pilotBin))^2);
        else
            Hdd = Ydd * conj(Pdd(1, pilotBin)) / (abs(Pdd(1, pilotBin))^2 + n0);
        end

        % Thresholding (fixed beta for reproducibility)
        [lp, vp] = find(abs(Hdd) >= threshold);

        chanEst.pathGains    = diag(Hdd(lp, vp));
        chanEst.pathDelays   = lp - 1;
        chanEst.pathDopplers = vp - pilotBin;

        reset(channel);

        % ---- Data generation ----
        Xdata = randi([0,1], log2(modorder)*M, N);
        Xgrid = helperModulator(Xdata, modorder);

        % ---- OTFS modulation ----
        txOut = helperOTFSmod(Xgrid, padLen, padType);

        % ---- Channel + noise ----
        dopplerOut2 = channel(txOut);
        if idealCE
            chOut = dopplerOut2;
        else
            chOut = awgn(dopplerOut2, SNRdB, 'measured');
        end

        % ---- LMMSE equalization in time domain ----
        G = getG(M, N, chanEst, padLen, padType);

        rxWindow = chOut(1:numSamps, 1);
        y_otfs = ((G'*G) + n0*eye(Meff*N)) \ (G' * rxWindow);

        % ---- OTFS demod + bit detect ----
        Xhat_otfs = helperOTFSdemod(y_otfs, M, padLen, 0, padType);
        XhatDataOTFS = helperDemodulator(Xhat_otfs, modorder);

        % SINR estimate
        SINR_now_W = 1/mean(abs(Xhat_otfs(:) - Xgrid(:)).^2);
        SINR_Slot = SINR_Slot + SINR_now_W;

        % Errors
        biterrors = xor(Xdata, XhatDataOTFS);
        erroresahora = sum(biterrors(:));
        nerrors = nerrors + erroresahora;

        % BLER: block error if this slot has any bit error
        nblkerrors = nblkerrors + (erroresahora > 0);
    end

    % Optional: plot Hdd once (first SNR point only)
    if channelresponse && (isnr == 1)
        figure;
        xa = -N/2 : N/2 - 1;
        ya = 0:M-1;
        mesh(xa, ya, abs(Hdd));
        view([-9.441 62.412]);
        title('Delay-Doppler Channel Response H_{dd} from Channel Sounding');
        xlabel('Normalized Doppler'); ylabel('Normalized Delay'); zlabel('Magnitude');
        savefig(gcf, fullfile(fullFolderPath, 'hdd_channel_response.fig'));
        exportgraphics(gcf, fullfile(fullFolderPath, 'hdd_channel_response.png'));
    end

    SINR_W = SINR_W + SINR_Slot;

    berOTFSvalues(isnr)  = nerrors/(numtxbits*numbslots);
    blerOTFSvalues(isnr) = nblkerrors/numbslots;
end

%% Final SINR over all SNR points
SINR_dB = pow2db(SINR_W/(numbslots*length(SNRdBvalues)));
fprintf('SINR_dB %.2f\n', SINR_dB);
save(fullfile(fullFolderPath, 'SINR_dB.mat'), 'SINR_dB');

%% Avoid log(0)
berOTFSvalues(berOTFSvalues==0) = 1e-10;
blerOTFSvalues(blerOTFSvalues==0) = 1e-10;

%% Plots + save
fig1 = figure;
semilogy(SNRdBvalues, berOTFSvalues, '-o', 'LineWidth', 2);
xlabel('SNR (dB)'); ylabel('BER');
title('OTFS — BER vs SNR (CDL)');
grid on; ylim([1e-4 1]);
savefig(fig1, fullfile(fullFolderPath, 'BER_vs_SNR.fig'));
exportgraphics(fig1, fullfile(fullFolderPath, 'BER_vs_SNR.png'));

fig2 = figure;
semilogy(SNRdBvalues, blerOTFSvalues, '-o', 'LineWidth', 2);
xlabel('SNR (dB)'); ylabel('BLER');
title('OTFS — BLER vs SNR (CDL)');
grid on; ylim([1e-1 1]);
savefig(fig2, fullfile(fullFolderPath, 'BLER_vs_SNR.fig'));
exportgraphics(fig2, fullfile(fullFolderPath, 'BLER_vs_SNR.png'));

save(fullfile(fullFolderPath, 'resultados.mat'));
fprintf('Resultados guardados en: %s\n', fullFolderPath);
