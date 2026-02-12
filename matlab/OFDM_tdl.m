% OFDM_tdl.m — OFDM over 3GPP TDL channel (reproducible)

clearvars -except CFG;  % keep CFG if run from run_all
close all; clc;

%% Paths
thisDir = fileparts(mfilename("fullpath"));
addpath(fullfile(thisDir,"utils"));

%% ---- Load central config if available ----
if exist("CFG","var")
    basePath = CFG.out.results_dir;

    M = CFG.M; N = CFG.N; df = CFG.df; fc = CFG.fc; %#ok<NASGU>
    SNRdBvalues = CFG.SNRdBvalues; numbslots = CFG.numbslots;

    modorder = CFG.modorder;
    idealCE  = CFG.idealCE;
    imagesturn = CFG.imagesturn;
    channelresponse = CFG.channelresponse;

    padType = char(CFG.padType); %#ok<NASGU>
    padLen  = CFG.padLen_tdl;

    tdlProfile = char(CFG.tdl.DelayProfile);
    DF = CFG.tdl.MaximumDopplerShift;
    DS = CFG.tdl.DelaySpread;
else
    basePath = fullfile(pwd,"results");

    M = 64; N = 30; df = 15e3; fc = 3e9; %#ok<NASGU>
    SNRdBvalues = linspace(0,20,12); numbslots = 100;

    modorder = 4; idealCE = false;
    imagesturn = false; channelresponse = true;

    padType = 'ZP'; %#ok<NASGU>
    padLen = 10;

    tdlProfile = "TDL-A";
    DF = 300;
    DS = 300e-9;
end
% ------------------------------------------

%% Derived
fsamp    = M*df;
Meff     = M + padLen;
numSamps = Meff * N; %#ok<NASGU>

%% Output folder
DS_ns = round(DS*1e9);
folderName = fullfile(basePath, sprintf("OFDM_TDL_%s_DF%dHz_DS%dns_padLen%d_mod%d", ...
    tdlProfile, round(DF), DS_ns, padLen, modorder));
if ~exist(folderName,'dir'); mkdir(folderName); end

%% Pilot transmission (full grid pilot) + pilot symbol
P = exp(1i*pi/4);
txOut0 = ofdmmod(P*ones(M,N), M, padLen);

%% Channel (no prompts)
channel = nrTDLChannel;
channel.SampleRate = fsamp;
channel.DelayProfile = tdlProfile;
channel.MaximumDopplerShift = DF;
channel.DelaySpread = DS;

%% Energy / counters
Es = computeSymbolEnergy(modorder);
numtxbits = log2(modorder)*M*N;

berOFDMvalues  = zeros(size(SNRdBvalues));
blerOFDMvalues = zeros(size(SNRdBvalues));

%% Main SNR loop
for isnr = 1:length(SNRdBvalues)
    SNRdB = SNRdBvalues(isnr);

    nerrors = 0;
    nblkerrors = 0;

    n0 = Es/(10^(SNRdB/10));

    for islot = 1:numbslots
        release(channel);
        channel.Seed = islot;
        reset(channel);

        % --- Channel estimation from pilots ---
        dopplerOut1 = channel(txOut0);
        if idealCE
            chOutPilot = dopplerOut1;
        else
            chOutPilot = awgn(dopplerOut1, SNRdB, 'measured');
        end

        Ypilot = ofdmdemod(chOutPilot(1:(M+padLen)*N,1), M, padLen);

        if idealCE
            Hofdm = Ypilot * conj(P) / (abs(P)^2);
        else
            Hofdm = Ypilot * conj(P) / (abs(P)^2 + n0);
        end

        reset(channel);

        % --- Data generation ---
        Xdata = randi([0,1], log2(modorder)*M, N);
        Xgrid = helperModulator(Xdata, modorder);

        % --- Tx/Rx through channel ---
        txOut = ofdmmod(Xgrid, M, padLen);
        dopplerOut2 = channel(txOut);

        if idealCE
            chOut = dopplerOut2;
        else
            chOut = awgn(dopplerOut2, SNRdB, 'measured');
        end

        rxWindow = chOut(1:(M+padLen)*N,1);
        Yofdm = ofdmdemod(rxWindow, M, padLen);

        % --- LMMSE equalization ---
        if idealCE
            Xhatofdm = conj(Hofdm).*Yofdm./(abs(Hofdm).^2);
        else
            Xhatofdm = conj(Hofdm).*Yofdm./(abs(Hofdm).^2 + n0);
        end

        % --- Demod + error counting ---
        XhatDataOFDM = helperDemodulator(Xhatofdm, modorder);

        biterrors = xor(Xdata, XhatDataOFDM);
        erroresahora = sum(biterrors(:));
        nerrors = nerrors + erroresahora;

        nblkerrors = nblkerrors + (erroresahora > 0);

        % Optional: constellation (rarely)
        if imagesturn && (mod(isnr,3)==0) && (islot==numbslots)
            figure; plot(real(Xhatofdm(:)), imag(Xhatofdm(:)), '.'); grid on;
            title(sprintf("OFDM constellation (SNR=%.1f dB)", SNRdB));
        end
    end

    berOFDMvalues(isnr)  = nerrors/(numtxbits*numbslots);
    blerOFDMvalues(isnr) = nblkerrors/numbslots;
end

%% Avoid log(0) issues
berOFDMvalues(berOFDMvalues==0) = 1e-10;
blerOFDMvalues(blerOFDMvalues==0) = 1e-10;

%% Plots + save
fig1 = figure;
semilogy(SNRdBvalues, berOFDMvalues, '-o', 'LineWidth', 2);
xlabel('SNR (dB)'); ylabel('BER');
title('OFDM — BER vs SNR (TDL)');
grid on; ylim([1e-4 1]);
savefig(fig1, fullfile(folderName, 'BER_vs_SNR.fig'));
exportgraphics(fig1, fullfile(folderName, 'BER_vs_SNR.png'));

fig2 = figure;
semilogy(SNRdBvalues, blerOFDMvalues, '-o', 'LineWidth', 2);
xlabel('SNR (dB)'); ylabel('BLER');
title('OFDM — BLER vs SNR (TDL)');
grid on; ylim([1e-1 1]);
savefig(fig2, fullfile(folderName, 'BLER_vs_SNR.fig'));
exportgraphics(fig2, fullfile(folderName, 'BLER_vs_SNR.png'));

save(fullfile(folderName, 'resultados.mat'));

fprintf("Resultados guardados en: %s\n", folderName);
