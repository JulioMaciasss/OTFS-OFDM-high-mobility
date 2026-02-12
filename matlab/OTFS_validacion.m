% OTFS_validacion.m — OTFS validation (AWGN + optional Doppler channel)
% Reproducible version (uses CFG when available), fixed BLER logic

clearvars -except CFG;
close all; clc;

%% Paths
thisDir = fileparts(mfilename("fullpath"));
addpath(fullfile(thisDir,"utils"));

%% ---- Load central config if available ----
if exist("CFG","var")
    M = CFG.M; N = CFG.N; df = CFG.df; fc = CFG.fc; %#ok<NASGU>
    padLen  = CFG.padLen_default;
    padType = char(CFG.padType);

    SNRdBvalues = CFG.SNRdBvalues;
    numbslots  = CFG.numbslots;

    basePath = CFG.out.results_dir;

    % Validation is AWGN by default
    channel_type = "awgn";
    useDoppler = false;

    % Doppler params (only used if useDoppler==true)
    chanParams.pathDelays   = CFG.doppler.pathDelays;
    chanParams.pathGains    = CFG.doppler.pathGains;
    chanParams.pathDopplers = CFG.doppler.pathDopplers;

    threshold = CFG.beta; %#ok<NASGU>  (only used if Doppler enabled)
else
    M = 64; N = 30; df = 15e3; fc = 3e9; %#ok<NASGU>
    padLen = 10; padType = 'ZP';

    SNRdBvalues = linspace(0,20,12);
    numbslots  = 100;

    basePath = fullfile(pwd,"results");

    channel_type = "awgn";
    useDoppler = false;

    chanParams.pathDelays   = [0 5 8];
    chanParams.pathGains    = [1 0.7 0.5];
    chanParams.pathDopplers = [0 -3 5];

    threshold = 0.05; %#ok<NASGU>
end
% ------------------------------------------

%% Derived
fsamp = M*df;
Meff = M + padLen;
numSamps = Meff * N;
T = (M+padLen)/(M*df); %#ok<NASGU>

chanParams.pathDopplerFreqs = chanParams.pathDopplers * (1/(N*T));

%% Output folder
folderName = fullfile(basePath, sprintf("OTFS_validation_%s_padLen%d", channel_type, padLen));
if ~exist(folderName,'dir'); mkdir(folderName); end

%% Pilot in delay–Doppler (single pilot)
pilotBin = floor(N/2)+1;
Pdd = zeros(M,N);
Pdd(1,pilotBi
