% config.m — central editable parameters for the whole repo (MATLAB R2024b)

%% Core simulation parameters (defaults from your scripts)
CFG.M = 64;
CFG.N = 30;
CFG.df = 15e3;
CFG.fc = 3e9;

CFG.padType = "ZP";

CFG.SNRdBvalues = linspace(0,20,12);
CFG.numbslots  = 100;

CFG.modorder = 4;        % QPSK
CFG.idealCE  = false;
CFG.imagesturn = false;
CFG.channelresponse = false;

%% Padding choices per experiment
CFG.padLen_default = 10; % validation + sparse Doppler
CFG.padLen_tdl     = 10; % you can set to 30 if DS is large
CFG.padLen_cdl     = 30; % your OTFS_CDL uses 30 [file:256]

%% Sparse delay–Doppler channel (indices)
CFG.doppler.pathDelays   = [0 5 8];
CFG.doppler.pathGains    = [1 0.7 0.5];
CFG.doppler.pathDopplers = [-2];  % matches your OTFS_SNRvariable example [file:258]

%% OTFS detector beta/threshold
CFG.beta = 0.05;   % common default; CDL script used 0.0305 [file:256]

%% 3GPP TDL / CDL settings (avoid input())
CFG.tdl.DelayProfile = "TDL-A";
CFG.tdl.MaximumDopplerShift = 300;   % Hz (e.g., 3, 30, 300)
CFG.tdl.DelaySpread = 300e-9;        % seconds (30e-9 / 300e-9 / 3000e-9)

CFG.cdl.DelayProfile = "CDL-A";
CFG.cdl.MaximumDopplerShift = 300;   % Hz
CFG.cdl.DelaySpread = 300e-9;        % seconds

%% Output (store INSIDE repo, not C:\...)
CFG.out.results_dir = fullfile(pwd, "..", "results");
CFG.out.figures_dir = fullfile(pwd, "..", "docs", "figures");

if ~exist(CFG.out.results_dir,"dir"); mkdir(CFG.out.results_dir); end
if ~exist(CFG.out.figures_dir,"dir"); mkdir(CFG.out.figures_dir); end
