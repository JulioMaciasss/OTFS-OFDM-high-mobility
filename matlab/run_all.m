% run_all.m — one entry point for all experiments
clear; close all; clc;

thisDir = fileparts(mfilename("fullpath"));
cd(thisDir);

addpath(thisDir);
addpath(fullfile(thisDir,"utils"));

run("config.m");  % creates CFG

fprintf("Running experiments...\nResults -> %s\n", CFG.out.results_dir);

% 1) Validations (fast)
run("OFDM_validacion.m");
run("OTFS_validacion.m");

% 2) Sparse Doppler channel (custom dopplerChannel)
run("OFDM_SNRvar.m");
run("OTFS_SNRvariable.m");

% 3) 3GPP channels
run("OFDM_tdl.m");
run("OTFS_TDL.m");
run("OTFS_CDL.m");

fprintf("DONE.\n");
