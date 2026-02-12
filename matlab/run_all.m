% run_all.m — single entry point

clear; close all; clc;

thisDir = fileparts(mfilename("fullpath"));
cd(thisDir);

addpath(thisDir);
addpath(fullfile(thisDir,"utils"));

run("config.m");   % creates CFG

fprintf("Running experiments (MATLAB R2024b)\n");

% Validations
run("OFDM_validacion.m");
run("OTFS_validacion.m");

% Sparse Doppler
run("OFDM_SNRvar.m");
run("OTFS_SNRvariable.m");

% 3GPP
run("OFDM_tdl.m");
run("OTFS_TDL.m");
run("OTFS_CDL.m");

fprintf("Done. Results in:\n%s\nFigures in:\n%s\n", CFG.out.results_dir, CFG.out.figures_dir);
