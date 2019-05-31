%
% Simulation_Fig_S13
%
% Simulation of a Poincaré oscillator network as a model of the Choroid Plexus
% with increasing couplings
%
% Reference: Myung J, Schmal C et al. (2018) The Choroid Plexus is an Important
% Circadian Clock Component. Nat Commun, in press.
%
% Written by Sungho Hong, Computational Neuroscience, OIST, Japan
% Correspondence: Sungho Hong (shhong@oist.jp)
%
% February 9, 2018
%

clear;

nx = 48;
ny = 16;

cfg = config_CPnet_box([nx, ny], 4);

cfg.lambda = 0.025;
cfg.coupling_scale = 0.0005;
cfg.nhrs = 24*30*96/12;
cfg.dt = 1;
cfg.n_sample = 20;
cfg.scale_max = 5;
cfg.scale_min = cfg.scale_max/tstop;
cfg.A = 1;

f = load('init_data.mat');
cfg.omega    = f.omega;
cfg.init_phi = f.init_phi;
cfg.init_amp = f.init_amp;
cfg.epsilon = 0;

% Untwisted (1) and twisted (2) case
cfgs = [cfg, cfg];
cfg(1).epsilon = 0;
cfg(2).epsilon = -0.02;

nn = 2;

res = cell(nn,1);
qss = cell(nn,1);

for i=1:nn
    [z, qs] = simCPex_cont_double(cfgs(i));
    res{i} = z;
    qss{i} = qs;
end

save('fig_s13_data.mat', 'cfgs', 'res', 'qss')