function [Wp] = WeaveProblemDef()
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Purpose: Define the weave problem, exactly as MicroProblemDef does for
%          NASA's tow problems. EDIT THIS FILE, then Run WeaveAnalysis.m
%          or WeaveSimDamage.m. You never need to run this file itself.
% Output:
% - Wp: struct of weave problem settings
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Micromechanics theory at the tow scale
Wp.Theory   = 'GMC';          % 'GMC' or 'HFGMC'

% -- Tow scale
Wp.towRUCid = 7;             % 26 for properties, 7 for damage runs (faster)
Wp.Vf       = 0.490;           % fibre volume fraction in the tow
Wp.Vi       = 0;              % interface volume fraction (0 = none)
Wp.Nf       = 19;             % fibre material number in GetConstitProps
Wp.Nm       = 20;             % matrix material number in GetConstitProps

% -- Weave scale: 401 plain, 402 twill 2/2, 403 8-harness satin
Wp.RUCid    = 401;

% -- Damage settings (used by WeaveSimDamage only)
Wp.NINC         = 3;          % number of load increments
Wp.NITER        = 4;          % max iterations per increment
Wp.DamageFactor = 1e-4;       % NASA default, do not change
Wp.epswMax      = [0.006; 0; 0; 0; 0; 0];   % final strain, along the warp
Wp.Figures      = true;      % true to draw the weave stress fields

% -- Laminate settings (used by WeaveLamDamage only)
Wp.Mat   = 220;               % first ply material number; uses Mat..Mat+3
Wp.NxMax = 500;                % final applied Nx [N/mm]
Wp.StrainControl = false;     % pull by strain, like your displacement test
Wp.epsxMax       = 0.009;    % final mid-plane strain
Wp.Orient = [0 90 0 90];      % ply angles
Wp.tply   = [0.15 0.15 0.15 0.15];

end
