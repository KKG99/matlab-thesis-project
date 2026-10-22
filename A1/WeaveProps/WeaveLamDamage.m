% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% PRESS RUN ON THIS FILE. Laminate progressive damage, weave plies, CLT.
% Settings are in WeaveProblemDef.m - edit that, not this.
% Slow: about a minute per load increment.
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear; clc;
addpath('Functions/Utilities');
addpath('Functions/WriteResults');
addpath('Functions/Micromechanics');
addpath('Functions/Margins');
addpath('Functions/Damage');
addpath('Functions/CLT');
addpath('WeaveProps');

[Wp] = WeaveProblemDef();
Wp.towRUCid = 7;

Geometry.Orient = Wp.Orient;
Geometry.plymat = repmat(Wp.Mat, 1, numel(Wp.Orient));
Geometry.tply   = Wp.tply;

disp(['Laminate Damage - [' num2str(Wp.Orient) '], weave RUC ' num2str(Wp.RUCid)]);

[Hist, LamResults, plyprops, states] = WeaveLamSimDamage(Wp, Geometry, Wp.NxMax, 20);

figure;
plot(Hist.ex, Hist.Nx, 'b-o', 'MarkerSize', 4);
xlabel('\epsilon_x'); ylabel('N_x  [N/mm]');
title('Cross-ply laminate, woven plies');
grid on;

% -- 3-D failed-block picture of every ply, like NASA's Subcell Failure panel
for i = 1:numel(states)
    PlotWeaveBlocks3D(states{i}, ...
        ['ply ' num2str(i) ' (' num2str(Wp.Orient(i)) 'deg) - failed blocks, final']);
end
