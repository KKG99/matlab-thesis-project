% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% PRESS RUN ON THIS FILE. Weave properties: fibre+matrix -> tow -> weave.
% Settings are in WeaveProblemDef.m - edit that, not this.
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear; clc;
addpath('Functions/Utilities');
addpath('Functions/WriteResults');
addpath('Functions/Micromechanics');
addpath('Functions/Margins');
addpath('Functions/Damage');
addpath('WeaveProps');

[Wp] = WeaveProblemDef();

disp(['Weave Problem - RUC ' num2str(Wp.RUCid) ', ' Wp.Theory ...
      ', tow RUC ' num2str(Wp.towRUCid) ', Vf = ' num2str(Wp.Vf)]);

[Weave, Tow] = RunWeave(Wp);

fprintf('\nweave E1 = %.1f   E2 = %.1f   E3 = %.1f   G12 = %.1f MPa\n', ...
    Weave.E1, Weave.E2, Weave.E3, Weave.G12);
fprintf('overall fibre Vf = %.6f\n', Weave.Vf);
