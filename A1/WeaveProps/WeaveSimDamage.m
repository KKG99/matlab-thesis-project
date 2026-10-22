% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% PRESS RUN ON THIS FILE. Weave progressive damage, stress-strain curve.
% Settings are in WeaveProblemDef.m - edit that, not this.
% This takes about a minute. It is not frozen.
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear; clc;
addpath('Functions/Utilities');
addpath('Functions/WriteResults');
addpath('Functions/Micromechanics');
addpath('Functions/Margins');
addpath('Functions/Damage');
addpath('WeaveProps');

[Wp] = WeaveProblemDef();
Wp.towRUCid = 7;        % damage runs use the faster tow RUCtthe s

disp(['Weave Damage - RUC ' num2str(Wp.RUCid) ', ' num2str(Wp.NINC) ' increments']);

[Hist, Weave, Tow, wruc] = SimWeaveDamage(Wp);

[pk, ip] = max(Hist.sig(:,1));
fprintf('\npeak sig11 = %.2f MPa at eps11 = %.6f\n', pk, Hist.eps(ip,1));

figure;
plot(Hist.eps(:,1), Hist.sig(:,1), 'b-o', 'MarkerSize', 4);
xlabel('\epsilon_{11}'); ylabel('\sigma_{11}  [MPa]');
title(['Weave RUC ' num2str(Wp.RUCid) ' - progressive damage']);
grid on;
