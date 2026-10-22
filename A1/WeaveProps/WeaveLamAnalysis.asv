function [LamResults, plyprops, Geometry] = WeaveLamAnalysis(Wp, Geometry)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Purpose: Laminate analysis with woven plies - SAME LOGIC and SAME
%          RESULTS as NASA's LamAnalysis.m: ply card -> CLT -> OutputLam.
%          OutputLam writes the full NASA report (laminate properties from
%          the ABD matrix, ply stresses, plots) to Output/Composite Text
%          Report-...txt, exactly like LamAnalysis; this function also
%          prints A, B, D and the full 6x6 ABD to the command window and
%          returns them in LamResults.ABD.
%          PRESS RUN with no arguments: Wp from WeaveProblemDef.m, laminate
%          defaults to [0,90,90,0], tply 0.15 mm.
% Results (like LamAnalysis):
%   LamResults.A / B / D : 3x3 extensional / coupling / bending stiffness
%   LamResults.ABD       : full 6x6 [A B; B D]
%   LamResults.EX...AlphXY : laminate engineering properties from ABD
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if nargin < 1 || isempty(Wp)
    Wp = WeaveProblemDef();
end

addpath('Functions/CLT');
addpath('Functions/Utilities');
addpath('Functions/WriteResults');
addpath('Functions/Micromechanics');
addpath('Functions/Margins');
addpath('WeaveProps');

%-----------------------------------------------------------------
% 1) Problem definition (self-contained; override via Geometry arg)
%-----------------------------------------------------------------
if nargin < 2 || isempty(Geometry)
    Geometry.Orient = [0, 90, 90, 0];
    Geometry.tply   = 0.15*ones(1, 4);
end
if iscell(Geometry), Geometry = Geometry{1}; end
Nplies = numel(Geometry.Orient);
if ~isfield(Geometry, 'tply') || numel(Geometry.tply) ~= Nplies
    Geometry.tply = 0.15*ones(1, Nplies);
end
Geometry.plymat(1:Nplies) = Wp.Mat;

Loads.Type  = [2 2 2 2 2 2];
Loads.Value = [200, 0, 0, 0, 0, 0];
Loads.DT = 0;
Loads.NM = 2;                 % force/moment identifier (laminate problem)
Loads.EK = 1;                 % strain/curvature identifier

OutInfo.Format = "txt";
OutInfo.MakePlots = true;
OutInfo.WriteMargins = false;
OutInfo.Name(1) = string(['WEAVE laminate RUC ' num2str(Wp.RUCid)]);

%-----------------------------------------------------------------
% 2) Ply card: our woven ply from the two-step GMC at Wp.Mat
%-----------------------------------------------------------------
s0 = WeaveDamageInit(Wp);
[~, ~, ~, ~, Cund0, ~] = WeaveFields(s0.Cb, s0.ab, s0.wruc.d, ...
    s0.wruc.h, s0.wruc.l, zeros(6, 1));
E0 = GetEffPropsFromC(Cund0);
plyprops = cell(1, 300);
plyprops{Wp.Mat}.used = true;
plyprops{Wp.Mat}.name = ['woven ply, RUC ' num2str(Wp.RUCid)];
plyprops{Wp.Mat}.E1 = E0.E1;  plyprops{Wp.Mat}.E2 = E0.E2;  plyprops{Wp.Mat}.E3 = E0.E3;
plyprops{Wp.Mat}.v12 = E0.v12; plyprops{Wp.Mat}.v13 = E0.v13; plyprops{Wp.Mat}.v23 = E0.v23;
plyprops{Wp.Mat}.G12 = E0.G12; plyprops{Wp.Mat}.G13 = E0.G13; plyprops{Wp.Mat}.G23 = E0.G23;
plyprops{Wp.Mat}.a1 = 0;  plyprops{Wp.Mat}.a2 = 0;  plyprops{Wp.Mat}.a3 = 0;

disp(['Problem #1 - ', char(OutInfo.Name(1))]);
disp(['   Per ply angle orientations   [', num2str(Geometry.Orient), ']']);
disp(['   Per ply thicknesses          [', num2str(Geometry.tply), ']']);
disp(['   Per ply material assignments [', num2str(Geometry.plymat), ']']);
fprintf('woven ply (RUC %d):  E1 = %.1f  E2 = %.1f  G12 = %.1f MPa,  v12 = %.4f\n', ...
    Wp.RUCid, E0.E1, E0.E2, E0.G12, E0.v12);

%-----------------------------------------------------------------
% 3) CLT + OutputLam, exactly as in LamAnalysis
%-----------------------------------------------------------------
[Geometry, LamResults] = CLT(plyprops, Geometry, Loads);
[OutInfo] = OutputLam(OutInfo, 1, Geometry, Loads, plyprops, LamResults);

%-----------------------------------------------------------------
% 4) Results in the command window: ABD + laminate properties
%-----------------------------------------------------------------
LamResults.ABD = [LamResults.A, LamResults.B; LamResults.B, LamResults.D];

fprintf('\n================ ABD MATRIX (same as NASA LamAnalysis report) ================\n');
fprintf(' A - extensional stiffness [N/mm] (force per unit width per strain):\n');
pm3(LamResults.A);
fprintf(' B - coupling stiffness [N] (zero for symmetric laminates):\n');
pm3(LamResults.B);
fprintf(' D - bending stiffness [N*mm] (moment per curvature):\n');
pm3(LamResults.D);
fprintf(' full 6x6 ABD = [A B; B D]:\n');
pm6(LamResults.ABD);

fprintf('================ LAMINATE PROPERTIES (derived from ABD) ================\n');
fprintf(' EX    = %.6g MPa   Young''s modulus in x (load) direction\n', LamResults.EX);
fprintf(' EY    = %.6g MPa   Young''s modulus in y direction\n', LamResults.EY);
fprintf(' GXY   = %.6g MPa   in-plane shear modulus\n', LamResults.GXY);
fprintf(' NuXY  = %.6g       major Poisson''s ratio\n', LamResults.NuXY);
fprintf(' AlphX = %.6g 1/degC   CTE in x\n', LamResults.AlphX);
fprintf(' AlphY = %.6g 1/degC   CTE in y\n', LamResults.AlphY);
fprintf(' NASA report file: %s.txt\n', OutInfo.OutFile);

end

%=========================================================================%
function pm3(M)
for r = 1:3
    fprintf('  %13.5e %13.5e %13.5e\n', M(r,:));
end
end

%=========================================================================%
function pm6(M)
for r = 1:6
    fprintf('  %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n', M(r,:));
end
end
