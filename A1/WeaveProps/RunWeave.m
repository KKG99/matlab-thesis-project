function [Weave, Tow, Sb, Eb, Sc, Ec, wruc] = RunWeave(Wp)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Purpose: Effective properties of a woven ply.
%          Scale 1 (fibre+matrix -> tow) with NASA's own RunMicro.
%          Scale 2 (tow -> weave)   with the two-step GMC of WeaveHomog.
% Input (struct Wp):
% - Wp.Theory   : 'GMC' or 'HFGMC'   (applies at the tow scale)
% - Wp.towRUCid : tow RUC id          (e.g. 26)
% - Wp.RUCid    : weave RUC id        (401 plain, 402 twill, 403 8HS)
% - Wp.Vf, Wp.Vi: tow fibre / interface volume fractions
% - Wp.Nf, Wp.Nm: fibre and matrix material numbers in GetConstitProps
% Optional:
% - Wp.epsw : 6x1 macro strain in weave axes (default: unit strain along x1)
% - Wp.Figures : true to draw the weave stress plots
% - Wp.TitleStr: title string for those plots
% Output:
% - Weave: struct with RUCid, Vf (overall fibre fraction), Cstar, CTE and
%          the 12 engineering constants in weave axes
% - Tow:   struct with the same constants in tow axes, plus As and RUC
% - Sb{b,g,k}, Eb{b,g,k}: block stress and strain, 6x1, weave axes
% - Sc{b,g}, Ec{b,g}:     column stress and strain, 6x1, weave axes
% - wruc:  the weave RUC struct
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Safety checks
if ~isfield(Wp,'Theory') || ~ismember(Wp.Theory, {'GMC','HFGMC'})
    error('WeaveAnalysis: Wp.Theory must be ''GMC'' or ''HFGMC''');
end
for f = {'RUCid','towRUCid','Vf','Vi','Nf','Nm'}
    if ~isfield(Wp, f{1})
        error(['WeaveAnalysis: Wp.' f{1} ' is required']);
    end
end

%-----------------------------------------------------------------
% 1) Constituents, exactly as NASA's GetEffProps does it
%-----------------------------------------------------------------
[constitprops] = GetConstitProps;
Constits.Fiber  = constitprops{Wp.Nf};
Constits.Matrix = constitprops{Wp.Nm};

%-----------------------------------------------------------------
% 2) Scale 1 : fibre + matrix -> tow, with NASA's own RunMicro
%-----------------------------------------------------------------
Vol.Vf = Wp.Vf;
Vol.Vi = Wp.Vi;
Tow.used = true;
Tow = RunMicro(Wp.Theory, 'Tow', Vol, Constits, Tow, Wp.towRUCid);

Ctow = Tow.Cstar;
atow = [Tow.a1; Tow.a2; Tow.a3; 0; 0; 0];

% -- Matrix pocket, built the way GMC.m builds a constituent
Cmat = GetCFromProps(Constits.Matrix);
amat = [Constits.Matrix.aL; Constits.Matrix.aT; Constits.Matrix.aT; 0; 0; 0];

%-----------------------------------------------------------------
% 3) Weave RUC
%-----------------------------------------------------------------
[wruc] = GetRUC(Constits, Wp.Vf, Wp.Vi, Wp.RUCid);

nx = numel(wruc.d);
ny = numel(wruc.h);
nz = numel(wruc.l);

if ~isequal(size(wruc.matsCh), [nx ny nz])
    error('WeaveAnalysis: RUC.matsCh must be %dx%dx%d', nx, ny, nz);
end

%-----------------------------------------------------------------
% 4) Blocks : every weave cell is the tow rotated into place, or matrix
%-----------------------------------------------------------------
Cb = cell(nx, ny, nz);
ab = cell(nx, ny, nz);
for b = 1:nx
    for g = 1:ny
        for k = 1:nz
            switch wruc.matsCh(b, g, k)
                case '5'      % warp tow, axis along x1, crimp in x1-x3
                    [Cb{b,g,k}, ab{b,g,k}] = TowCfull(Ctow, atow, 0,  wruc.crimp(b,g,k));
                case '6'      % weft tow, axis along x2, crimp in x2-x3
                    [Cb{b,g,k}, ab{b,g,k}] = TowCfull(Ctow, atow, 90, wruc.crimp(b,g,k));
                case 'M'      % matrix pocket
                    Cb{b,g,k} = Cmat;
                    ab{b,g,k} = amat;
                otherwise
                    error(['WeaveAnalysis: bad matsCh = ' wruc.matsCh(b,g,k)])
            end
        end
    end
end

% -- Overall fibre volume fraction, from the RUC geometry
towVol = 0;
for b = 1:nx
    for g = 1:ny
        for k = 1:nz
            if wruc.matsCh(b,g,k) ~= 'M'
                towVol = towVol + wruc.d(b)*wruc.h(g)*wruc.l(k);
            end
        end
    end
end
VfWeave = Wp.Vf * towVol / (sum(wruc.d)*sum(wruc.h)*sum(wruc.l));

%-----------------------------------------------------------------
% 5) Scale 2 : tow -> weave, two-step GMC, plus the block fields
%-----------------------------------------------------------------
if isfield(Wp, 'epsw')
    epsw = Wp.epsw;
else
    epsw = [1; 0; 0; 0; 0; 0];      % default: unit strain along the warp
end
[Sb, Eb, Sc, Ec, Cweave, aweave] = WeaveFields(Cb, ab, wruc.d, wruc.h, wruc.l, epsw);

% -- Stress distribution on the weave RUC grid
if isfield(Wp, 'Figures') && Wp.Figures
    if isfield(Wp, 'TitleStr')
        PlotWeaveFields(wruc, Sc, Sb, Wp.TitleStr);
    else
        PlotWeaveFields(wruc, Sc, Sb, ['weave RUC ' num2str(Wp.RUCid)]);
    end
end

Weave.RUCid = Wp.RUCid;
Weave.Vf    = VfWeave;
Weave.Cstar = Cweave;
Weave.CTE   = aweave;
[Weave] = AddEngProps(Weave, Cweave, aweave);

fprintf('\n            E1        E2        E3      v12      v13      v23      G12       G13       G23      CTE1     CTE2     CTE3\n');
fprintf('tow   %9.1f %9.1f %9.1f %8.4f %8.4f %8.4f %9.1f %9.1f %9.1f %8.2e %8.2e %8.2e\n', ...
    Tow.E1, Tow.E2, Tow.E3, Tow.v12, Tow.v13, Tow.v23, Tow.G12, Tow.G13, Tow.G23, ...
    Tow.a1, Tow.a2, Tow.a3);
fprintf('weave %9.1f %9.1f %9.1f %8.4f %8.4f %8.4f %9.1f %9.1f %9.1f %8.2e %8.2e %8.2e\n', ...
    Weave.E1, Weave.E2, Weave.E3, Weave.v12, Weave.v13, Weave.v23, ...
    Weave.G12, Weave.G13, Weave.G23, Weave.a1, Weave.a2, Weave.a3);
fprintf('\nweave RUCid = %d   overall fibre Vf = %.6f   (tow Vf = %.4f)\n\n', ...
    Weave.RUCid, Weave.Vf, Wp.Vf);

end

%=========================================================================
function [s] = AddEngProps(s, C, a)
[E] = GetEffPropsFromC(C);
s.E1  = E.E1;   s.E2  = E.E2;   s.E3  = E.E3;
s.v12 = E.v12;  s.v13 = E.v13;  s.v23 = E.v23;
s.G12 = E.G12;  s.G13 = E.G13;  s.G23 = E.G23;
s.a1 = a(1);    s.a2 = a(2);    s.a3 = a(3);
end
