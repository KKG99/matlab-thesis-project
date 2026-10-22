function [Hist, LamResults, plyprops, states] = WeaveLamSimDamage(Wp, Geometry, NxMax, NINC)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Purpose: Progressive-damage response of a laminate whose plies are woven.
%          CLT gives the ply strains, WeavePlyDamage degrades each ply at
%          the weave scale, then CLT is re-run. This replaces NASA's
%          LamSimDamage, which would damage at the tow scale instead.
% Input:
% - Wp: weave definition (Theory, towRUCid, RUCid, Vf, Vi, Nf, Nm, Mat,
%       DamageFactor optional)
% - Geometry: struct with Orient, plymat, tply (NASA LamProblemDef format)
% - NxMax: final applied Nx [N/mm]
% - NINC:  number of load increments
% Output:
% - Hist: .Nx .ex .ey .sigply (Nx2 per ply, material axes)
% - LamResults: last CLT result
% - plyprops: ply property cards at the final damaged state
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if ~isfield(Wp,'DamageFactor'), Wp.DamageFactor = 1e-4; end
if ~isfield(Wp,'NITER'), Wp.NITER = 3; end
if ~isfield(Wp,'Figures'), Wp.Figures = true; end
if ~isfield(Wp,'StrainControl'), Wp.StrainControl = false; end
if ~isfield(Wp,'epsxMax'), Wp.epsxMax = 0.006; end
for f = {'Mat','RUCid','towRUCid','Vf','Vi','Nf','Nm','Theory'}
    if ~isfield(Wp, f{1}), error(['WeaveLamSimDamage: Wp.' f{1} ' is required']); end
end

Nplies = numel(Geometry.Orient);

%-----------------------------------------------------------------
% One damage state per ply - plies are oriented differently, so they
% see different strains and must not share accumulated damage.
%-----------------------------------------------------------------
states = cell(1, Nplies);
for i = 1:Nplies
    states{i} = WeaveDamageInit(Wp);
end

%-----------------------------------------------------------------
% Undamaged ply card, in the shape NASA's Qbar expects
%-----------------------------------------------------------------
% -- One material card PER PLY. The 0 and 90 degree plies see different
%    strains, so they damage differently and cannot share a card.
%    Cards occupy Wp.Mat .. Wp.Mat+Nplies-1.
s0 = states{1};
E0 = GetEffPropsFromC(WeaveUndamagedC(s0));
plyprops = cell(1, 300);
for i = 1:Nplies
    m = Wp.Mat + i - 1;
    Geometry.plymat(i) = m;
    plyprops{m}.used = true;
    plyprops{m}.name = ['woven ply ' num2str(i)];
    plyprops{m}.E1 = E0.E1;  plyprops{m}.E2 = E0.E2;  plyprops{m}.E3 = E0.E3;
    plyprops{m}.v12 = E0.v12; plyprops{m}.v13 = E0.v13; plyprops{m}.v23 = E0.v23;
    plyprops{m}.G12 = E0.G12; plyprops{m}.G13 = E0.G13; plyprops{m}.G23 = E0.G23;
    plyprops{m}.a1 = 0;  plyprops{m}.a2 = 0;  plyprops{m}.a3 = 0;
end

% -- Current full 6x6 ply stiffness, one per ply, for plane-stress recovery
Cply = cell(1, Nplies);
Cund = WeaveUndamagedC(s0);
for i = 1:Nplies, Cply{i} = Cund; end

Hist.Nx    = zeros(NINC, 1);
Hist.ex    = zeros(NINC, 1);
Hist.ey    = zeros(NINC, 1);
Hist.sigply = zeros(NINC, 2*Nplies);
Hist.ndamply = zeros(NINC, Nplies);
Hist.nsubply = zeros(NINC, Nplies);
Eply = zeros(NINC, Nplies);
Sbply = cell(1, Nplies);
damfig = 0;

for INC = 1:NINC
    if Wp.StrainControl
        % -- Prescribe mid-plane strain in x (EK=1), leave the rest force-free (NM=2)
        LoadsI.Type  = [1, 2, 2, 2, 2, 2];
        LoadsI.Value = [Wp.epsxMax*INC/NINC, 0, 0, 0, 0, 0];
    else
        % -- Prescribe force resultant Nx (NM=2)
        LoadsI.Type  = [2, 2, 2, 2, 2, 2];
        LoadsI.Value = [NxMax*INC/NINC, 0, 0, 0, 0, 0];
    end
    LoadsI.DT = 0;
    LoadsI.NM = 2;                        % force/moment identifier
    LoadsI.EK = 1;                        % strain/curvature identifier

    for ITER = 1:Wp.NITER
        [Geometry, LamResults] = CLT(plyprops, Geometry, LoadsI);

        changed = false;
        for i = 1:Nplies
            % -- Ply strain in material axes, averaged over the thickness
            % -- CLT returns only the 3 in-plane components (MCstrain is
            %    3x2N), so recover the 3 out-of-plane ones from plane stress:
            %    sigma3 = sigma23 = sigma13 = 0  ->  eps_n = -Cnn\(Cnt*eps_t)
            t3 = [1 2 6];  n3 = [3 4 5];
            e3 = 0.5*(LamResults.MCstrain(:, 2*i-1) + LamResults.MCstrain(:, 2*i));
            epsply = zeros(6, 1);
            epsply(t3) = e3(:);
            Ccur = Cply{i};
            epsply(n3) = -Ccur(n3, n3) \ (Ccur(n3, t3) * epsply(t3));

            [C, ~, states{i}, Sb_i, ~, ndam_i, nsub_i] = WeavePlyDamage(states{i}, epsply, 2);
            Sbply{i} = Sb_i;
            Cply{i} = C;
            E = GetEffPropsFromC(C);
            m = Wp.Mat + i - 1;
            Hist.ndamply(INC, i) = ndam_i;
            Hist.nsubply(INC, i) = nsub_i;
            Eply(INC, i) = E.E1;

            if abs(E.E1 - plyprops{m}.E1) > 1e-9*max(E.E1, 1)
                changed = true;
            end
            plyprops{m}.E1 = E.E1;  plyprops{m}.E2 = E.E2;
            plyprops{m}.v12 = E.v12; plyprops{m}.G12 = E.G12;
            Hist.sigply(INC, 2*i-1) = LamResults.MCstress(1, 2*i-1);
            Hist.sigply(INC, 2*i)   = LamResults.MCstress(1, 2*i);
        end

        if ~changed && ITER > 1, break; end
    end

    if Wp.StrainControl
        Hist.Nx(INC) = LamResults.NM(1);      % reaction force for the imposed strain
    else
        Hist.Nx(INC) = NxMax*INC/NINC;
    end
    Hist.ex(INC) = LamResults.strain(1);
    Hist.ey(INC) = LamResults.strain(2);

    fprintf('  INC %3d/%3d  Nx = %10.3f  ex = %11.6e  E1 = %10.1f\n', ...
        INC, NINC, Hist.Nx(INC), Hist.ex(INC), plyprops{Wp.Mat}.E1);
    fprintf('               ply E1: 0deg %.1f   90deg %.1f\n', ...
        plyprops{Wp.Mat}.E1, plyprops{Wp.Mat+1}.E1);
    fprintf('               damaged blocks = %d   damaged tow subcells = %d\n', ...
        sum(Hist.ndamply(INC,:)), sum(Hist.nsubply(INC,:)));

    if Wp.Figures
        damfig = LamDamageWrite(damfig, Hist, Eply, INC, Nplies, states, Sbply, Geometry.Orient, sum(Geometry.tply));
    end
end

end

%=========================================================================%
% Live 3-panel figure for the laminate, mirroring the weave driver.
%   1) Nx vs ex            2) ply E1 vs increment   3) damaged blocks per ply
%=========================================================================%
function [damfig] = LamDamageWrite(damfig, Hist, Eply, INC, Nplies, states, Sbply, Orient, ttotal)
% Original 3 panels (stress-strain, ply stiffness, block damage) + a von Mises
% panel of the currently most-stressed ply, with x1/x2 axis directions.
% ttotal = total laminate thickness [mm], to convert Nx [N/mm] to avg stress [MPa].
if INC == 1
    damfig = figure('units','normalized','outerposition',[0 0.5 1 0.5]);
else
    figure(damfig);
end

subplot(2,2,1);
plot(Hist.ex(1:INC), Hist.Nx(1:INC)/ttotal, '-o', 'MarkerSize', 3);
xlabel('\epsilon_x'); ylabel('\sigma_{avg}  [MPa]');
title('stress vs. strain'); grid on; drawnow;

subplot(2,2,2);
plot(1:INC, Eply(1:INC, :)/1000, '-o', 'MarkerSize', 3);
xlabel('increment'); ylabel('ply E1  [GPa]');
title('ply stiffness vs. increment'); grid on; drawnow;

subplot(2,2,3);
bar([Hist.ndamply(INC, :), Hist.nsubply(INC, :)]);
xlabel('ply'); ylabel('count');
legend('damaged blocks', 'damaged tow subcells', 'Location', 'northwest');
title(['damage, INC = ',num2str(INC)]); drawnow;

[vm, kp] = vm_top(states, Sbply);
if ~isempty(vm)
    nx = states{kp}.nx; ny = states{kp}.ny;
    [X, Y] = meshgrid(0:nx, 0:ny);
    subplot(2,2,4);
    colormap(jet);
    pcolor(X, Y, pad2(vm.', ny, nx)); shading interp;
    colorbar; axis image;
    xlabel('x_1 (warp, load)'); ylabel('x_2 (weft)');
    title(['von Mises, ply ' num2str(kp) ' (' num2str(Orient(kp)) 'deg), INC = ' num2str(INC)]);
    drawnow;
end
end

%=========================================================================%
function [vm, kp] = vm_top(states, Sbply)
best = -1; kp = 1; vm = [];
for p = 1:numel(states)
    if isempty(Sbply{p}), continue; end
    v = vm_of(states{p}, Sbply{p});
    if max(v(:)) > best, best = max(v(:)); vm = v; kp = p; end
end
end

%=========================================================================%
function [v] = vm_of(st, Sb)
nx = st.nx; ny = st.ny; nz = st.nz; wruc = st.wruc;
v = zeros(nx, ny);
for b = 1:nx
    for g = 1:ny
        kk = nz;
        for k = nz:-1:1
            c = wruc.matsCh(b,g,k);
            if c == '5' || c == '6', kk = k; break; end
        end
        s = Sb{b,g,kk};
        v(b,g) = sqrt((s(1)-s(2))^2 + (s(2)-s(3))^2 + (s(1)-s(3))^2 + ...
                    6*s(4)^2 + 6*s(5)^2 + 6*s(6)^2)/sqrt(2);
    end
end
end

%=========================================================================%
function [P] = pad2(M, ny, nx)
P = zeros(ny+1, nx+1);
P(1:ny, 1:nx) = M;
P(ny+1, :) = P(ny, :);
P(:, nx+1) = P(:, nx);
end

%=========================================================================%
function [C] = WeaveUndamagedC(state)
[~, ~, ~, ~, C, ~] = WeaveFields(state.Cb, state.ab, state.wruc.d, ...
    state.wruc.h, state.wruc.l, zeros(6, 1));
end
