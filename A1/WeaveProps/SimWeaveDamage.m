function [Hist, Weave, Tow, wruc] = SimWeaveDamage(Wp)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Purpose: Progressive-damage tensile response of a woven ply, driven and
%          drawn the same way as NASA's MicroSimDamage:
%            - ramp the load in NINC increments (default 100, like NASA),
%            - live 3-panel figure: stress-strain / von Mises / block failure,
%            - console progress "* INC = i out of N *" each step.
% Input (struct Wp) - everything RunWeave takes, plus:
% - Wp.epswMax : 6x1 final macro strain in weave axes
% - Wp.NINC    : load increments (default 100, same as NASA)
% - Wp.NITER   : max equilibrium iterations per increment (default 5)
% - Wp.DamageFactor : stiffness reduction, NASA default 1e-4
% - Wp.Figures : true for the live plots (default true here, like NASA)
% Output:
% - Hist: .eps .sig .ndam per increment
% - Weave, Tow, wruc: as from RunWeave, at the final state
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if ~isfield(Wp,'DamageFactor'), Wp.DamageFactor = 1e-4; end
if ~isfield(Wp,'NITER'),        Wp.NITER = 5;         end
if ~isfield(Wp,'NINC'),         Wp.NINC = 100;        end
if ~isfield(Wp,'epswMax'), error('SimWeaveDamage: Wp.epswMax is required'); end
if ~isfield(Wp,'Figures'), Wp.Figures = true; end

state = WeaveDamageInit(Wp);
nx = state.nx; ny = state.ny; nz = state.nz;
wruc = state.wruc;

% -- Totals for "damaged out of total" reporting
totBlk = nx*ny*nz;                                   % all weave blocks
towBlk = nnz(ismember(wruc.matsCh, '56'));           % blocks that are tows
totSub = towBlk * state.Tow.RUC.NB * state.Tow.RUC.NG;  % tow subcells in all tow blocks

N = Wp.NINC;
Hist.eps  = zeros(N+1, 6);
Hist.sig  = zeros(N+1, 6);
Hist.ndam = zeros(N+1, 1);
Hist.nsub = zeros(N+1, 1);
strainhist = zeros(N+1,1); stresshist = zeros(N+1,1);

damfig = 0;
ndamPrev = -1;

for INC = 0:N
    pct  = INC/N;
    epsw = Wp.epswMax * pct;

    if INC == 0
        % -- Undamaged baseline, no damage pass needed
        [Sb, Eb, Sc, ~, Cweave, aweave] = WeaveFields(state.Cb, state.ab, ...
            wruc.d, wruc.h, wruc.l, epsw);
        ndam = 0; nsub = 0;
    else
        [Cweave, aweave, state, Sb, Sc, ndam, nsub] = WeavePlyDamage(state, epsw, Wp.NITER);
    end

    sig = Cweave * epsw;
    ic = INC + 1;
    Hist.eps(ic, :) = epsw.';
    Hist.sig(ic, :) = sig.';
    Hist.ndam(ic) = ndam;
    Hist.nsub(ic) = nsub;
    strainhist(ic) = epsw(1);
    stresshist(ic) = sig(1);

    disp(' ');
    disp('*****************************************************');
    disp(['* INC = ',num2str(INC),'  out of ',num2str(N)]);
    disp('*****************************************************');
    fprintf('   eps11 = %9.6f  sig11 = %10.3f  blocks = %4d/%d  subcells = %5d/%d\n', ...
        epsw(1), sig(1), ndam, totBlk, nsub, totSub);

    if Wp.Figures
        newfail = (ndam ~= ndamPrev);
        damfig = WeaveDamageWrite(damfig, state, Sb, strainhist(1:ic), ...
                   stresshist(1:ic), INC, newfail, ndam, nsub);
        ndamPrev = ndam;
    end
end

Weave.RUCid = Wp.RUCid;
Weave.Cstar = Cweave;
Weave = AddEngProps(Weave, Cweave, aweave);
Tow = state.Tow;
wruc = state.wruc;

if Wp.Figures
    PlotWeaveBlocks3D(state, ...
        ['Weave RUC ' num2str(Wp.RUCid) ' - failed blocks, final']);
end

end

%=========================================================================%
% Live 3-panel figure, mirroring NASA's DamageWrite layout.
%=========================================================================%
function [damfig] = WeaveDamageWrite(damfig, state, Sb, strainhist, stresshist, INC, newfail, ndam, nsub)
nx = state.nx; ny = state.ny; nz = state.nz;
wruc = state.wruc;
totBlk = nx*ny*nz;
totSub = nnz(ismember(wruc.matsCh, '56')) * state.Tow.RUC.NB * state.Tow.RUC.NG;

if INC == 0
    damfig = figure('units','normalized','outerposition',[0 0.5 1 0.5]);
else
    figure(damfig);
end

% -- Panel 1: stress vs strain, always updated
subplot(1,3,1);
plot(strainhist, stresshist/1, '-o', 'MarkerSize', 3);
title('stress vs. strain'); xlabel('strain'); ylabel('stress (MPa)','rotation',90);
grid on; drawnow;

if INC == 0 || newfail
    % -- Build a top-view map: tow = 1 (blue), matrix = 2 (green), failed = 3 (red)
    map = 2*ones(nx, ny);
    vm  = zeros(nx, ny);
    for b = 1:nx
        for g = 1:ny
            kk = nz;
            for k = nz:-1:1
                c = wruc.matsCh(b,g,k);
                if c == '5' || c == '6', kk = k; map(b,g) = 1; break; end
            end
            s = Sb{b,g,kk};
            vm(b,g) = sqrt((s(1)-s(2))^2 + (s(2)-s(3))^2 + (s(1)-s(3))^2 + ...
                        6*s(4)^2 + 6*s(5)^2 + 6*s(6)^2)/sqrt(2);
            if wruc.matsCh(b,g,kk) == 'M'
                if state.mDam(b,g,kk), map(b,g) = 3; end
            else
                if any(state.mSt{b,g,kk}(:) == '123'), map(b,g) = 3; end
            end
        end
    end

    % -- Pad to (ny+1)x(nx+1) for pcolor
    Mp = padmap(map.', ny, nx);
    Vp = padmap(vm.',  ny, nx);
    [X, Y] = meshgrid(0:nx, 0:ny);

    % -- Panel 2: von Mises
    subplot(1,3,2);
    colormap(jet);
    pcolor(X, Y, Vp); shading interp;
    colorbar; axis image;
    xlabel('x_1 (warp, load)'); ylabel('x_2 (weft)');
    title(['von Mises stress, x_3 = thru-thk, INC = ',num2str(INC)],'FontSize',12);
    drawnow;

    % -- Panel 3: block failure
    subplot(1,3,3);
    colormap(jet);
    pcolor(X, Y, Mp); shading faceted;
    caxis([1 3]); axis image;
    xlabel('x_1 (warp, load)'); ylabel('x_2 (weft)');
    title(['Block failure (red), INC=',num2str(INC),':  blocks=',num2str(ndam),'/',num2str(totBlk),...
           '  subcells=',num2str(nsub),'/',num2str(totSub)],'FontSize',11);
    drawnow;
end
end

%=========================================================================%
function [P] = padmap(M, ny, nx)
P = zeros(ny+1, nx+1);
P(1:ny, 1:nx) = M;
P(ny+1, :) = P(ny, :);
P(:, nx+1) = P(:, nx);
end

%=========================================================================%
function [s] = AddEngProps(s, C, a)
[E] = GetEffPropsFromC(C);
s.E1 = E.E1;  s.E2 = E.E2;  s.E3 = E.E3;
s.v12 = E.v12; s.v13 = E.v13; s.v23 = E.v23;
s.G12 = E.G12; s.G13 = E.G13; s.G23 = E.G23;
s.a1 = a(1);   s.a2 = a(2);   s.a3 = a(3);
end
