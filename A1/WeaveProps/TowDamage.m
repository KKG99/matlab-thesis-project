function [Cdam, adam, nfail, matsChDam] = TowDamage(tow, epsTow, DamageFactor, matsCh0)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Purpose: Degrade a tow whose constituents have failed, using NASA's OWN
%          damage law (the RUC branch of DamageMicro.m):
%            * a subcell fails only when its MoS < 0,
%            * and only the subcells within tol of the MOST negative margin
%              fail on a given call (worst-first),
%            * already-failed subcells stay failed (irreversible),
%            * the tow is re-homogenised with GMC only if something failed.
%          The outer iteration loop (WeavePlyDamage) re-solves and calls this
%          again, so damage accumulates gradually exactly like NASA's tow.
% Input:
% - tow: struct returned by RunMicro (Fiber, Matrix, RUC, Cstar, a1..a3, micro)
% - epsTow: 6x1 strain in TOW axes (engineering shear strains)
% - DamageFactor: stiffness reduction for a failed constituent (NASA 1e-4)
% - matsCh0: current damaged state of RUC.matsCh ('1'/'2'/'3' stay damaged)
% Output:
% - Cdam, adam: degraded tow stiffness / CTE in tow axes
% - nfail: number of degraded subcells (incl. previously degraded)
% - matsChDam: resulting matsCh, pass back in next time
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Constits.Fiber     = tow.Fiber;
Constits.Matrix    = tow.Matrix;
Constits.FiberDam  = DamageConstit(tow.Fiber,  DamageFactor);
Constits.MatrixDam = DamageConstit(tow.Matrix, DamageFactor);

% -- Subcell fields, NASA's own function
Res = MicroFields(epsTow, 0, tow, struct());

% -- Margins: Tsai-Wu only, NASA's default
CritOn = [0 0 0 1];
z6 = zeros(6, 1);

if nargin < 4 || isempty(matsCh0)
    matsCh0 = tow.RUC.matsCh;
end
RUCd = tow.RUC;
RUCd.matsCh = matsCh0;

NB = tow.RUC.NB; NG = tow.RUC.NG;

% -- MoS of every not-yet-failed subcell
MoSgrid = zeros(NB, NG);
failed0 = false(NB, NG);
for b = 1:NB
    for g = 1:NG
        if any(matsCh0(b, g) == '123')
            failed0(b, g) = true;
            continue;
        end
        if tow.RUC.matsCh(b, g) == 'F'
            allow = tow.Fiber.allowables;
        else
            allow = tow.Matrix.allowables;
        end
        M = FCriteria(Res.MicroFields(b,g).stress, Res.MicroFields(b,g).strain, ...
                      z6, z6, allow, CritOn, true, struct());
        MoSgrid(b, g) = M.MinMoS;
    end
end

% -- min MoS over unfailed subcells (NASA: starts at 0)
minMoS = 0;
for b = 1:NB
    for g = 1:NG
        if ~failed0(b, g), minMoS = min(minMoS, MoSgrid(b, g)); end
    end
end

% -- NASA tolerance: only fail subcells this close to the minimum
tol = 0.001;

nfail = 0;
for b = 1:NB
    for g = 1:NG
        if failed0(b, g)
            nfail = nfail + 1;
            continue;
        end
        if (minMoS < 0) && (MoSgrid(b, g) - minMoS - tol < 0)
            if tow.RUC.matsCh(b, g) == 'F'
                RUCd.matsCh(b, g) = '1';
            else
                RUCd.matsCh(b, g) = '2';
            end
            nfail = nfail + 1;
        end
    end
end

matsChDam = RUCd.matsCh;

if nfail == 0
    Cdam = tow.Cstar;
    adam = [tow.a1; tow.a2; tow.a3; 0; 0; 0];
else
    [Cdam, CTEs] = GMC(Constits, RUCd);
    adam = CTEs;
end

end
