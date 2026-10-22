function [state] = WeaveDamageInit(Wp)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Purpose: Build the undamaged weave damage state for WeavePlyDamage.
%          Call once, then reuse for every strain the laminate driver asks
%          about, so the tow homogenisation is not repeated.
% Input (struct Wp): Theory, towRUCid, RUCid, Vf, Vi, Nf, Nm
%                    DamageFactor optional (NASA default 1e-4)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if ~isfield(Wp,'DamageFactor'), Wp.DamageFactor = 1e-4; end
for f = {'RUCid','towRUCid','Vf','Vi','Nf','Nm','Theory'}
    if ~isfield(Wp, f{1}), error(['WeaveDamageInit: Wp.' f{1} ' is required']); end
end

[cp] = GetConstitProps;
Constits.Fiber  = cp{Wp.Nf};
Constits.Matrix = cp{Wp.Nm};

Vol.Vf = Wp.Vf;  Vol.Vi = Wp.Vi;
Tow.used = true;
Tow = RunMicro(Wp.Theory, 'Tow', Vol, Constits, Tow, Wp.towRUCid);
atow = [Tow.a1; Tow.a2; Tow.a3; 0; 0; 0];

[wruc] = GetRUC(Constits, Wp.Vf, Wp.Vi, Wp.RUCid);
nx = numel(wruc.d);  ny = numel(wruc.h);  nz = numel(wruc.l);

state.Wp = Wp;  state.Tow = Tow;  state.Constits = Constits;  state.wruc = wruc;
state.nx = nx;  state.ny = ny;  state.nz = nz;
state.Cmat0 = GetCFromProps(Constits.Matrix);
state.amat  = [Constits.Matrix.aL; Constits.Matrix.aT; Constits.Matrix.aT; 0; 0; 0];
state.CmatD = GetCFromProps(DamageConstit(Constits.Matrix, Wp.DamageFactor));
state.atow = atow;

state.Cb   = cell(nx, ny, nz);
state.ab   = cell(nx, ny, nz);
state.Teps = cell(nx, ny, nz);
state.mSt  = cell(nx, ny, nz);
state.mDam = false(nx, ny, nz);

for b = 1:nx
    for g = 1:ny
        for k = 1:nz
            switch wruc.matsCh(b, g, k)
                case '5'
                    [state.Cb{b,g,k}, state.ab{b,g,k}, state.Teps{b,g,k}] = ...
                        TowCfull(Tow.Cstar, atow, 0,  wruc.crimp(b,g,k));
                    state.mSt{b,g,k} = Tow.RUC.matsCh;
                case '6'
                    [state.Cb{b,g,k}, state.ab{b,g,k}, state.Teps{b,g,k}] = ...
                        TowCfull(Tow.Cstar, atow, 90, wruc.crimp(b,g,k));
                    state.mSt{b,g,k} = Tow.RUC.matsCh;
                case 'M'
                    state.Cb{b,g,k} = state.Cmat0;
                    state.ab{b,g,k} = state.amat;
                    state.Teps{b,g,k} = eye(6);
                    state.mSt{b,g,k} = [];
                otherwise
                    error(['WeaveDamageInit: bad matsCh = ' wruc.matsCh(b,g,k)])
            end
        end
    end
end

end
