function [Cweave, aweave, state, Sb, Sc, ndam, nsub] = WeavePlyDamage(state, epsw, NITER)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Purpose: Damaged weave stiffness at one given strain, using NASA's damage
%          law at BOTH scales:
%            - tow blocks  -> TowDamage (NASA DamageMicro law, worst-first),
%            - matrix blk  -> same law: fail only blocks within tol of the
%                             most negative margin, MoS < 0.
%          Damage is permanent (state carries mSt and mDam forward), and the
%          iteration loop re-solves until no new failure, so damage builds up
%          gradually with load, exactly like NASA's tow damage.
% Input:  state (WeaveDamageInit or previous call), epsw 6x1, NITER
% Output: Cweave, aweave, state, Sb, Sc, ndam
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if nargin < 3 || isempty(NITER), NITER = 4; end
DF  = state.Wp.DamageFactor;
tol = 0.001;                       % NASA DamageMicro tolerance
wruc = state.wruc;
nx = state.nx;  ny = state.ny;  nz = state.nz;
CritOn = [0 0 0 1];
z6 = zeros(6, 1);

for ITER = 1:NITER
    [Sb, Eb, Sc, ~, Cweave, aweave] = ...
        WeaveFields(state.Cb, state.ab, wruc.d, wruc.h, wruc.l, epsw);

    % -- min MoS over the NOT-yet-damaged matrix blocks (NASA starts at 0)
    minMoSm = 0;
    for b = 1:nx
        for g = 1:ny
            for k = 1:nz
                if wruc.matsCh(b,g,k) == 'M' && ~state.mDam(b,g,k)
                    M = FCriteria(Sb{b,g,k}, Eb{b,g,k}, z6, z6, ...
                        state.Constits.Matrix.allowables, CritOn, true, struct());
                    minMoSm = min(minMoSm, M.MinMoS);
                end
            end
        end
    end

    newdam = 0;
    for b = 1:nx
        for g = 1:ny
            for k = 1:nz
                switch wruc.matsCh(b, g, k)
                    case {'5','6'}
                        epsTow = state.Teps{b,g,k} * Eb{b,g,k};
                        [Cd, ad, ~, mcNew] = TowDamage(state.Tow, epsTow, DF, state.mSt{b,g,k});
                        if ~isequal(mcNew, state.mSt{b,g,k}), newdam = newdam + 1; end
                        state.mSt{b,g,k} = mcNew;
                        % -- FIX: TowDamage returns the degraded tow C in TOW axes;
                        %    rotate back to weave axes with the same orientation
                        %    TowCfull applied at init ('5' -> 0deg, '6' -> 90deg, + crimp)
                        if wruc.matsCh(b,g,k) == '5', thz = 0; else, thz = 90; end
                        [Cw2, aw2] = TowCfull(Cd, ad, thz, wruc.crimp(b,g,k));
                        state.Cb{b,g,k} = Cw2;  state.ab{b,g,k} = aw2;
                    case 'M'
                        if ~state.mDam(b,g,k)
                            M = FCriteria(Sb{b,g,k}, Eb{b,g,k}, z6, z6, ...
                                state.Constits.Matrix.allowables, CritOn, true, struct());
                            if (minMoSm < 0) && (M.MinMoS - minMoSm - tol < 0)
                                state.mDam(b,g,k) = true;
                                newdam = newdam + 1;
                            end
                        end
                        if state.mDam(b,g,k), state.Cb{b,g,k} = state.CmatD;
                        else,                 state.Cb{b,g,k} = state.Cmat0; end
                        state.ab{b,g,k} = state.amat;
                end
            end
        end
    end
    if newdam == 0 && ITER > 1, break; end
end

% -- Converged pass so the returned stiffness matches the final damage state
[Sb, Eb, Sc, ~, Cweave, aweave] = ...
    WeaveFields(state.Cb, state.ab, wruc.d, wruc.h, wruc.l, epsw);

ndam = 0; nsub = 0;
for b = 1:nx
    for g = 1:ny
        for k = 1:nz
            if wruc.matsCh(b,g,k) == 'M'
                if state.mDam(b,g,k), ndam = ndam + 1; end
            else
                d = sum(ismember(state.mSt{b,g,k}(:), '123'));  % damaged tow subcells
                nsub = nsub + d;
                if d > 0, ndam = ndam + 1; end                   % block counts as damaged
            end
        end
    end
end

end
