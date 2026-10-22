function [Cweave, aweave, Acol, Ccol, Atile] = WeaveHomog(Cb, ab, dx, dy, dz)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Purpose: Two-step GMC homogenisation of a weave RUC.
%          Bednarcyk, NASA/CR-2000-210370.
% Input:
% - Cb{b,g,k}: 6x6 stiffness of weave block (b,g,k) in weave axes
%              (Voigt order [11 22 33 23 13 12], engineering shear strains)
% - ab{b,g,k}: 6x1 CTE of the same block
% - dx, dy, dz: cell sizes along x1, x2, x3
% Output:
% - Cweave: 6x6 effective weave stiffness, in weave axes
% - aweave: 6x1 effective weave CTE
% - Acol{b,g}: step-A concentration tensors, as returned by GMCaniso
% - Ccol{b,g}: 6x6 effective column stiffness
% - Atile: step-B concentration tensors (stored in PERMUTED axes, see below)
%
% STEP A: each through-thickness column (NB = 1, NG = nz layers).
%         GMC axes coincide with weave axes, so no permutation is needed.
%         With NB = 1 the GMC equations reduce exactly to classical
%         laminate theory for a stack in x3.
% STEP B: the in-plane tile of columns, prismatic in x3. GMC treats its own
%         x1 as the prismatic direction, so weave x1 and x3 are swapped for
%         this call: GMC x1 <- weave x3, GMC x2 <- weave x2, GMC x3 <- weave x1.
%         The matching Voigt permutation pv = [3 2 1 6 5 4] is its own
%         inverse, so the same map converts in and back out.
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

pv = [3 2 1 6 5 4];

nx = numel(dx);
ny = numel(dy);
nz = numel(dz);

% -- Safety checks
if ~isequal(size(Cb), [nx ny nz]) || ~isequal(size(ab), [nx ny nz])
    error('WeaveHomog: Cb/ab must be %dx%dx%d cell arrays', nx, ny, nz);
end

%-----------------------------------------------------------------
% STEP A : homogenise each through-thickness column
%-----------------------------------------------------------------
Ccol  = cell(nx, ny);
acolC = cell(nx, ny);
Acol  = cell(nx, ny);

for b = 1:nx
    for g = 1:ny
        R = struct();
        R.NB = 1;
        R.NG = nz;
        R.h = dy(g);        % GMC x2 = weave x2, one cell wide
        R.l = dz;           % GMC x3 = weave x3, the layers
        for k = 1:nz
            R.C{1,k} = Cb{b,g,k};
            R.alpha{1,k} = ab{b,g,k};
        end
        [Ccol{b,g}, acolC{b,g}, A] = GMCaniso(R);
        Acol{b,g} = A;
    end
end

%-----------------------------------------------------------------
% STEP B : homogenise the in-plane tile of columns
%-----------------------------------------------------------------
R = struct();
R.NB = ny;                  % GMC x2 = weave x2  -> index g
R.NG = nx;                  % GMC x3 = weave x1  -> index b
R.h = dy;
R.l = dx;
for b = 1:nx
    for g = 1:ny
        R.C{g,b} = Ccol{b,g}(pv, pv);
        R.alpha{g,b} = acolC{b,g}(pv);
    end
end

[Ct, at, Atile] = GMCaniso(R);

% -- Undo the permutation
Cweave = Ct(pv, pv);
aweave = at(pv);

end
