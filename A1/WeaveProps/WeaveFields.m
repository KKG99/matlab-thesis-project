function [Sb, Eb, Sc, Ec, Cweave, aweave] = WeaveFields(Cb, ab, dx, dy, dz, epsw)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Purpose: Strain and stress in every block of a weave RUC under a macro
%          strain, using the concentration tensors of the two-step GMC.
% Input:
% - Cb{b,g,k}: 6x6 block stiffness in weave axes
% - ab{b,g,k}: 6x1 block CTE
% - dx, dy, dz: cell sizes along x1, x2, x3
% - epsw: 6x1 macro strain in weave axes (engineering shear strains)
% Output:
% - Sb{b,g,k}: 6x1 block stress in weave axes
% - Eb{b,g,k}: 6x1 block strain in weave axes
% - Sc{b,g}:   6x1 stress of the whole through-thickness column (b,g)
% - Ec{b,g}:   6x1 strain of that column
% - Cweave, aweave: effective weave stiffness and CTE
%
% Chain: macro strain -> (step B) column strain -> (step A) block strain.
% Step B works in permuted axes (GMC x1 <- weave x3), so its concentration
% tensor is applied in permuted axes and the result mapped back with
% pv = [3 2 1 6 5 4]. Step A needs no permutation.
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

pv = [3 2 1 6 5 4];

nx = numel(dx);
ny = numel(dy);
nz = numel(dz);

if ~isequal(size(Cb), [nx ny nz])
    error('WeaveFields: Cb must be %dx%dx%d', nx, ny, nz);
end
if numel(epsw) ~= 6
    error('WeaveFields: epsw must be 6x1');
end

% -- Homogenise, keeping the concentration tensors
[Cweave, aweave, Acol, Ccol, Atile] = WeaveHomog(Cb, ab, dx, dy, dz);

% -- Macro strain in step-B (permuted) axes
epswp = epsw(pv);

Sb = cell(nx, ny, nz);
Eb = cell(nx, ny, nz);
Sc = cell(nx, ny);
Ec = cell(nx, ny);
for b = 1:nx
    for g = 1:ny
        % -- Step B : macro -> column. Atile is indexed (g, b), see WeaveHomog
        epscol = Atile(g, b).A * epswp;
        epscol = epscol(pv);              % back to weave axes

        Ec{b,g} = epscol;
        Sc{b,g} = Ccol{b,g} * epscol;

        % -- Step A : column -> blocks. No permutation
        for k = 1:nz
            Eb{b,g,k} = Acol{b,g}(1, k).A * epscol;
            Sb{b,g,k} = Cb{b,g,k} * Eb{b,g,k};
        end
    end
end

end
