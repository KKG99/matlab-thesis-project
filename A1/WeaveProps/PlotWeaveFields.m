function hfig = PlotWeaveFields(wruc, Sc, Sb, TitleStr, hfig)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Purpose: Pseudocolor plots of the weave stress field on the weave RUC grid.
%          No fibres are shown - the tow is the smallest unit here.
%          Figure 1: step-B column stresses, all 6 components, on the tile.
%          Figures 2-7: step-A block stresses, one figure per component,
%                       one panel per layer through the thickness.
% Input:
% - wruc: weave RUC struct (needs d, h, l)
% - Sc{b,g}:   6x1 column stress
% - Sb{b,g,k}: 6x1 block stress
% - TitleStr: string for the figure titles
% - hfig: optional 7x1 vector of existing figure handles to redraw into.
%         Pass the output back in to update the same figures in place,
%         which is what the damage loop needs. Omit to create new ones.
% Output:
% - hfig: the 7x1 vector of figure handles used
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if nargin < 5 || isempty(hfig)
    hfig = zeros(7, 1);
end

nx = numel(wruc.d);
ny = numel(wruc.h);
nz = numel(wruc.l);

comp = {'\sigma_{11}','\sigma_{22}','\sigma_{33}', ...
        '\sigma_{23}','\sigma_{13}','\sigma_{12}'};

% -- Cell edge coordinates (cells are NOT all the same width)
xe = [0 cumsum(wruc.d)];
ye = [0 cumsum(wruc.h)];

%-----------------------------------------------------------------
% Figure 1 : step-B column stresses on the tile
%-----------------------------------------------------------------
if hfig(1) == 0, hfig(1) = figure('Name','Step B column stresses'); else, figure(hfig(1)); clf; end
for c = 1:6
    M = zeros(nx, ny);
    for b = 1:nx
        for g = 1:ny
            M(b, g) = Sc{b, g}(c);
        end
    end
    subplot(2, 3, c);
    TilePlot(xe, ye, M, [comp{c} '  [MPa]']);
end
sgtitle(['Step B - column stresses, ' TitleStr]);

%-----------------------------------------------------------------
% Figures 2-7 : step-A block stresses, one per component, 4 layers
%-----------------------------------------------------------------
for c = 1:6
    if hfig(c+1) == 0, hfig(c+1) = figure('Name',['Step A ' comp{c}]); else, figure(hfig(c+1)); clf; end
    for k = 1:nz
        M = zeros(nx, ny);
        for b = 1:nx
            for g = 1:ny
                M(b, g) = Sb{b, g, k}(c);
            end
        end
        subplot(2, 2, k);
        TilePlot(xe, ye, M, ['layer ' num2str(k)]);
    end
    sgtitle(['Step A - ' comp{c} ' by layer [MPa], ' TitleStr]);
end

end

%=========================================================================
function TilePlot(xe, ye, M, ttl)
% Pseudocolor of M on the weave tile. pcolor is used rather than imagesc
% because the cell widths alternate between tow width and gap.
% M is (nx, ny), i.e. M(b, g) with b -> x1 and g -> x2, so it is transposed
% for pcolor. pcolor with matrix coordinates needs X, Y and C all the same
% size and draws (rows-1) by (cols-1) patches, so C is padded by one row
% and one column, which it then ignores.
[Xg, Yg] = meshgrid(xe, ye);
C = zeros(numel(ye), numel(xe));
C(1:numel(ye)-1, 1:numel(xe)-1) = M.';
pcolor(Xg, Yg, C);
shading flat;
colorbar;
axis equal tight;
xlabel('x_1 (warp)');
ylabel('x_2 (weft)');
title(ttl);
end
