function hfig = PlotLamVonMises3D(states, Sbply, Orient, TitleStr, hfig)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Purpose: Draw the WHOLE laminate as one 3-D block grid and colour every
%          block by its von Mises stress, so you can see exactly where the
%          highest stresses concentrate through all plies and all layers.
%          Plies are stacked in z (ply 1 at the bottom). One colourbar for
%          the whole laminate; the black arrow marks the load direction (+x1).
% Input:
% - states: cell of per-ply damage states (wruc, nx, ny, nz)
% - Sbply : cell of per-ply block stresses Sb{b,g,k}
% - Orient: ply angles (for the title)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if nargin < 4 || isempty(TitleStr), TitleStr = 'Laminate von Mises (whole stack)'; end
nx = states{1}.nx; ny = states{1}.ny; nzl = states{1}.nz;
Nplies = numel(states);
wruc = states{1}.wruc;

xe = [0, cumsum(wruc.d)];
ye = [0, cumsum(wruc.h)];
le = [0, cumsum(wruc.l)];
H  = le(end);                       % one ply thickness

if nargin < 5 || isempty(hfig), hfig = figure; else, figure(hfig); end
clf(hfig); hold on;

vmax = 0;
% -- first pass: global max for the colour scale
for i = 1:Nplies
    for b = 1:nx
        for g = 1:ny
            for k = 1:nzl
                s = Sbply{i}{b,g,k};
                vmax = max(vmax, vm(s));
            end
        end
    end
end
if vmax == 0, vmax = 1; end

% -- second pass: draw every block, coloured by vm
for i = 1:Nplies
    for b = 1:nx
        for g = 1:ny
            for k = 1:nzl
                s = Sbply{i}{b,g,k};
                v = vm(s);
                z0 = (i-1)*H + le(k);
                z1 = (i-1)*H + le(k+1);
                vbox(xe(b), xe(b+1), ye(g), ye(g+1), z0, z1, v);
            end
        end
    end
end

colormap(jet); caxis([0 vmax]); colorbar;
axis equal; axis tight; grid on; view(35, 25);
xlabel('x_1 (load)'); ylabel('x_2'); zlabel('through-thickness (plies)');
% -- load direction arrow along +x1
quiver3(0, 0, Nplies*H*1.05, 0.6*xe(end), 0, 0, 0, 'k', 'LineWidth', 2);
title(TitleStr);
hold off; drawnow;
end

%=========================================================================%
function [v] = vm(s)
v = sqrt((s(1)-s(2))^2 + (s(2)-s(3))^2 + (s(1)-s(3))^2 + ...
    6*s(4)^2 + 6*s(5)^2 + 6*s(6)^2)/sqrt(2);
end

%=========================================================================%
function vbox(x0,x1,y0,y1,z0,z1,val)
X = [x0 x1 x1 x0  x0 x1 x1 x0];
Y = [y0 y0 y1 y1  y0 y0 y1 y1];
Z = [z0 z0 z0 z0  z1 z1 z1 z1];
F = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
patch('Faces',F,'Vertices',[X(:) Y(:) Z(:)], ...
    'FaceColor','flat','CData',val*ones(8,1), ...
    'EdgeColor','k','LineWidth',0.3,'FaceAlpha',0.95);
end
