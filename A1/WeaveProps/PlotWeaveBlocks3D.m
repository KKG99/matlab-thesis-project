function hfig = PlotWeaveBlocks3D(state, TitleStr, hfig)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Purpose: Draw the weave RUC as its real 3-D block grid (e.g. 4x4x4 = 64
%          blocks for RUC 401), coloured by material and with FAILED blocks
%          in red - the 3-D analogue of NASA's "Subcell Failure" panel.
%          warp '5' = blue, weft '6' = green, matrix 'M' = light grey,
%          failed   = red (drawn opaque so it stands out).
% Input:
% - state: damage state from WeaveDamageInit / WeavePlyDamage / SimWeaveDamage
%          (needs wruc, nx, ny, nz, mSt, mDam)
% - TitleStr (optional), hfig (optional existing figure handle)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if nargin < 2 || isempty(TitleStr), TitleStr = 'Weave RUC - failed blocks'; end
wruc = state.wruc;
nx = state.nx; ny = state.ny; nz = state.nz;

% -- Cumulative edge coordinates
xe = [0, cumsum(wruc.d)];
ye = [0, cumsum(wruc.h)];
ze = [0, cumsum(wruc.l)];

if nargin < 3 || isempty(hfig), hfig = figure; else, figure(hfig); end
clf(hfig); hold on;

for b = 1:nx
    for g = 1:ny
        for k = 1:nz
            c = wruc.matsCh(b,g,k);
            if c == 'M'
                failed = state.mDam(b,g,k);
            else
                failed = any(state.mSt{b,g,k}(:) == '123');
            end

            if failed
                col = [0.85 0 0];  a = 1.0;      % red, opaque
            elseif c == '5'
                col = [0.10 0.35 0.90]; a = 0.55; % warp blue
            elseif c == '6'
                col = [0.10 0.75 0.40]; a = 0.55; % weft green
            else
                col = [0.85 0.85 0.85]; a = 0.25; % matrix, faint
            end

            drawbox(xe(b), xe(b+1), ye(g), ye(g+1), ze(k), ze(k+1), col, a);
        end
    end
end

axis equal; axis tight; grid on; view(35, 25);
xlabel('x_1 (warp)'); ylabel('x_2 (weft)'); zlabel('x_3 (layers)');
title(TitleStr);
hold off;
drawnow;
end

%=========================================================================%
function drawbox(x0,x1,y0,y1,z0,z1,col,alpha)
X = [x0 x1 x1 x0  x0 x1 x1 x0];
Y = [y0 y0 y1 y1  y0 y0 y1 y1];
Z = [z0 z0 z0 z0  z1 z1 z1 z1];
F = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
patch('Faces',F,'Vertices',[X(:) Y(:) Z(:)], ...
    'FaceColor',col,'FaceAlpha',alpha,'EdgeColor','k','LineWidth',0.5);
end
