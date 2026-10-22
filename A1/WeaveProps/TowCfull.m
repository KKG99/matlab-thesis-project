function [Cw, aw, Teps] = TowCfull(Cl, al, ThetaZ, PhiY)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Purpose: Rotate a tow stiffness and CTE from tow axes into weave axes.
%          The tow is transversely isotropic about its own axis 1.
% Input:
% - Cl: 6x6 tow stiffness in tow axes
%       (Voigt order [11 22 33 23 13 12], engineering shear strains)
% - al: 6x1 tow CTE in tow axes
% - ThetaZ: rotation about weave x3 [deg]  (0 = warp, 90 = weft)
% - PhiY:   crimp angle [deg]             (+ up, - down, 0 flat)
% Output:
% - Cw: 6x6 stiffness in weave axes
% - aw: 6x1 CTE in weave axes
% - Teps: 6x6 map taking a strain in weave axes to the same strain in tow
%         axes (eps_tow = Teps*eps_weave). Needed to check the tow's own
%         allowables against the strain it actually sees.
%
% Tow axes expressed in weave coordinates:
%   e1 = [cos(phi)cos(theta), cos(phi)sin(theta), sin(phi)]
%   e2 = [-sin(theta), cos(theta), 0]
%   e3 = e1 x e2
%
% The transformation matrices are built numerically by rotating the six
% Voigt basis states, so no transformation-matrix convention is assumed.
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

th = ThetaZ*pi/180;
ph = PhiY  *pi/180;

e1 = [cos(ph)*cos(th), cos(ph)*sin(th), sin(ph)];
e2 = [-sin(th), cos(th), 0];
e3 = cross(e1, e2);
a  = [e1; e2; e3];          % a(i,j) = (tow axis i) . (weave axis j)

% -- Voigt index -> tensor index pair
P = [1 1; 2 2; 3 3; 2 3; 1 3; 1 2];

% -- Build Tsig (stress) and Teps (engineering strain), 6x6
Tsig = zeros(6);
Teps = zeros(6);
for k = 1:6
    % stress basis: sigma_pq = 1
    S = zeros(3); S(P(k,1),P(k,2)) = 1; S(P(k,2),P(k,1)) = 1;
    Sr = a*S*a.';
    % strain basis: eps_voigt(k) = 1  ->  eps_pq = 1 (normal) or 0.5 (shear)
    E = zeros(3);
    if P(k,1) == P(k,2)
        E(P(k,1),P(k,2)) = 1;
    else
        E(P(k,1),P(k,2)) = 0.5;  E(P(k,2),P(k,1)) = 0.5;
    end
    Er = a*E*a.';
    for m = 1:6
        Tsig(m,k) = Sr(P(m,1),P(m,2));
        if P(m,1) == P(m,2)
            Teps(m,k) = Er(P(m,1),P(m,2));
        else
            Teps(m,k) = 2*Er(P(m,1),P(m,2));
        end
    end
end

% -- sigma_tow = C_tow * eps_tow,  with  sigma_tow = Tsig*sigma_weave and
%    eps_tow = Teps*eps_weave   =>   C_weave = Tsig \ C_tow * Teps
Cw = Tsig \ (Cl * Teps);

% -- CTE transforms like a (symmetric) second-rank tensor, i.e. like strain
aw = Teps \ al;

end
