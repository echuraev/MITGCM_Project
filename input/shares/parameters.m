function h = parameters()
%parameters Calculats area parameters
% This function returns a structure hander which has the whole set of
% parameters using for calculations and results depiction.

    h.deb = 1;

    %** Domain dimmension parameters
    %*  begin
    % Dimensions of the grid
    h.nx=1000; %2560; % - 0.0025;
    h.ny=1;
    h.nz=212; %  - 0.00125
              % 106 - 0.0025
    %* end

    %** Experiment description
    %* begin
    % |                               |  |--------|______________________|
    % |-------------------------------|  |  V+h2  |   r2       h2        |
    % |     r2      h2                |  |        |----------------------|
    % |-------------------------------|  | r2 h2l |                      |
    % |                               |->|        |dH     r1 h1r         |
    % |          V1                   |  |--------|                      |
    % |     r1      h1                |  | r1 h1l |         V1           |
    % |_______________________________|  |_______________________________|
    % where h, r - level and density of corresponding part of liquid, dH -
    % difference of r1 water levels.
    %
    % h2l can be calculated from volume V, which was added for level changing.
    % h2l = V/(r2*Lx*Lw), where V - volume in liters
    %
    % The hydrostatic balance after water adding:
    % r2*h2l = r2*h2 + r1*dH => dH = r2*(h2l - h2)/r1
    %
    % Just calculate a volume before water adding and after that.
    % V1 - volume of water with density r1. It is constant.
    %  /
    %  | V1 = h1*Lx*Lw
    % < 
    %  | V1 = h1l*Lg*Lw + (dH + h1l)*(Lx - Lg)*Lw
    %  \
    %
    % Express h1l from the system above.
    % h1*Lx = h1l*Lg + (dH + h1l)*(Lx - Lg)
    % h1*Lx = h1l*Lg + h1l*(Lx - Lg) + dH*(Lx - Lg)
    % h1l = (h1*Lx - dH*(Lx - Lg))/Lx
    %
    % h1l = h1 - dH*(1 - Lg/Lx)
    % h1r = h1l + dH
    %
    % Water depth to the left and to the right of the gate:
    % Hl = h1l + h2l = 0.2679
    % Hr = h1l + dH + h2 = 0.2642
    %
    % We'll decrease h2l by (Hr - Hl) and then add it to the left surface
    % elevation hight.
    %
    %* end

    %** Area parameters calculating using algorithm in section above.
    %* begin
    % Size of X domain
    h.Lx=6.4;
    % Size of Y domain
    h.Ly=0.4;
    %* Length of the domain to the left of the gate
    h.Lg = 0.6; % meters
    %* Width of the domain
    h.Lw = 0.4; % meters
    % Density of the layers
    h.r1 = 1047.0;
    h.r2 = 1022.0;
    % Reference dencity
    h.r0 = 1024.8;
    % Volume of added water
    h.V=38;
    % Initial water levels
    h.h1 = 0.2;
    h.h2 = 0.05;
    % Depth of the added liquid.
    h.h2add = (h.V/h.r2)/(h.Lg*h.Lw);
    % h2add  = 0.158; - In article
    % h2add  = 0.1549; - We have
    % Target h2 levels calculating.
    h.h2r = h.h2;
    h.h2l = h.h2 + h.h2add;
    %h2r = 0.05; % meters - In article
    %h2r = 0.05; % meters - We have
    %h2l = 0.21; % meters - In article
    %h2l = 0.2049; % meters - We have (It gives less amplitude of a soliton)
    % Delta of r1 levels calculating.
    h.dH = h.r2*(h.h2l - h.h2)/h.r1;
    % r1 level to the left and to the right of gate calculating.
    h.h1l = h.h1 - h.dH*(1 - h.Lg/h.Lx);
    h.h1r = h.h1l + h.dH;
    % The target depth of water to the laft and to the rigtht of gate
    % calculating.
    h.Hl = h.h2l + h.h1l;
    h.Hr = h.h2r + h.h1r; % or Hr = h2 + h1r;
    h.H = h.Hl;
    % Common level is H, but surface to the left of the gate is 
    % increased on (Hl - Hr). And depth of the left pycnocline depth will be
    % decreased. It's very debatable! Lets keep h2l as calculated for the time being!
    h.dSurf = h.Hl - h.Hr;
    % h2l = h2l - dSurf;
    % Rest of parameters setting
    % Pycnocline thickness
    %dh=0.015; % meters - In articlce
    %dh=0.01; % meters - We have
    % NOTE! dh is just half of thickness, that means the real thickness is 2*dh.
    % in the experiment 0.02 - 0.35.
    h.dh = 0.01;
    % Gate thickness
    % But the real thickness is 0.005.
    % dg=0.005; % - In article
    % dg=0.005; % - We have
    % NOTE! The real thickness of the interafce is two times thicker! The real
    % gates are 0.005 meters + mixing give us approximately 0.01 m.
    h.dg = 0.005; % meters
    % Pycnoclines depth
    h.zpl = -h.h2l;
    h.zpr = -h.h2r;
    % Middle and delta of density
    h.r_m = (h.r1 + h.r2)/2;
    h.dr = (h.r1 - h.r2)/2;
    % Acceleration
    h.g = 9.81;
    h.nu = 1.0e-5;
    %* end

    %** Additional constants calculations.
    %* begin
    % internal long wave velocity
    h.c0 = sqrt(h.g*((h.r1 - h.r2)/(h.r1))*((h.h1r*h.h2r)/(h.h1r + h.h2r)));
    h.c_max = sqrt(h.g*(h.h1r + h.h2r)*(h.r1 - h.r2)/((sqrt(h.r1) + sqrt(h.r2))^2));
    % Internal soliton amplitude
    h.a = (h.h1r*sqrt(h.r2) - h.h2r*sqrt(h.r1))/(sqrt(h.r1) + sqrt(h.r2));
    % Wave Reynolds number
    h.ReW = h.c0*h.H/h.nu; % ~2.6*10e4
    %* end
end

