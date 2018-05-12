function [x] = meshgen_x(p)

switch p.xmesh_type
    % Linearly spaced
    case 1
        x = linspace(0, p.xmax, p.pp+p.pii+p.pn);
    
    % Linearly spaced, more points at interfaces
    case 2
        x = [linspace(0, p.tp-p.tint, p.pp),...
            linspace(p.tp-p.tint+p.deltax, p.tp, p.pint),...
            linspace(p.tp+p.deltax, p.tp+p.tint, p.pint),...
            linspace(p.tp+p.tint+p.deltax, p.tp+p.ti-p.tint-p.deltax, p.pii),...
            linspace(p.tp+p.ti-p.tint, p.tp+p.ti-p.deltax, p.pint),...
            linspace(p.tp+p.ti, p.tp+p.ti+p.tint-p.deltax, p.pint),...
            linspace(p.tp+p.ti+p.tint, p.xmax, p.pn)];

    % Linearly spaced, more points at interfaces and electrodes
    case 3
        x = [linspace(0, p.te, p.pepe),...
            linspace(p.te+p.deltax, p.tp-p.tint, p.pp),...
            linspace(p.tp-p.tint+p.deltax, p.tp, p.pint),...
            linspace(p.tp+p.deltax, p.tp+p.tint, p.pint),...
            linspace(p.tp+p.tint+p.deltax, p.tp+p.ti-p.tint-p.deltax, p.pii),...
            linspace(p.tp+p.ti-p.tint, p.tp+p.ti-p.deltax, p.pint),...
            linspace(p.tp+p.ti, p.tp+p.ti+p.tint-p.deltax, p.pint),...
            linspace(p.tp+p.ti+p.tint, p.xmax-p.te-p.deltax, p.pn),...
            linspace(p.xmax-p.te, p.xmax, p.pepe);];
    
    % More points in SCR and interfaces
    case 4
        x = [linspace(0, p.tp-p.tscr, p.pp),...
            linspace(p.tp-p.tscr+(p.tscr/p.pscr), p.tp-p.tint, p.pscr),...
            linspace(p.tp-p.tint+(p.tint/p.pint), p.tp, p.pint),...
            linspace(p.tp+(p.tint/p.pint), p.tp+p.tint, p.pint),...
            linspace(p.tp+p.tint+(p.tint/p.pint), p.tp+p.tscr, p.pscr),...
            linspace(p.tp+p.tscr+(p.tscr/p.pscr), p.tp+p.ti-p.tscr-(p.tscr/p.pscr), p.pii),...
            linspace(p.tp+p.ti-p.tscr, p.tp+p.ti-p.tint-(p.tint/p.pint), p.pscr),...
            linspace(p.tp+p.ti-p.tint, p.tp+p.ti-(p.tint/p.pint), p.pint),...
            linspace(p.tp+p.ti, p.tp+p.ti+p.tint-(p.tint/p.pint), p.pint),...
            linspace(p.tp+p.ti+p.tint, p.tp+p.ti+p.tscr-(p.tscr/p.pscr), p.pscr),...
            linspace(p.tp+p.ti+p.tscr, p.xmax, p.pn);];
    otherwise
       error('DrIFtFUSION:xmesh_type', [mfilename ' - xmesh_type not recognized'])
end

end