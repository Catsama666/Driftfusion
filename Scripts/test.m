% A test script for comparing the voltage step measurement of a single schottky junction  
% in perovskite solar cell
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
%% Start code
tic
initialise_df

%% Input parameters
params_filepath = './PEM_workshop_Input_files/Schottky.csv';     % Filepath to the parameters file
output_filename = 'schottky';                                    % Filename for output file

tdwell = 1e-8; % TDWELL - the relaxation time after jumping
mobseti = 0; % MOBSETI - determines whether the ion mobility is switched on after the jump
Int = 0; % INT - Intensity during stabilisation phase
stabilise = 0; % STABILISE - Determines whether the time step is increased such that the
% device is at steady-state at the end of the solution
accelerate = 0; % ACCELERATE - Accelerate the ions to be at the same mobility as the
% electrons (for easy access to steady state solutions)
Vjump = 0; % VJUMP - the voltage to be jumped to
Vstep = 0.5; % Vstep - the step voltage add each time
% Vinitial = Vjump+Vstep; % Vinitial - the start voltage of the loop;
Vmax = 1; % Vmax - the max Vstep voltage applied
n = 1; 
m = 1;
%% In a single loop, Load in parameters

par = pc(params_filepath);

%% Obtain equilibrium solution
sol_eq = equilibrate(par, 1);
sol_relax(n) = sol_eq.el;

while(n<=Vmax/Vstep)
  n=n+1;
  Vjump=Vjump+Vstep;
  sol_relax(n) = jumptoV(sol_relax(n-1), Vjump, tdwell, mobseti, Int, stabilise, accelerate);
  


end

while(m<n)
  m=m+1;
  t = sol_relax(m).t;
  [J, j, xmesh] = dfana.calcJ(sol_relax(m));
  ppos = getpointpos(sol_relax(m).x(end), xmesh);
  deltaQ = cumtrapz(t, J.tot(:, ppos));
  figure(m+100)
  plot(t, J.n(:, ppos),t, J.p(:, ppos),t, J.a(:, ppos),t, J.c(:, ppos), t, J.disp(:,ppos), t, J.tot(:, ppos));
  legend('Jn', 'Jp', 'Ja', 'Jc', 'Jdisp', 'Jtotal')
  xlabel('time [s]');
  ylabel('J [A cm^{-2}]');
  set(legend,'FontSize',16);
  set(legend,'EdgeColor',[1 1 1]);
  
  dfplot.ELx(sol_eq.el,sol_eq.el.t(end));
  hold on;
  dfplot.ELx(sol_relax(m),sol_relax(m).t(end));
  hold off;

end





% %% Load in parameters
% par = pc(params_filepath);
% 
% %% Obtain equilibrium solution
% sol_eq = equilibrate(par, 1);
% 
% %% Call function to obtain equilibrium and perform voltage step 500mV for 1e-8s
% sol_relax = jumptoV(sol_eq.el, Vjump, tdwell, mobseti, Int, stabilise, accelerate);
% t = sol_relax(n).t;
% [J, j, xmesh] = dfana.calcJ(sol_relax(n));
% ppos = getpointpos(sol_relax(n).x(end), xmesh);
% deltaQ = cumtrapz(t, J.tot(:, ppos));
% 
% %% Plot the current density diagram, the energy diagram after voltage step
% 
% while(n<=(Vmax-Vinitial)/Vstep+1)
%   n=n+1;
%   Vjump=Vjump+Vstep;
%   sol_relax(n) = jumptoV(sol_relax(n-1), Vjump, tdwell, mobseti, Int, stabilise, accelerate);
%   t = sol_relax(n).t;
%   [J, j, xmesh] = dfana.calcJ(sol_relax(n));
%   ppos = getpointpos(sol_relax(n).x(end), xmesh);
%   deltaQ = cumtrapz(t, J.tot(:, ppos));
%   figure(n+100)
%   plot(t, J.n(:, ppos),t, J.p(:, ppos),t, J.a(:, ppos),t, J.c(:, ppos), t, J.disp(:,ppos), t, J.tot(:, ppos));
%   legend('Jn', 'Jp', 'Ja', 'Jc', 'Jdisp', 'Jtotal')
%   xlabel('time [s]');
%   ylabel('J [A cm^{-2}]');
%   set(legend,'FontSize',16);
%   set(legend,'EdgeColor',[1 1 1]);
% 
%   
%   
% %   dfplot.ELx(sol_eq.el,sol_eq.el.t(end));
% %   hold on;
% %   dfplot.ELx(sol_relax(n),sol_relax(n).t(end));
% %   hold off;
% 
% end














