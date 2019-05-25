classdef dfplot
    % Plotting class - contains methods for plotting
    % DFPLOT.ELX = Energy level diagrams and charge carrier densities
    % DFPLOT.JT = Currents as a function of time
    % DFPLOT.JX = Total currents as a function of position
    % DFPLOT.JDDX = Drift and diffusion currents as a function of position
    % Plotting functions that are a funciton of position can accept a time
    % array as the second argument- the procedure will loop and plot the
    % solution at multiple times.
    % The third optional argument defines the x-range.

    methods (Static)

        function ELx(varargin)
            % Energy Level diagram, and charge densities plotter
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]

            % tarr is a time time array for the time you wish to plot
            if length(varargin) == 1
                sol = varargin{1};
                tarr = sol.t(end);
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)]*1e7;    % converts to nm
            elseif length(varargin) == 2
                sol = varargin{1};
                tarr = varargin{2};
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)]*1e7;    % converts to nm
            elseif length(varargin) == 3
                sol = varargin{1};
                tarr = varargin{2};
                xrange = varargin{3};
                pointtype = 't';
            end

            % Call dfana to obtain band energies and QFLs
            [u,t,x,par,dev,n,p,a,V] = dfana.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = dfana.QFLs(sol);

            xnm = x*1e7;    % x in nm for plotting

            for i = 1:length(tarr)
                % find the tarr
                p1 = find(sol.t <= tarr(i));
                p1 = p1(end);

                % Band Diagram
                FH1 = figure(1);
                PH1 = subplot(3,1,1);
                plot (xnm, Efn(p1,:), '--', xnm, Efp(p1,:), '--', xnm, Ecb(p1, :), xnm, Evb(p1 ,:));
                hold on

                % Final Charge Densities
                PH2 = subplot(3,1,2);
                semilogy(xnm, n(p1, :), xnm, p(p1, :));
                hold on

                % Ionic space charge density
                PH3 = subplot(3,1,3);
                plot(xnm, a(p1,:)-par.dev.Nion);
                hold on
            end

            figure(1)
            subplot(3,1,1);
            legend('E_{fn}', 'E_{fp}', 'CB', 'VB');
            set(legend,'FontSize',12);
            xlabel('Position [nm]');
            ylabel('Energy [eV]');
            xlim([xrange(1), xrange(2)]);
            set(legend,'FontSize',12);
            set(legend,'EdgeColor',[1 1 1]);
            hold off

            subplot(3,1,2);
            ylabel('Carrier density [cm-3]')
            legend('\itn', '\itp')
            xlabel('Position [nm]')
            xlim([xrange(1), xrange(2)]);
            ylim([1e0, 1e20]);
            set(legend,'FontSize',12);
            set(legend,'EdgeColor',[1 1 1]);
            hold off

            subplot(3,1,3);
            ylabel('Mobile ionic charge density [cm-3]');
            xlabel('Position [nm]');
            xlim([xrange(1), xrange(2)]);
            set(legend,'FontSize',12);
            set(legend,'EdgeColor',[1 1 1]);
            hold off

        end

        function Jt(sol, pos)
            % Currents as a function of time
            % POS = the readout position
            
            t = sol.t;
            [j, J] = dfana.calcJ(sol);

            figure(2);
            plot(t, J.n(:, pos),t, J.p(:, pos),t, J.a(:, pos),t, J.disp(:,pos), t, J.tot(:, pos));
            legend('Jn', 'Jp', 'Ja', 'Jdisp', 'Jtotal')
            xlabel('time [s]');
            ylabel('J [A cm^{-2}]');
            set(legend,'FontSize',16);
            set(legend,'EdgeColor',[1 1 1]);
        end

        function Jx(varargin)
            % Plots the currents
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]

            if length(varargin) == 1
                sol = varargin{1};
                tarr = sol.t(end);
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)]*1e7;    % converts to nm
            elseif length(varargin) == 2
                sol = varargin{1};
                tarr = varargin{2};
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)]*1e7;    % converts to nm
            elseif length(varargin) == 3
                sol = varargin{1};
                tarr = varargin{2};
                xrange = varargin{3};
                pointtype = 't';
            end

            xnm = sol.x*1e7;
            [j, J] = dfana.calcJ(sol);

            for i = 1:length(tarr)
                % find the time
                p1 = find(sol.t <= tarr(i));
                p1 = p1(end);

                % electron and hole currents as function of position from continuity
                figure(3)
                plot(xnm, J.n(p1, :), xnm, J.p(p1, :), xnm, J.a(p1, :), xnm, J.disp(p1, :), xnm, J.tot(p1, :))
                hold on
            end

            xlabel('Position [nm]')
            ylabel('J [A]')
            legend('Jn', 'Jp', 'Ja', 'Jdisp', 'Jtot')
            xlim([xrange(1), xrange(2)])
            hold off
        end

        function JV(JV, option)
            % JV - a solution from doJV
            % OPTION - 1 = dark only, 2 = light only, 3 = dark & light
            % JV is a structure containing dark and illuminated JVs

            if option == 1 || option == 3
                [j, J.dk.f] = dfana.calcJ(JV.dk.f);
                Vapp.dk.f = dfana.calcVapp(JV.dk.f);
                [j, J.dk.r] = dfana.calcJ(JV.dk.r);
                Vapp.dk.r = dfana.calcVapp(JV.dk.r);

                figure(4)
                plot(Vapp.dk.f, J.dk.f.tot(:,end), '--', Vapp.dk.r, J.dk.r.tot(:,end));
                hold on

            end

            if option == 2 || option == 3

                [j, J.ill.f] = dfana.calcJ(JV.ill.f);
                Vapp.ill.f = dfana.calcVapp(JV.ill.f);
                [j, J.ill.r] = dfana.calcJ(JV.ill.r);
                Vapp.ill.r = dfana.calcVapp(JV.ill.r);

                figure(4)
                plot(Vapp.ill.f, J.ill.f.tot(:,end),'--')%, 'Color', [0, 0.4470, 0.7410]);
                hold on
                plot(Vapp.ill.r, J.ill.r.tot(:,end));%,'Color', [0, 0.4470, 0.7410]);
            end

            figure(4)
            %ylim([-30e-3, 10e-3]);
            xlabel('Applied voltage [V]')
            ylabel('Current density [Acm-2]');
            hold off

        end

        function Jddx(varargin)
            % figure(5)
            % drift and diffusion currents as a function of position

            if length(varargin) == 1
                sol = varargin{1};
                tarr = sol.t(end);
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)]*1e7;    % converts to nm
            elseif length(varargin) == 2
                sol = varargin{1};
                tarr = varargin{2};
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)]*1e7;    % converts to nm
            elseif length(varargin) == 3
                sol = varargin{1};
                tarr = varargin{2};
                xrange = varargin{3};
                pointtype = 't';
            end

            % get drift and diffusion currents
            Jdd = dfana.Jddxt(sol);
            xnm = sol.x*1e7;

            figure(5)
            for i = 1:length(tarr)
                % find the time
                p1 = find(sol.t <= tarr(i));
                p1 = p1(end);

                plot(xnm, Jdd.ndiff(p1,:), xnm, Jdd.ndrift(p1,:), xnm, Jdd.pdiff(p1,:),...
                    xnm, Jdd.pdrift(p1,:), xnm, Jdd.adiff(p1,:), xnm, Jdd.adrift(p1,:));
                hold on
            end
            xlabel('Position [nm]')
            ylabel('Current density [Acm-2]')
            legend('n,diff', 'n,drift', 'p,diff', 'p,drift', 'a,diff', 'a,drift')
            hold off
        end

        function Voct(sol)
            Voc = dfana.Voct(sol);
            figure(6)
            plot(sol.t, Voc)
            xlabel('Time [s]')
            ylabel('Voc [V]')
        end

        function PLt(sol)
            PL = dfana.PLt(sol);
            figure(7)
            plot(sol.t, PL)
            xlabel('Time [s]')
            ylabel('PL [cm-5]')
        end
        
        function Vappt(sol)
            % Difference in potential between the left and right boundary
            par = sol.par
            
            if par.JV == 2
                Vapp = par.Vapp_func(par.Vapp_params, sol.t);
            else
                Vapp = -(sol.u(:,end,4)-sol.u(:,1,4)-sol.par.Vbi);
            end
            Vcell = -(sol.u(:,end,4)-sol.u(:,1,4)-sol.par.Vbi);
            
            figure(8)
            plot(sol.t, Vapp, sol.t, Vcell);
            xlabel('Time [s]')
            ylabel('V [V]')
            legend('Vapp', 'Vcell')
        end
        
        function JVapp(sol, pos)
            [j, J] = dfana.calcJ(sol);
            Vapp = -(sol.u(:,end,4)-sol.u(:,1,4)-sol.par.Vbi);
            figure(9)
            plot(Vapp, J.n(:, pos),Vapp, J.p(:, pos),Vapp, J.a(:, pos),Vapp, J.disp(:,pos), Vapp, J.tot(:, pos));
            legend('Jn', 'Jp', 'Ja', 'Jdisp', 'Jtotal')
            xlabel('Vapp [V]');
            ylabel('J [A cm^{-2}]');
            set(legend,'FontSize',16);
            set(legend,'EdgeColor',[1 1 1]);
        end
            
        function logJVapp(sol, pos)
            % plot the log of the mod J
            [j, J] = dfana.calcJ(sol);
            Vapp = dfana.calcVapp(sol);

            figure(10)
            semilogy(Vapp, abs(J.tot(:,pos)), Vapp, abs(J.n(:,pos)), Vapp, abs(J.p(:,pos)), Vapp, abs(J.a(:,pos)), Vapp, abs(J.disp(:,pos)));
            xlabel('Vapp [V]');
            ylabel('|J| [A cm^{-2}]');
            legend('Jtot', 'Jn', 'Jp', 'Ja', 'Jdisp')
            set(legend,'FontSize',16);
            set(legend,'EdgeColor',[1 1 1]);
        end
            
        function logJVapp3D(sol, pos, ylogon)
        
            t = sol.t;
            [j, J] = dfana.calcJ(sol);
            Vapp = dfana.calcVapp(sol)';
            Jtot=J.tot(:, pos);
        
            figure(11)
            surface('XData', [Vapp Vapp],             ... % N.B.  XYZC Data must have at least 2 cols
            'YData', [abs(Jtot) abs(Jtot)],             ...
            'ZData', [t' t'], ...
            'CData', [t' t'],             ...
            'FaceColor', 'none',        ...
            'EdgeColor', 'interp',      ...
            'Marker', 'none','LineWidth',1);
            s1 = gca;
            xlabel('Vapp [V]');
            ylabel('|J| [A cm^{-2}]');

            if ylogon
                set(s1,'YScale','log');
            else
                set(s1,'YScale','linear');
            end
            hold off
            
        end
        
        function xmesh(sol)
           figure(11)
           plot(sol.x)
           xlabel('Position [cm]')
           ylabel('Point')
        end
            
        
        % multiplot 1
        function mp1(varargin)

            % tarr is a time time array for the time you wish to plot
            if length(varargin) == 1
                solstruct = varargin{1};
                tarr = solstruct.t(end);
                pointtype = 't';
            elseif length(varargin) == 2
                solstruct = varargin{1};
                tarr = varargin{2};
                pointtype = 't';
            elseif length(varargin) == 3
                solstruct = varargin{1};
                pointtype = varargin{2};
                tarr = varargin{3};
            end

        end

    end

end
