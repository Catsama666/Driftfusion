% An example script for comparing the current-voltage curves for two different
% architects of perovskite solar cell
clear all
close all

initialise_df

% Create parameters objects for Spiro/MAPI/TiO2 and PEDOT/MAPI/PCBM devices
par.tio2 = pc('Input_files/spiro_mapi_tio2.csv');
par.pcbm = pc('Input_files/ptpd_mapi_pcbm.csv');

% Find equilibrium solutions
soleq.tio2 = equilibrate(par.tio2);
soleq.pcbm = equilibrate(par.pcbm);

% Perform current-voltage scans to 1.2 V
JV.tio2.kscan50em3 = doJV(soleq.tio2.ion, 50e-3, 100, 1, 1, 0, 1.2, 3);
JV.pcbm.kscan50em3 = doJV(soleq.pcbm.ion, 50e-3, 100, 1, 1, 0, 1.2, 3);

% plot the current voltage curve
dfplot.JV(JV.tio2.kscan50em3,3)

figure(4)
hold on

dfplot.JV(JV.pcbm.kscan50em3,3)
legend('tio2-dk-f', 'tio2-dk-r', 'tio2-1sun-f', 'tio2-1sun-r',...
        'pcbm-dk-f', 'pcbm-dk-r', 'pcbm-1sun-f', 'pcbm-1sun-r')
xlim([0, 1.2])
ylim([-30e-3, 10e-3])
hold off

% plot the energy level diagram and carrier densities for the tio2 device at
% 1 V (t= 20s) during the illuminated forward scan
dfplot.ELx(JV.tio2.kscan50em3.ill.f, 20)

% plot the currents as a function of position in the PCBM device at 0.5 V
% (t = 10 s) during the illuminated forward scan
dfplot.Jx(JV.pcbm.kscan50em3.ill.f, 10)

% Save the workspace- this is commented out as the filepath should lead to
% a folder on your computer. It is not recommended to store large files in
% the Github repository as they will be cached even when deleted leading to
% a large folder sizes
% save('/Users/Username/Data/temp.mat')
