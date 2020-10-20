clear all
close all

% Script to run and plot the feedback model dynamics and eigenvalues
% 
% Copyright 2020 Cristina Parigini
% 
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
% 
%      http://www.apache.org/licenses/LICENSE-2.0
% 
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

% parameters
fileName = 'feedbackModel';
load(fullfile(pwd, fileName))
Tspan0 = [0 100]; Tspan1 = Tspan0(end)*[1/2 1];
ns0 = 1e3; % ref state

% run simulations
xn0 = ns0/3*[1; 1; 1];
[t_hms1, n_hms1, nt_hms1, xeig_hms1, minRate_hms1] = fbDisruption(L, T, Lfb, Tfb, Tspan0, xn0);
nf_hms1 = n_hms1(end,:)';  % [79.5113, 265.4589, 655.0298]
% B
xn0 = [200; 500; 520];
[t_hms2, n_hms2, nt_hms2, xeig_hms2] = fbDisruption(L, T, Lfb, Tfb, Tspan1, xn0);
nf_hms2 = n_hms2(end,:)';
% C
xn0 = [10; 150; 700];
[t_hms3, n_hms3, nt_hms3, xeig_hms3] = fbDisruption(L, T, Lfb, Tfb, Tspan1, xn0);
nf_hms3 = n_hms3(end,:)';

% output folder
if ~exist(fullfile(pwd, 'fig'),'dir')
   mkdir(fullfile(pwd, 'fig'))
end

% set axis properties and legend for each curve
tf = 1/minRate_hms1; xax = 'Time [1/$\alpha_{min}$]'; XLim = [0 10]; % xlabel
ns = nt_hms1(end); yax_n = '$N^x_a/N^H_{a,0}$'; YLim_n = [0.8 1.4]; % ylabel (total numbe of cells)
yax_eig = '$\mu_a$'; YLim_eig = [-0.18 0.18]; % ylabel (max eigenvalue)
xleg_hms = 'H';
xleg_1 = 'P_1'; 
xleg_2 = 'P_2'; 

% TOTAL NUMBER OF CELLS
hfig1 = figure; setFigureProp(hfig1)
hold on; grid on
hh = plot((t_hms1-Tspan1(1))/tf, nt_hms1/ns);
h1 = plot((t_hms2-Tspan1(1))/tf, nt_hms2/ns);
h2 = plot((t_hms3-Tspan1(1))/tf, nt_hms3/ns);
xlabel(xax); ylabel(yax_n)
legend([hh, h1, h2], ...
    {xleg_hms, xleg_1, xleg_2}, 'Location', 'northeast')
set(gca, 'XLim', XLim, 'YLim', YLim_n)

% EIGENVALUE
hfig2 = figure; setFigureProp(hfig2)
hold on; grid on
hh = plot((t_hms1-Tspan1(1))/tf, max(xeig_hms1,[],2));
h1 = plot((t_hms2-Tspan1(1))/tf, max(xeig_hms2,[],2));
h2 = plot((t_hms3-Tspan1(1))/tf, max(xeig_hms3,[],2));
xlabel(xax); ylabel(yax_eig)
set(gca, 'XLim', XLim, 'YLim', YLim_eig)
legend([hh, h1, h2], ...
    {xleg_hms, xleg_1, xleg_2}, 'Location', 'northeast')

% save figures
saveas(hfig1, fullfile(pwd, 'fig', 'fbExample_steadyState'))
saveas(hfig2, fullfile(pwd, 'fig', 'fbExample_eigenvalue'))