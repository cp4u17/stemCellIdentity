clear all
close all

% Script to generate the model, solve the ODEs and plot the dynamics.
%  icase = 1 - homeostasis
%  icase = 2 - growing 
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

% MODEL SETTINGS
% --------------
icase = 1;
Tspan = [0 25]; X0 = zeros(12,1); X0(1:2) = 100;
switch icase
    case 1 % homeostasis
        figTg = '_hms';
    case 2 % growing
        figTg = '_grw';
end
% model definition (1 develop; 2 renewing; 3 committed)
N = 0; % initialization
ND = 2; % SCC - subcritical
lambda1 = 0.8; omega21 = 1; omega13 = 1;
TD = [omega21 2 1; omega13 1 3; omega13 1 5]; 
LD = [lambda1 1 1 2]; 
N = N + ND;
NR1 = 2; % SCC - critical
lambda1 = 1; lambda2 = 1; omega21 = 2; omega12 = 1; gamma = 1;
TR1 = [omega21 2 1; omega12 1 2]; TR1(:,2:end) = TR1(:,2:end) + N;
TR1 = [TR1; gamma 1+N 7]; 
LR1 = [lambda2 2 2 2]; LR1(:,2:end) = LR1(:,2:end) + N;
LR1 = [LR1; lambda1 1+N 1+N 8];
N = N + NR1;
NR2 = 2; % SCC - critical
lambda1 = 1; lambda2 = 1; omega21 = 2; omega12 = 1; gamma = 1;
TR2 = [omega21 2 1; omega12 1 2]; TR2(:,2:end) = TR2(:,2:end) + N;
TR2 = [TR2; gamma 1+N 8]; 
switch icase
    case 2 % changed transition (from 5-8 to 5-3)
        TR2(end,3) = 3;
end
LR2 = [lambda2 2 2 2]; LR2(:,2:end) = LR2(:,2:end) + N;
LR2 = [LR2; lambda1 1+N 1+N 10]; 
N = N + NR2;
NC1 = 1; % SCC - subcritical
lambda1 = 0.4; gamma = 1;

TC1 = [gamma 1+N 12]; 
LC1 = [lambda1 1 1 1]; LC1(:,2:end) = LC1(:,2:end) + N;
N = N + NC1;
NC2 = 2; % SCC - subcritical
lambda1 = 0.6; omega21 = 2; omega12 = 1; gamma = 1; gamma2 = 1;
TC2 = [omega21 2 1; omega12 1 2]; TC2(:,2:end) = TC2(:,2:end) + N;
TC2 = [TC2; gamma 1+N 12; gamma2 2+N 12]; 
LC2 = [lambda1 1 2 2]; LC2(:,2:end) = LC2(:,2:end) + N;
N = N + NC2;
NC3 = 2; % SCC - subcritical
lambda1 = 0.2; omega21 = 2; omega12 = 1; gamma = 1;
TC3 = [omega21 2 1; omega12 1 2]; TC3(:,2:end) = TC3(:,2:end) + N;
TC3 = [TC3; gamma 1+N 12]; 
LC3 = [lambda2 2 1 2]; LC3(:,2:end) = LC3(:,2:end) + N;
N = N + NC2;
% all SCC
N = N+1; % add death
T = [TD; TR1; TR2; TC1; TC2; TC3];
L = [LD; LR1; LR2; LC1; LC2; LC3];

% SOLVE ODEs
% ----------
[A, J] = TL2AJ(T, L, N);
odeFun = @(t, xn) J*xn;
options = odeset('Events', @(t,y) myEvents(t,y,1:1:N-1, 1e5), 'AbsTol', 1e-8, 'RelTol', 1e-8);
[time, nout] = ode45(odeFun, Tspan, X0, options);

% PLOT
% ----
indxD = 1:ND;
indxR = ND+(1:NR1+NR2);
indxC = ND+NR1+NR2+(1:NC1+NC2+NC3);
% normalization
ns = sum(X0); % sum(nout(end,[indxD indxR indxC])); % final number of cells
ts = 1/min([T(:,1); L(:,1)]);

% figure
if ~exist(fullfile(pwd,'fig'), 'dir')
    mkdir(fullfile(pwd,'fig'))
end
hfig = figure; setFigureProp(hfig)
hold on; grid on
hD = plot(time/ts, sum(nout(:,indxD), 2)/ns, '--', 'color', [0.5020    0.5020    0.5020]);
hR = plot(time/ts, sum(nout(:,indxR), 2)/ns, 'b');
hC = plot(time/ts, sum(nout(:,indxC), 2)/ns, 'k');
ht = plot(time/ts, sum(nout(:,[indxD indxR indxC])/ns, 2), 'color', [0 0.5 0], 'LineWidth', 3);
xlabel('Time [1/$\alpha_{min}$]'); 
ylabel('$n_x/n_{T,0}$');
legend('D', 'S', 'C', 'T', 'Location', 'northwest')
saveas(hfig, fullfile(pwd, 'fig', ['ExDynNcells' figTg]))