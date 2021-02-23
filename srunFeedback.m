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
caseID = 'unstable';

switch caseID
    case 'stable'
        D = 3; ns = 1e3; K = 1; N = 2;
        r11 = 1; r12 = 1; r13 = 0;
        r31 = 0; r32 = 2; r33 = 0;
        lambda1 = @(n) 0.005 + 1.02./(K^N + (n/ns).^N);
        lambda3 = @(n) 0.005 + 1.01./(K^N + (n/ns).^N);
        omega13 = @(n) 0.005 + 1.12./(K^N + (n/ns).^N);
        omega23 = @(n) 0.005 + 1.85./(K^N + (n/ns).^N);
        omega21 = @(n) 0.005 + 1.52*(n/ns).^N./(K^N + (n/ns).^N);
        delta1 = @(n) 0.005 + 2.68*(n/ns).^N./(K^N + (n/ns).^N);
        delta2 = @(n) 0.005 + 1.64*(n/ns).^N./(K^N + (n/ns).^N);
        % matrix A elements
        a11 = @(n) (r11-1)*lambda1(n) - omega13(n) - delta1(n);
        a22 = @(n) -omega21(n) -omega23(n) - delta2(n);
        a33 = @(n) (r33-1)*lambda3(n) - omega21(n);
        a21 = @(n) r12*lambda1(n);
        a31 = @(n) r13*lambda1(n) + omega13(n);
        a12 = @(n) omega21(n);
        a32 = @(n) omega23(n);
        a13 = @(n) r31*lambda3(n);
        a23 = @(n) r32*lambda3(n);
        AFnc = @(n) [a11(n) a12(n) a13(n); a21(n) a22(n) a23(n); a31(n) a32(n) a33(n)];
        % reference timescale
        alpha_ref = @(ns) min([lambda1(ns) lambda3(ns) omega13(ns) omega21(ns) omega23(ns) delta1(ns) delta2(ns)]);
    case 'unstable'
        D = 2; ns = 1e3; K = 1; 
        r11 = 0.09; r12 = 0.11; r21 = 0.47; r22 = 0.53;
        lambda1 = @(n) 0.1 + 0.68./(K^50+ (n/ns).^50);
        lambda2 = @(n) 0.005 + 13.54./(K^2 + (n/ns).^2);
        omega12 = @(n) 0.001 + 0.18./(K^5 + (n/ns).^5);
        omega21 = @(n) 1 + 0.87*(n/ns)./(K + (n/ns));
        delta1 = @(n) 0.92*ones(size(n));
        % matrix A elements
        a11 = @(n) (2*r11-1)*lambda1(n) - omega12(n) - delta1(n);
        a22 = @(n) (2*r22-1)*lambda2(n) - omega21(n);
        a21 = @(n) 2*r12*lambda1(n) + omega12(n);
        a12 = @(n) 2*r21*lambda2(n) + omega21(n);
        AFnc = @(n) [a11(n) a12(n); a21(n) a22(n)];
        % reference timescale
        alpha_ref = @(ns) min([lambda1(ns) lambda2(ns) omega12(ns) omega21(ns) delta1(ns)]);
end

% dynamics and eigenvalues
dynFnc = @(t, n) AFnc(sum(n))*n;
eigFnc = @(n) max(real(eig(AFnc(n))));
xs2 = fsolve(@(n) dynFnc([], n), ns*ones(D,1));
xns = [zeros(D,1) xs2]; ns = sum(xns);
tscale = 1/alpha_ref(ns(2));
nscale = ns(2);
muscale = alpha_ref(ns(2));

% Jacobian (linearized dynamics)
Tol = 1e-5;
J = zeros(D, D, length(ns));
xeig = zeros(D, length(ns));
for ii = 1:length(ns)
    a11p = numjac(@(t, x) a11(x), 0, ns(ii), a11(ns(ii)), Tol,[]); a11p(isnan(a11p)) = 0;
    a12p = numjac(@(t, x) a12(x), 0, ns(ii), a12(ns(ii)), Tol,[]); a12p(isnan(a12p)) = 0;
    a21p = numjac(@(t, x) a21(x), 0, ns(ii), a21(ns(ii)), Tol,[]); a21p(isnan(a21p)) = 0;
    a22p = numjac(@(t, x) a22(x), 0, ns(ii), a22(ns(ii)), Tol,[]); a22p(isnan(a22p)) = 0;
    if D == 3
        a13p = numjac(@(t, x) a13(x), 0, ns(ii), a13(ns(ii)), Tol,[]); a13p(isnan(a13p)) = 0;
        a23p = numjac(@(t, x) a23(x), 0, ns(ii), a23(ns(ii)), Tol,[]); a23p(isnan(a23p)) = 0;
        a31p = numjac(@(t, x) a31(x), 0, ns(ii), a31(ns(ii)), Tol,[]); a31p(isnan(a11p)) = 0;
        a32p = numjac(@(t, x) a32(x), 0, ns(ii), a32(ns(ii)), Tol,[]); a32p(isnan(a12p)) = 0;
        a33p = numjac(@(t, x) a33(x), 0, ns(ii), a33(ns(ii)), Tol,[]); a33p(isnan(a13p)) = 0;
        J(:,:,ii) = AFnc(ns(ii)) + [a11p a12p a13p; a21p a22p a23p; a31p a32p a33p]*xns(:,ii)*ones(1,D);
    else
        J(:,:,ii) = AFnc(ns(ii)) + [a11p a12p; a21p a22p]*xns(:,ii)*ones(1,D);
    end
    xeig(:,ii) = eig(J(:,:,ii));
end

% Integration and plot of the dynamics
Tspan = [0 4*tscale];
% x0 = [xns(:,1) xns(:,2)*[0.8 1 1.2]];
% xleg = {'N_a(0) = 0', 'N_a(0) = 0.8N_a^*', 'N_a(0) = N_a^*', 'N_a(0) = 1.2N_a^*'};
x0 = [xns(:,2)*[0.8 1.2 1]];
xleg = {'{\bf n}_a(t=0) = 0.8 {\bf n}_a^*', '{\bf n}_a(t=0) = 1.2 {\bf n}_a^*', '{\bf n}_a(t=0) = {\bf n}_a^*'};
xcol = lines(length(xleg));

% set figures settings
hdyn = figure; setFigureProp(hdyn); hold on; grid on
xlabel('Time $[1/\alpha_{min}]$'); ylabel('$N_a/N_a^*$')
% plot(Tspan/tscale, ns(2)/nscale*[1 1], 'k--', 'linewidth', 3); % plot N*
heig = figure; setFigureProp(heig); hold on; grid on
xlabel('Time $[1/\alpha_{min}]$'); ylabel('$\mu_a [\alpha_{min}]$')
% phase plot
if D == 2
    hdyn1 = figure; setFigureProp(hdyn1); hold on; grid on
    xlabel('$n_1/N_a^*$'); ylabel('$n_2/N_a^*$')
end
h1 = zeros(1, size(x0,2)); h2 = zeros(1, size(x0,2));

% loop on initial conditions
eventFnc = @(t, y) myEvents(t, y, 1:D, 1e5);
options = odeset('AbsTol', 1e-10, 'RelTol', 1e-10, 'Events', eventFnc);
for ii = 1:size(x0,2)
    [t,y,te,ye,ie] = ode45(dynFnc, Tspan, x0(:,ii)+1e-4, options);
    % N VS time
    figure(hdyn)
    h1(ii) = plot(t/tscale, sum(y, 2)/nscale, 'linewidth',2, 'Color', xcol(ii,:));
    % eigenvalue
    figure(heig)
    muA = zeros(1,length(t));
    for jj = 1:length(t)
        muA(jj) = eigFnc(sum(y(jj,:)));
    end
    h2(ii) = plot(t/tscale, muA/muscale, 'linewidth',2, 'Color', xcol(ii,:));
    if D == 2 % n1 VS n2
        figure(hdyn1)
        plot(y(:,1)/nscale, y(:,2)/nscale, 'linewidth',2, 'Color', xcol(ii,:));
    end
end
% set leged
figure(hdyn); set(gca, 'FontSize', 20)
legend(h1([1 3 2]), xleg([1 3 2]), 'FontSize', 12)
figure(heig); set(gca, 'FontSize', 20)
legend(h2([1 3 2]), xleg([1 3 2]), 'FontSize', 12)
% axis([0 4 -0.8 0.8])
if D == 2 
    figure(hdyn1); set(gca, 'FontSize', 10)
end

% save figures
saveas(hdyn, fullfile(pwd, 'fig', ['fbExample_steadyState_' caseID]))
saveas(heig, fullfile(pwd, 'fig', ['fbExample_eigenvalue_' caseID]))
if D == 2 
    saveas(hdyn1, fullfile(pwd, 'fig', ['fbExample_phasePlot_' caseID]))
end