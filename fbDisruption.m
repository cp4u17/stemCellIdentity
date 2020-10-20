function [time, nout, ntout, xeig, minRate] = fbDisruption(Lin, Tin, Lfb, Tfb, Tspan, xn0)

% This function build the model and integrate the ODEs for the Feedback
% Disruption model. Feedback disruption is applied to all cells.
%
% See notes_feedbackTestCase for more details.

% parameters
indxAS = [1 2 3]; N = 4;
[T, L] = fbDisruption_model(Lin, Tin, Lfb, Tfb);
T = symb2input(T); L = symb2input(L);
[~, J] = TL2AJ(T, L, N);
if isnumeric(J)
    odeFun = @(t, xn) J*xn;
else
    odeFun = @(t, xn) J(sum(xn(indxAS)))*xn;
end

% integration of ODEs
% -------------------
options = odeset('Events', @(t,y) myEvents(t,y,indxAS, 1e6), 'AbsTol', 1e-8, 'RelTol', 1e-8);

% integration of the dynamics
[time, nout] = ode45(odeFun, Tspan, [xn0; 0], options);
nout = nout(:,indxAS);
ntout = sum(nout,2);

% evaluate eigenvalues
if isnumeric(J)
    xeig = repmat(eig(J)', length(time), 1);
else
    xeig = NaN(length(time),3);
    for ii = 1:length(time)
        Jii = J(sum(ntout(ii)));
        xeig(ii,:) = eig(Jii(indxAS,indxAS));
    end
end

% minimum rate
if isnumeric(T)
    Tend = T;
else
    Tend = T(ntout(end));
end
if isnumeric(L)
    Lend = L;
else
    Lend = L(ntout(end));
end
minRate = min([Tend(:,1); Lend(:,1)]);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Built-in functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Tout, Lout] = fbDisruption_model(L, T, Lfb, Tfb)

% This function returns the input matrix T (cell transitions) and L (cell
% divisions) and the rates functions for the Feedback Disruption model.
% If a given constant rate (lambda1, omega21, omega13, omega14) is empty,
% feedback is applied as K/n or Kn. 
% 
% See notes_feedbackTestCase for more details.

% model definition
syms n real

% include feedback models
Lout = [];
for ii = 1:size(L, 1)
    if Lfb(ii,2) == 1
        pii = @(n) fb_k2n(Lfb(ii,1), 1, n);
    elseif Lfb(ii,2) == -1
        pii = @(n) fb_kn(Lfb(ii,1), 1, n);
    else
        pii = @(n) fb_none(L(ii,1));
    end
    Lout = [Lout; pii(n) L(ii,2:end)];
end
Tout = [];
for ii = 1:size(T, 1)
    if Tfb(ii,2) == 1
        pii = @(n) fb_k2n(Tfb(ii,1), 1, n);
    elseif Tfb(ii,2) == -1
        pii = @(n) fb_kn(Tfb(ii,1), 1, n);
    else
        pii = @(n) fb_none(T(ii,1));
    end
    Tout = [Tout; pii(n) T(ii,2:end)];
end

return

function input = symb2input(input_symb)

% convert input matrix from symbolic to function/double

if isnumeric(input_symb)
    input = input_symb;
else
    input = matlabFunction(input_symb);
    if nargin(input) == 0
        input = input();
    end
end

function r = fb_k2n(r0, neq, n)

r = r0*neq/n;

return

function r = fb_kn(r0, neq, n)

r = r0/neq*n;

return

function r = fb_none(r0)

r = r0;

return