function [A, J, symFlag, Tf, Lf, T, L] = TL2AJ(T, L, N)
% This function buld the A and J matrices starting from T and L.
% Assumptions:
% - in T only state transitions, no repeated
% - in L division in two states
%
% T     IN: Transition options. The first column contains the rates, the
%           second/third the incoming/outcoming states (M,3)
% L     IN: Division options. The first column contains the rates, the
%           second the incoming states and the third and fourth the outcoming states (L,4)
% N     IN: total number of states
% A     OUT: Adjacency (N,N)
% J     OUT: Jacobian matrix (N,N)
% symFlag OUT:
% Tf    OUT:
% Lf    OUT:
% T     OUT:
% L     OUT:
% 
% Cristina Parigini, 14/03/2020
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

% initialize variables in case of constant/not constant rates
if ~isnumeric(T) || ~isnumeric(L)
    syms n real
    symFlag = true;
else
    symFlag = false;
end
if ~isnumeric(T)
    Tf = T(n); Tf(:,2:end) = []; 
    T = T(NaN); % indexing does not depend on input
    At = sym('At', [N,N]); At(:) = 0; % initialization
else
    if ~isempty(T)
        Tf = T(:,1);
    else
        Tf = [];
    end
    At = zeros(N,N);
end
if ~isnumeric(L)
    Lf = L(n); Lf(:,2:end) = [];
    L = L(NaN); % indexing does not depend on input
    Al = sym('Al', [N,N]); Al(:) = 0; % initialization
    Dl = sym('Dl', [N,N]); Dl(:) = 0; % initialization
else
    if ~isempty(L)
        Lf = L(:,1);
    else
        Lf = [];
    end
    Al = zeros(N,N);
    Dl = zeros(N,N);
end

% Adjacency matrix
if ~isempty(T)
    At(sub2ind([N N], T(:,3), T(:,2))) = Tf;
end
if ~isempty(L)
    for il = 1:size(L,1) % loop in case of repetition (e.g. 1->1+1 and 1->1+2)
        indx = sub2ind([N N], L(il,3), L(il,2));
        Al(indx) = Al(indx) + Lf(il);
        indx = sub2ind([N N], L(il,4), L(il,2));
        Al(indx) = Al(indx) + Lf(il);
    end
end
% sum of transition and division parts
A = At + Al;

% Jacobian
Dt = diag(sum(At));
if ~isempty(L)
    if length(unique(L(:,2))) == length(L(:,2))
        Dl(sub2ind([N N], L(:,2), L(:,2))) = Lf;
    else
        for ii = 1:length(L(:,2))
            indx = sub2ind([N N], L(ii,2), L(ii,2));
            Dl(indx) = Dl(indx) + Lf(ii);
        end
    end
end
D = Dt + Dl;
J = A-D;

% convert to a matlab function in case of feedback
if symFlag
    A = matlabFunction(A);
    J = matlabFunction(J);
    if isnumeric(Tf)
        Tf = @(x) Tf();
    else
        Tf = matlabFunction(Tf);
    end
    if isnumeric(Lf)
        Lf = @(x) Lf();
    else
        Lf = matlabFunction(Lf);
    end
end

end