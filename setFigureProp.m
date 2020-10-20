function setFigureProp(hf)

% Change fig properties here

% set figure properties
set(hf, ...
    'defaultaxesfontname', 'Times New Roman',...
    'defaultTextInterpreter', 'LaTeX', ...
    'defaultaxesfontsize',20,...
    'defaultlinelinewidth',2,...
    ... % 'defaultlinecolor','k',...
    ... % 'DefaultAxesColorOrder',[0 0 0],...
    ... % 'DefaultAxesLineStyleOrder','-|--|:|-.',...
    'DefaultAxesBox','on',...
    'DefaultAxesXcolor', [0, 0, 0],...
    'DefaultAxesYcolor', [0, 0, 0],...
    'DefaultAxesZcolor', [0, 0, 0],...
    'DefaultAxesXGrid', 'on',...
    'DefaultAxesYGrid', 'on',...
    'DefaultAxesZGrid', 'on');