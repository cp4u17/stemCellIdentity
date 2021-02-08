function [value,isterminal,direction] = myEvents(~,y,indxAS, ntMax)
% stop when ntMax cells are reached (hardcoded param)
 value = sum(y(indxAS)) - ntMax;
 isterminal = 1;
 direction = 0;
return