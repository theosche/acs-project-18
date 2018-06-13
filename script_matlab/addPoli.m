%--------------------------------------------------------------------------
% synthax: function [p] = addPoli(p1,p2)
% 
% Description: add two polynomials
%
% Date: 13.06.18
%--------------------------------------------------------------------------
function [p] = addPoli(p1,p2)
p = [p1,zeros(1,max(length(p1),length(p2))-length(p1))]+...
    [p2,zeros(1,max(length(p1),length(p2))-length(p2))];
end