function j_max = get_j_max(B, l_max)
%GET_J_MAX   Computes the maximal j.
%
%   j_max = get_j_max(B, l_max)
%
% Inputs:
%   B - the parameter
%   l_max - the maximal l
%
% Outputs:
%   j_max - the maximal j
%
% Author: Minjie Fan, 2015

j_max = 0;
while ceil(B^(j_max-1))<=l_max
    j_max = j_max+1;
end
j_max = j_max-1;
if ceil(B^(j_max-1))==B^(j_max-1) && ceil(B^(j_max-1))==l_max
    j_max = j_max-1;
end

end