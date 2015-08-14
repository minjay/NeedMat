function Nside = get_Nside(B, j)
%GET_NSIDE   Computes Nside.
%
%   Nside = get_Nside(B, j)
% Inputs:
%   B - the parameter
%   j - the frequency
%
% Outputs:
%   Nside - the Nside
%
% Author: Minjie Fan, 2015

Nside = 1;
lb = floor(B^(j+1));
while 2*Nside<lb
    Nside = Nside*2;
end

end