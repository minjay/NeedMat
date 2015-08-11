function Nside = get_Nside(B, j)

Nside = 1;
lb = floor(B^(j+1));
while 2*Nside<lb
    Nside = Nside*2;
end

end