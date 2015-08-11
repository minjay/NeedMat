function j_max = get_j_max(B, l_max)

j_max = 0;
while ceil(B^(j_max-1))<=l_max
    j_max = j_max+1;
end
j_max = j_max-1;
if ceil(B^(j_max-1))==B^(j_max-1)
    j_max = j_max-1;
end

end