function bl_vector = get_bl_vector(B, j_max, l_max)

bl_vector = zeros(j_max+1, l_max);
for j = 0:j_max
    for l = 1:l_max
        bl_vector(j+1, l) = fun_b(l/B^j, B)*(2*l+1)/(4*pi);
    end
end

end