function b_vector = get_b_vector(B, j_max, l_max)

b_vector = zeros(j_max+1, l_max);
for j = 0:j_max
    for l = 1:l_max
        b_vector(j+1, l) = fun_b(l/B^j, B);
    end
end

end