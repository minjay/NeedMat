function psi = spneedlet_eval_fast(B, j, bl_vector, P, dist, sqrt_lambda)

l_min = ceil(B^(j-1));
l_max = floor(B^(j+1));

n = length(dist);

psi = zeros(n, 1);

for i = 1:n
    psi(i) = sum(bl_vector(j+1, l_min:l_max).*P(i, l_min+1:l_max+1));
end

psi = psi*sqrt_lambda;

end