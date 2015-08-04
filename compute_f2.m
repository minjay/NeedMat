function f2_value = compute_f2( u )

neg_u = -abs(u);

f1 = @(x) exp(-1./(1-x.^2));

f2_value = integral(f1, -1, neg_u)/integral(f1, -1, 1);

%\psi(u)=1-\psi(-u) if u>0
if u>0
    f2_value = 1-f2_value;
end

end