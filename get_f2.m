function f2_value = get_f2( u )

f1 = @(x) exp(-1./(1-x.^2));

f2_value = integral(f1, -1, u)/integral(f1, -1, 1);

end