function f3_value = get_f3(x, B)

if x<0
    disp('x is not in the domain of f3!')
elseif x<=1/B
    f3_value = 1;
elseif x<=1
    f3_value = get_f2(1-2*B/(B-1)*(x-1/B));
else
    f3_value = 0;
end

end