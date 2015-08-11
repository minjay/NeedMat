function b_value = fun_b( x, B )


b_value =  get_f3(x/B, B)-get_f3(x, B);
if b_value<0
    b_value = 0;
end
b_value = sqrt(b_value);

end
