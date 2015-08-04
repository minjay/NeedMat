function f3x = compute_f3( x, B )
%F3 Evaluate the function f3 at x with the parameter B.

if x<0
    disp('x is not in the domain of f3!')
elseif x<=1/B
    f3x = 1;
elseif x<=1
    f3x = compute_f2(1-2*B/(B-1)*(x-1/B));
else
    f3x = 0;
end

end