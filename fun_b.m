function bx = fun_b( x, B )
%FUN_B Evaluate the function b at x with the parameter B.

bx = sqrt( compute_f3(x/B, B)-compute_f3(x, B) );
