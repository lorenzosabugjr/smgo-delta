function f = fn_ariafar2( x_e )
    % from "ADMMBO: Bayesian optimization with unknown constraints using ADMM"   
    
    f = 0.5*sin(2*pi*(x_e(1)^2 - 2*x_e(2))) + x_e(1) + 2*x_e(2) - 1.5;
end