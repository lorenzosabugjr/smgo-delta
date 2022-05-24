function f = fn_ariafar3( x_e )
    % from "ADMMBO: Bayesian optimization with unknown constraints using ADMM"   
    
    f = -x_e(1)^2 - x_e(2)^2 + 1.5;
end