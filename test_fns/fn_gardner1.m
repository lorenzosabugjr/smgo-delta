function f = fn_gardner1( x_e )
    % from "ADMMBO: Bayesian optimization with unknown constraints using ADMM"   
    
    f = sin(x_e(1)) + x_e(2);    
end