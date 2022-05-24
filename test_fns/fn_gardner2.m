function f = fn_gardner2( x_e )
    % from "ADMMBO: Bayesian optimization with unknown constraints using ADMM"   
    
    f = -sin(x_e(1))*sin(x_e(2)) - 0.95;    
end