function f = fn_salomon( x_e )
    % Salomon function
    % http://benchmarkfcns.xyz/benchmarkfcns/salomonfcn.html
    % f(x*) = 0 at x* = (0, 0, ...)    

    if istable( x_e )
        x = table2array( x_e );
    else
        x = x_e;
    end
       
    f = 1 - cos( 2 * pi * sqrt( sum( x.^2 ) ) ) + 0.1 * sqrt( sum( x.^ 2 ) );
       
end