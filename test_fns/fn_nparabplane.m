function f = fn_nparabplane( x_e )
    % simple parabolic function 

    if istable( x_e )
        x = table2array( x_e );
    else
        x = x_e;
    end

    f = -max( 3 - sum( ( x + 1.5 ) .^ 2 ), sum( x ) );
end