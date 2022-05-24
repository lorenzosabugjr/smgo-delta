function f = fn_g08( x )
    f = -( sin( 2 * pi * x( 1 ) ) ^ 3 * sin( 2 * pi * x( 2 ) ) ) / ( x( 1 ) ^ 3 * ( x( 1 ) + x( 2 ) ) );
    g1 = x( 1 ) ^ 2 - x( 2 ) + 1;
    g2 = 1 - x( 1 ) + ( x( 2 ) - 4 ) ^ 2;
    f = [ f -g1 -g2 ];
end