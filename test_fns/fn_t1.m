function f = fn_t1( x )
    f = x( 1 ) + x( 2 );
    g1 = 0.5 * sin( 2 * pi * ( x(1)^2 - 2*x(2) ) ) + x(1) + 2*x(2) - 1.5;
    g2 = -x(1)^2 - x(2)^2 + 1.5;
    f = [ f g1 g2 ];
end