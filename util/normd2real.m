function f = normd2real( x, bounds )
    f = x .* ( bounds( : , 2 ) - bounds( : , 1 ) ) + bounds( : , 1 );
end
