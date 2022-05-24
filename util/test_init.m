if test_fn == FN_G04    
    options.objfun       = @( x )fn_g04( x );
    options.objfunb      = @( x )fn_g04b( x );
    options.conlen       = 6;
    options.bounds       = [ 78 102;
                             33 45;
                             27 45;
                             27 45;
                             27 45  ];
    opt_dim              = 5;
elseif test_fn == FN_G05    
    options.objfun       = @( x )fn_g05( x );
    options.objfunb      = @( x )fn_g05b( x );
    options.conlen       = 5;
    options.bounds       = [ 0 1200;
                             0 1200;
                             -0.55 0.55;
                             -0.55 0.55 ];
    opt_dim              = 4;
elseif test_fn == FN_G06    
    options.objfun       = @( x )fn_g06( x );
    options.objfunb      = @( x )fn_g06b( x );
    options.conlen       = 2;
    options.bounds       = [ 13 100;
                             0 100 ];
    opt_dim              = 2;
elseif test_fn == FN_G08    
    options.objfun       = @( x )fn_g08( x );
    options.objfunb      = @( x )fn_g08b( x );
    options.conlen       = 2;
    options.bounds       = [ 0 10;
                             0 10 ];
    opt_dim              = 2;
elseif test_fn == FN_G09  
    options.objfun       = @( x )fn_g09( x );
    options.objfunb      = @( x )fn_g09b( x );
    options.conlen       = 4;
    options.bounds       = ones( 7, 1 ) * [ -10 10 ];
    opt_dim              = 7;
elseif test_fn == FN_G10    
    options.objfun       = @( x )fn_g10( x );
    options.objfunb      = @( x )fn_g10b( x );
    options.conlen       = 6;
    options.bounds       = [ 100 10000;
                             1000 10000;
                             1000 10000;
                             10 1000;
                             10 1000;
                             10 1000;
                             10 1000;
                             10 1000 ];
    opt_dim              = 8;
elseif test_fn == FN_G12    
    options.objfun       = @( x )fn_g12( x );
    options.objfunb      = @( x )fn_g12b( x );
    options.conlen       = 1;
    options.bounds       = ones( 3, 1 ) * [ 0 9 ];
    opt_dim              = 3;
elseif test_fn == FN_G23    
    options.objfun       = @( x )fn_g23( x );
    options.objfunb      = @( x )fn_g23b( x );
    options.conlen       = 2;
    options.bounds       = [ 0 300;
                             0 300;
                             0 100;
                             0 200;
                             0 100;
                             0 300;
                             0 100;
                             0 200;
                             0.01 0.03 ];
    opt_dim              = 9;
elseif test_fn == FN_G24    
    options.objfun       = @( x )fn_g24( x );
    options.objfunb      = @( x )fn_g24b( x );
    options.conlen       = 2;
    options.bounds       = [ 0 3;
                             0 4 ];
    opt_dim              = 2;
elseif test_fn == FN_T1    
    options.objfun       = @( x )fn_t1( x );
    options.conlen       = 2;
    options.bounds       = [ 0 1;
                             0 1 ];
    opt_dim              = 2;
elseif test_fn == FN_STYB_10D
    options.conlen = 0;
    options.bounds = repmat( [ -5 5 ], 10, 1 ) ;
    options.objfun = @( c ) fn_styblinski( c );
    opt_dim              = 10;
end
