function out = opt_smgo_d( options )
% SMGO-Delta
% - uses global Lipschitz constants for f and g
% - uses SM function estimate of g for fulfillment/violation of the
%   constraint
% - includes multiple constraints

% Arguments to this function:
% f  - test function/experiment, accepts a function handle
% g  - constraint function/s, accepts an array of function handles, where array length is the number of constraints
% x0 - starting point (N x 1) column matrix, can be [].
% bounds - lower and upper bounds. must be an (N x 2) matrix
% max_iter - maximum iterations for a run

% (Optional) tunable optimization parameters:
% - alpha - expected improvement threshold
% - beta  - weighting between lower bound and central approximation
% - B     - number of candidate points generated across segments between candidate points
% - phi   - age multiplier

%% Acquiring function arguments/options
max_iter = getfield( options, 'maxiter' );
f        = getfield( options, 'objfun' );
bounds   = getfield( options, 'bounds' );
[ D , ~ ] = size( bounds );

if isfield( options, 'sobolsize' )
    sbl_size = options.sobolsize;
else
    sbl_size = 500;
end
sbl_seq = sobolset( D );
sbl_seq = sbl_seq( 1:sbl_size, : )';

if isfield( options, 'confun' )
    g = options.confun;
else
    g = {};
end

if isfield( options, 'conlen' )
    g_len = options.conlen;
else
    g_len = 0;
end

if isfield( options, 'startpt' )
    x0 = options.startpt;
else
    x0 = rand( D, 1 );
    disp( 'Starting from a random start point.' );
end

if isfield( options, 'alpha' )
    alpha = options.alpha;
else
    alpha = 0.0025;
end

if isfield( options, 'beta' )
    beta = options.beta;
else
    beta = 0.5;
end

if isfield( options, 'delta' )
    delta = options.delta;
else
    delta = 0.1;
end

if isfield( options, 'phi' )
    phi = options.phi;
else
    phi = 1e-5;
end

if isfield( options, 'B' )
    B = options.B;
else
    B = 5;
end

diam  = sqrt( D );   % diameter of search space
if ~exist( 'g_len', 'var' )
    g_len = length( g ); % number of constraints (g is an array of function handles)
end

%% Define initial variable values

% Lipschitz constants
fgam_0    = 1e-6;   ggam_0    = 1e-6 * ones( g_len, 1 );
fgam      = fgam_0; ggam      = ggam_0;
fgam_prev = fgam_0; ggam_prev = ggam_0;

% Noise estimations
feps_0    = 0;      geps_0    = zeros( g_len, 1 );
feps      = feps_0; geps      = geps_0;

MODE_EXPLOIT = 1;
MODE_EXPLORE = 2;

X_n = [];
mode_hist   = NaN( 1, max_iter );
mode_prev   = MODE_EXPLORE;
opt_z       = Inf;
opt_z_hist  = NaN( 1, max_iter );
opt_x       = NaN( D, 1 );
opt_x_hist  = NaN( D, max_iter );
calc_t_hist = NaN( 1, max_iter );
eps_hist    = NaN( 1 + g_len, max_iter );
gam_hist    = NaN( 1 + g_len, max_iter );

% trust region parameters
tr_coeff = 0.5;                tr_exp   = 0;
tr_hist  = NaN( 1, max_iter ); tr_exp_0 = 0;

% additional rows for X_n
XN_ROW_FVAL = D + 1;
XN_ROW_GVAL = D + 2;

%% Construct the database of segments and candidate points
% Database of Mode psi: from existing points, we calculate the candidate 
%   points (cdpts), and the value of uncertainty. Starting as empty
db_cdpt = [];

% Array structure:
% - rows 1 to D: the midpoint coordinate
% - row D+1: the value of upper bound cone vertex at that cdpt
% - row D+2: the value of upper bound cone height at that cdpt
% - row D+3: the value of lower bound cone vertex at that cdpt
% - row D+4: the value of lower bound cone height at that cdpt
% - row D+5: the value of uncertainty (lambda) at that cdpt
CDPT_ROW_F_UB_VTX = D + 1;
CDPT_ROW_F_UB_HGT = D + 2;
CDPT_ROW_F_LB_VTX = D + 3;
CDPT_ROW_F_LB_HGT = D + 4;
CDPT_ROW_F_LAMBDA = D + 5;

% The succeeding rows of Mode psi describe similar information, but for each
%   constraint. As all constraints are black-box too, we estimate the bounds
%   and the uncertainty (Pi_s)
CDPT_ROW_G_INFO = D + 6;
CDPT_G_UB_VTX   = 0;
CDPT_G_UB_HGT   = 1;
CDPT_G_LB_VTX   = 2;
CDPT_G_LB_HGT   = 3;
CDPT_G_LAMBDA   = 4;
CDPT_G_EST      = 5;

% Sobol sequence-based exploitation array declaration
SBL_ROW_FUB    = D + 1;
SBL_ROW_FLB    = D + 2;
SBL_ROW_G_INFO = D + 3;
SBL_ROW_GUB    = 0;
SBL_ROW_GLB    = 1;


%% loop to max iterations
for iter = 1:max_iter    
    if iter == 1
        x_n = x0;        
    end
       
    if isfield( options, 'confun' )
        z_n = f( normd2real( x_n, bounds ) );
        c_n = NaN(g_len, 1);
        for g_i = 1:g_len
            c_n( g_i ) = g{ g_i }( normd2real( x_n, bounds ) );
        end    
    else
        zc = f( normd2real( x_n, bounds ) ); 
        z_n = zc( 1 );
        c_n = NaN(g_len, 1);
        for g_i = 1:g_len
            c_n( g_i ) = zc( 1 + g_i );
        end    
    end
    X_n = [ X_n [ x_n; z_n; c_n ] ];          % stacking the column vector into collection of samples  
    
    calc_time = tic;
    if iter > 1        
        %% Calculating noise
        
        % collecting which points are within radius 'rad' from each other
        % increment 'rad' when there are no points inside the ball of any other
        for rad = 0.1:0.1:1.0
            ball_rad = rad * diam/2;
            ball_idx = zeros( iter );
            for k = 1:iter
                ball_idx( k, : ) = sum( ( X_n( 1:D, : ) - X_n( 1:D, k ) * ones( 1, iter ) ) .^ 2 ) <= ball_rad ^ 2;
                ball_idx( k, k ) = 0;
            end
            if sum( sum( ball_idx ) ) > 0
                break;
            end
        end
        ball_idx = logical( ball_idx );
        
        % collecting the noise values
        ns_lkup = NaN( 1, iter );
        for j = 1:iter
            if sum( ball_idx( j, : ) ) > 0
                ns_lkup( j ) = max( abs( X_n( XN_ROW_FVAL, ball_idx( j, : ) ) - X_n( XN_ROW_FVAL, j ) ) );
            end
        end
        feps = sum( ns_lkup( ns_lkup > 0 ) ) / ( 2 * sum( sum( ball_idx ) ) );
        
        for g_i = 1:g_len
            XN_ROW_GVAL_I = XN_ROW_GVAL + g_i - 1;            
            ns_lkup = NaN( 1, iter );
            if sum( ball_idx( j, : ) ) > 0
                ns_lkup( j ) = max( abs( X_n( XN_ROW_GVAL_I, ball_idx( j, : ) ) - X_n( XN_ROW_GVAL_I, j ) ) );
            end            
            geps( g_i ) = sum( ns_lkup( ns_lkup > 0 ) ) / ( 2 * sum( sum( ball_idx ) ) );
        end 
        
        %% Generate/update Lipschitz constants
        fgam_prev = fgam; 
        fgam_db = zeros( iter, 1 );
        for i = 1:iter
            elg_idx = abs( X_n( XN_ROW_FVAL, : ) - X_n( XN_ROW_FVAL, i ) ) > 2 * feps;            
            if sum( elg_idx ) > 0
                fgam_db( i ) = max( ( abs( X_n( XN_ROW_FVAL, elg_idx ) - X_n( XN_ROW_FVAL, i ) ) - 2*feps ) ./ ...
                                      sqrt( sum( ( X_n( 1:D, elg_idx ) - X_n( 1:D, i ) * ones( 1, sum( elg_idx ) ) ) .^ 2 ) ) );
            else 
                fgam_db( i ) = 1e-6;
            end            
        end           
        fgam = max( fgam_db );
        
        ggam_prev = ggam;
        for g_i = 1:g_len
            ggam_db = zeros( iter, 1 );
            XN_ROW_GVAL_I = XN_ROW_GVAL + g_i - 1;
            for i = 1:iter
                elg_idx = abs( X_n( XN_ROW_GVAL_I, : ) - X_n( XN_ROW_GVAL_I, i ) ) > 2 * geps( g_i );
                if sum( elg_idx ) > 0
                    ggam_db( i ) = max( ( abs( X_n( XN_ROW_GVAL_I, elg_idx ) - X_n( XN_ROW_GVAL_I, i ) ) - 2*geps( g_i ) ) ./ ...
                                          sqrt( sum( ( X_n( 1:D, elg_idx ) - X_n( 1:D, i ) * ones( 1, sum( elg_idx ) ) ) .^ 2 ) ) );
                else 
                    ggam_db( i ) = 1e-6;
                end
            end
            ggam( g_i ) = max( ggam_db );
        end
    end
    
    eps_hist( : , iter ) = [ feps; geps ];
    gam_hist( : , iter ) = [ fgam; ggam ];
    %% Updating the current best point
    if sum( c_n >= 0 ) == g_len && z_n < opt_z % build exploitation table only if all constraints are fulfilled
        opt_z = z_n;
        opt_x = x_n;        
        opt_z_new = true;
    else
        opt_z_new = false;
    end    
    opt_z_hist( iter )     = opt_z;
    opt_x_hist( : , iter ) = opt_x;
    
    %% iteratively update the midpoint database db_cdpt
    % - saves the vertex value of the upper (lower) bound
    % - also saves the hypercone height value
    % - only the height value is scaled
    % - if the new incoming upper (lower) bound is tighter than existing
    %   - replace the vertex value
    %   - replace the height value
    if iter ~= 1
        db_cdpt( [CDPT_ROW_F_UB_HGT CDPT_ROW_F_LB_HGT], : ) = ...
            ( fgam / fgam_prev ) .* db_cdpt( [CDPT_ROW_F_UB_HGT CDPT_ROW_F_LB_HGT], : );
        
        for g_i = 1:g_len
            CDPT_G_I = CDPT_ROW_G_INFO + 6*(g_i-1);
            db_cdpt( CDPT_G_I + [CDPT_G_UB_HGT CDPT_G_LB_HGT], : ) = ...
                ( ggam( g_i ) / ggam_prev( g_i ) ) .* db_cdpt( CDPT_G_I + [CDPT_G_UB_HGT CDPT_G_LB_HGT], : );
        end
        [ ~ , cdpts ] = size( db_cdpt );
        
        % calculating if i should update the f upper (lower) bounds
        new_cone_hgt = fgam * sqrt( sum( ( db_cdpt( 1:D, : ) - x_n * ones( 1, cdpts ) ) .^ 2 ) );
        new_f_ub     = z_n + feps + new_cone_hgt; 
        new_f_lb     = z_n - feps - new_cone_hgt;
        f_ub_update  = ( db_cdpt( CDPT_ROW_F_UB_VTX, : ) + db_cdpt( CDPT_ROW_F_UB_HGT, : ) + feps) > new_f_ub;
        f_lb_update  = ( db_cdpt( CDPT_ROW_F_LB_VTX, : ) - db_cdpt( CDPT_ROW_F_LB_HGT, : ) - feps) < new_f_lb;
        
        % replacing the f upper (lower) bounds data if needed
        db_cdpt( CDPT_ROW_F_UB_VTX, f_ub_update ) = z_n;
        db_cdpt( CDPT_ROW_F_UB_HGT, f_ub_update ) = new_cone_hgt( :, f_ub_update );
        db_cdpt( CDPT_ROW_F_LB_VTX, f_lb_update ) = z_n;
        db_cdpt( CDPT_ROW_F_LB_HGT, f_lb_update ) = new_cone_hgt( :, f_lb_update );
        
        db_cdpt( end-1, : ) = min( db_cdpt( end-1, : ), sqrt( sum( ( db_cdpt( 1:D, : ) - x_n * ones( 1, cdpts ) ) .^ 2 ) ) );
        % incrementing age of existing candidate points
        db_cdpt( end, : ) = db_cdpt( end, : ) + 1;
        
        for g_i = 1:g_len
            CDPT_G_I = CDPT_ROW_G_INFO + 6*(g_i-1);
            % calculating if i should update the g upper (lower) bounds
            new_cone_hgt = ggam( g_i ) * sqrt( sum( ( db_cdpt( 1:D, : ) - x_n * ones( 1, cdpts ) ) .^ 2 ) );
            new_g_ub     = c_n( g_i ) + geps( g_i ) + new_cone_hgt; 
            new_g_lb     = c_n( g_i ) - geps( g_i ) -  new_cone_hgt;
            g_ub_update  = new_g_ub < ( db_cdpt( CDPT_G_I + CDPT_G_UB_VTX, : ) + db_cdpt( CDPT_G_I + CDPT_G_UB_HGT, : ) ) + geps( g_i );
            g_lb_update  = new_g_lb > ( db_cdpt( CDPT_G_I + CDPT_G_LB_VTX, : ) - db_cdpt( CDPT_G_I + CDPT_G_LB_HGT, : ) ) - geps( g_i );

            % replacing the g upper (lower) bounds data if needed
            db_cdpt( CDPT_G_I + CDPT_G_UB_VTX, g_ub_update ) = c_n( g_i );
            db_cdpt( CDPT_G_I + CDPT_G_UB_HGT, g_ub_update ) = new_cone_hgt( :, g_ub_update );
            db_cdpt( CDPT_G_I + CDPT_G_LB_VTX, g_lb_update ) = c_n( g_i );
            db_cdpt( CDPT_G_I + CDPT_G_LB_HGT, g_lb_update ) = new_cone_hgt( :, g_lb_update );
        end
    end 
    
    %% introduce additional cdpts to db_cdpt
    % Generate candidate points: calculate end points from x_n in the cardinal directions
    db_cdpt_end1 = [ x_n * ones( 1, 2*D );
                     zeros( 5 + 6*g_len + 2, 2*D ) ];
    n = size( db_cdpt_end1, 1 );
    db_cdpt_end1( 1:(2*n+1):end )     = 0;
    db_cdpt_end1( (n+1):(2*n+1):end ) = 1;

    if iter > 1
        db_cdpt_end1 = [ db_cdpt_end1 [ X_n( 1:D, 1:end-1 ); zeros( 5 + 6*g_len + 2, iter-1 ) ] ];
    end

    % draw candidate points along the directions from x_n to the cdpt_ends
    if iter == 1
        db_cdpt_iter = [ sbl_seq; zeros( 5 + 6*g_len + 2, sbl_size ) ]; % TODO: use another Sobol sequence for exploration part
    else
        db_cdpt_iter = [];
    end
    for cdpt_w = ( 1:(B-1) ) / B
        db_cdpt_iter = [ db_cdpt_iter ( cdpt_w*db_cdpt_end1 + (1-cdpt_w)*( [ x_n; zeros( 5 + 6*g_len + 2, 1 ) ] ) ) ];
    end
        
    [ ~, new_cdpts ] = size( db_cdpt_iter );
    for new_cdpt = 1:new_cdpts
        % calculate the upper and lower bounds w.r.t. f
        % - *_ulb_lookup stores the locations of the points x_n, and the
        %   heights of the cones from x_n's
        f_ulb_lookup = [ X_n( XN_ROW_FVAL, 1:iter ); ...
                         fgam * sqrt( sum( ( X_n( 1:D, 1:iter ) - db_cdpt_iter( 1:D, new_cdpt ) * ones( 1, iter ) ) .^ 2 ) ) ];
        [ ~, ub_idx ] = min( f_ulb_lookup( 1, : ) + f_ulb_lookup( 2, : ) ); % calculating individual upper bound values
        [ ~, lb_idx ] = max( f_ulb_lookup( 1, : ) - f_ulb_lookup( 2, : ) ); % calculating individual lower bound values        
        db_cdpt_iter( CDPT_ROW_F_UB_VTX:CDPT_ROW_F_LB_HGT, new_cdpt ) = [ f_ulb_lookup( 1:2, ub_idx );
                                                                          f_ulb_lookup( 1:2, lb_idx ) ];
                                                                      
        % calculate the upper and lower bounds w.r.t. g
        for g_i = 1:g_len
            XN_ROW_GVAL_I = XN_ROW_GVAL + g_i - 1;
            CDPT_G_I      = CDPT_ROW_G_INFO + 6*(g_i-1);
            g_ulb_lookup = [ X_n( XN_ROW_GVAL_I, 1:iter ); ...
                             ggam( g_i ) * sqrt( sum( ( X_n( 1:D, 1:iter ) - db_cdpt_iter( 1:D, new_cdpt ) * ones( 1, iter ) ) .^ 2 ) ) ];
            [ ~, ub_idx ] = min( g_ulb_lookup( 1, : ) + g_ulb_lookup( 2, : ) ); % calculating individual upper bound values
            [ ~, lb_idx ] = max( g_ulb_lookup( 1, : ) - g_ulb_lookup( 2, : ) ); % calculating individual lower bound values        
            db_cdpt_iter( CDPT_G_I + (CDPT_G_UB_VTX:CDPT_G_LB_HGT), new_cdpt ) = [ g_ulb_lookup( 1:2, ub_idx );
                                                                                   g_ulb_lookup( 1:2, lb_idx ) ];
        end
        db_cdpt_iter( end-1, new_cdpt ) = min( sqrt( sum( ( X_n( 1:D, 1:end ) - db_cdpt_iter( 1:D, new_cdpt ) * ones( 1, iter ) ) .^ 2 ) ) );
    end
    db_cdpt = [ db_cdpt db_cdpt_iter ];
    
    db_cdpt( CDPT_ROW_F_LAMBDA, : ) = ( db_cdpt( CDPT_ROW_F_UB_VTX, : ) + db_cdpt( CDPT_ROW_F_UB_HGT, : ) ) ... % lambda is upper bound ...
                                    - ( db_cdpt( CDPT_ROW_F_LB_VTX, : ) - db_cdpt( CDPT_ROW_F_LB_HGT, : ) ) + 2 * feps;    % minus the lower bound    
    for g_i = 1:g_len
        CDPT_G_I = CDPT_ROW_G_INFO + 6*(g_i-1);
        db_cdpt( CDPT_G_I + CDPT_G_LAMBDA, : ) = ...
            ( db_cdpt( CDPT_G_I + CDPT_G_UB_VTX, : ) + db_cdpt( CDPT_G_I + CDPT_G_UB_HGT, : ) ) ... 
            - ( db_cdpt( CDPT_G_I + CDPT_G_LB_VTX, : ) - db_cdpt( CDPT_G_I + CDPT_G_LB_HGT, : ) ) + 2 * geps( g_i );
        db_cdpt( CDPT_G_I + CDPT_G_EST, : )    = ...
            ( ( db_cdpt( CDPT_G_I + CDPT_G_UB_VTX, : ) + db_cdpt( CDPT_G_I + CDPT_G_UB_HGT, : ) ) ... 
            + ( db_cdpt( CDPT_G_I + CDPT_G_LB_VTX, : ) - db_cdpt( CDPT_G_I + CDPT_G_LB_HGT, : ) ) ) / 2;
    end
    
    % generate/update trust region hyperbox
    if ( mode_prev == MODE_EXPLOIT ) && opt_z_new
        % if not already the largest measure, then enlarge the trust region hyperbox
        if tr_exp > 0
            tr_exp = tr_exp - 1;
        end
    else
        % shrink the trust region hyperbox
        if tr_exp < 10
            tr_exp = tr_exp + 1;
        else
            tr_exp = tr_exp_0;
        end
    end    
    tr_bounds = ones( D, 1 ) * [ 0 1 ] * ( tr_coeff ^ tr_exp ) + opt_x * ones( 1, 2 ) * ( 1 - tr_coeff ^ tr_exp );
    tr_hist( iter ) = tr_exp;
    
    % create candidate points inside tr_bounds, decided by location of sblset points
    % sblset points are D x num matrix
    sbl_cdpt = [ sbl_seq * ( tr_coeff ^ tr_exp ) + tr_bounds( : , 1 ) * ones( 1, sbl_size ); 
                 zeros( 2 + 2*g_len, sbl_size ) ];
    for sbl_i = 1:sbl_size
        sbl_cdpt( SBL_ROW_FUB, sbl_i ) = min( (X_n(XN_ROW_FVAL,:) + feps) + fgam*sqrt(sum((sbl_cdpt( 1:D, sbl_i )*ones(1,iter) - X_n(1:D,:)).^2)) );  
        sbl_cdpt( SBL_ROW_FLB, sbl_i ) = max( (X_n(XN_ROW_FVAL,:) - feps) - fgam*sqrt(sum((sbl_cdpt( 1:D, sbl_i )*ones(1,iter) - X_n(1:D,:)).^2)) );
        for g_i = 1:g_len
            SBL_G_I = SBL_ROW_G_INFO + 2*(g_i-1);
            sbl_cdpt( SBL_G_I + SBL_ROW_GUB, sbl_i ) = min( (X_n(XN_ROW_FVAL,:) + geps(g_i)) + ggam(g_i)*sqrt(sum((sbl_cdpt( 1:D, sbl_i )*ones(1,iter) - X_n(1:D,:)).^2)) );  
            sbl_cdpt( SBL_G_I + SBL_ROW_GLB, sbl_i ) = max( (X_n(XN_ROW_FVAL,:) - geps(g_i)) - ggam(g_i)*sqrt(sum((sbl_cdpt( 1:D, sbl_i )*ones(1,iter) - X_n(1:D,:)).^2)) );
        end
    end
    
    %% Algorithm (1) exploitation
    mode1_ok = 0;
    
    %% - Processing and choosing exploitation point from CDPT database
    [ ~ , cdpt_len ] = size( db_cdpt );
    cdpt_vld_idx     = true( 1, cdpt_len );
    % filtering due to constraints satisfaction
    for g_i = 1:g_len
        CDPT_G_I     = CDPT_ROW_G_INFO + 6*(g_i-1);
        cdpt_vld_idx = cdpt_vld_idx & ( db_cdpt( CDPT_G_I + CDPT_G_EST, : ) >= 0 );
    end    
    % filtering due to trust region
    cdpt_vld_idx = cdpt_vld_idx & prod( db_cdpt( 1:D, : ) > tr_bounds( : , 1 ) );
    cdpt_vld_idx = cdpt_vld_idx & prod( db_cdpt( 1:D, : ) < tr_bounds( : , 2 ) );
    
    cdpt_vld = db_cdpt( : , cdpt_vld_idx );
    if ~isempty( cdpt_vld )
        w_min = ( ( cdpt_vld( CDPT_ROW_F_UB_VTX, : ) + cdpt_vld( CDPT_ROW_F_UB_HGT, : ) ) + ...
                  ( cdpt_vld( CDPT_ROW_F_LB_VTX, : ) - cdpt_vld( CDPT_ROW_F_LB_HGT, : ) ) ) / 2;
        w_vld = -w_min + beta*cdpt_vld( CDPT_ROW_F_LAMBDA, : );
        [ ~, mdpt_idx ] = max( w_vld );

        % the 'while' snippet below is a last-resort way to avoid multiple samples 
        % at the same location, but practically should not run at all!
        cdpt_tmp = cdpt_vld( 1:D, mdpt_idx );
        % if chosen point x_n is too close to an existing sample, choose another
        while ~isempty( cdpt_tmp ) && min( sum( ( X_n( 1:D, : ) - cdpt_tmp * ones( 1, iter ) ) .^ 2 ) ) < 1e-9
    %         disp( 'Duplicate sampling detected. Choosing another exploitation point.' );
            cdpt_vld( : , mdpt_idx ) = [];
            w_vld( : , mdpt_idx )    = [];
            [ ~, mdpt_idx ]          = max( w_vld );        
            cdpt_tmp                 = cdpt_vld( 1:D, mdpt_idx );            
        end
    else
        cdpt_tmp = [];
    end
    
    %% Processing and choosing exploitation point from Sobol-generated candidate points
    sbl_vld_idx   = true( 1, sbl_size );
    % filtering due to constraints satisfaction
    for g_i = 1:g_len
        SBL_G_I = SBL_ROW_G_INFO + 2*(g_i-1);
        sbl_vld_idx = sbl_vld_idx & ( ( sbl_cdpt( SBL_G_I + SBL_ROW_GUB, : ) + sbl_cdpt( SBL_G_I + SBL_ROW_GLB, : ) ) >= 0 );
    end
    sbl_vld = sbl_cdpt( : , sbl_vld_idx );
    if ~isempty( sbl_vld )
        [ ~, mdpt_idx ] = max( -( sbl_vld( SBL_ROW_FUB, : ) + sbl_vld( SBL_ROW_FLB, : ) ) / 2 + beta*( sbl_vld( SBL_ROW_FUB, : ) - sbl_vld( SBL_ROW_FLB, : ) ) );

        % the 'while' snippet below is a last-resort way to avoid multiple samples
        % at the same location but practically should not run at all!
        sbl_tmp = sbl_vld( 1:D, mdpt_idx );
        % if chosen point x_n is too close to an existing sample, choose another
        while ~isempty( sbl_tmp ) && min( sum( ( X_n( 1:D, : ) - sbl_tmp * ones( 1, iter ) ) .^ 2 ) ) < 1e-9 
            disp( 'Duplicate sampling detected. Choosing another exploitation point.' );
            sbl_vld( : , mdpt_idx ) = [];
            [ ~, mdpt_idx ]         = min( sbl_vld( end, : ) );        
            sbl_tmp                 = sbl_vld( 1:D, mdpt_idx );            
        end
    else
        sbl_tmp = [];
    end
    
    % just selecting if the candidate point from db_cdpt or sbl_cdpt should be chosen
    if ~isempty( cdpt_tmp ) && ~isempty( sbl_tmp )
        if max( (X_n(XN_ROW_FVAL,:) - feps) - fgam*sqrt(sum((cdpt_tmp*ones(1,iter) - X_n(1:D,:)).^2)) ) < ...
           max( (X_n(XN_ROW_FVAL,:) - feps) - fgam*sqrt(sum((sbl_tmp*ones(1,iter) - X_n(1:D,:)).^2)) )
            x_n_tmp = cdpt_tmp;
        else
            x_n_tmp = sbl_tmp;
        end
    elseif ~isempty( cdpt_tmp )
        x_n_tmp = cdpt_tmp;        
    elseif ~isempty( sbl_tmp )
        x_n_tmp = sbl_tmp;
    else
        x_n_tmp = [];
    end
    
    if ~isempty( x_n_tmp )
        if opt_z < Inf && max( (X_n(XN_ROW_FVAL,:) - feps) - fgam*sqrt(sum((x_n_tmp*ones(1,iter) - X_n(1:D,:)).^2)) ) < opt_z - alpha * fgam
            x_n = x_n_tmp;
            mode_hist( iter ) = MODE_EXPLOIT;
            mode_prev         = MODE_EXPLOIT;
            mode1_ok          = true;
            disp( [ num2str( iter ), ' Exploitation' ] );
        end
    end
    
    %% Algorithm (2) exploration
    if ~mode1_ok
        [ ~ , cdpt_len ] = size( db_cdpt );
        cdpt_vld_idx     = true( 1, cdpt_len );
        w_bst            = ones( 1, cdpt_len );
        w_unc            = zeros( 1, cdpt_len );        
        w_viol           = ones( 1, cdpt_len );
        for g_i = 1:g_len
            CDPT_G_I     = CDPT_ROW_G_INFO + 6*(g_i-1);
            g_i_vld      = ( db_cdpt( CDPT_G_I + CDPT_G_EST, : ) >= 0 );
            w_bst        = w_bst .* (g_i_vld+1);
            w_unc        = w_unc + ( db_cdpt( CDPT_G_I + CDPT_G_LAMBDA, : ) )/ggam( g_i ); %  - 2 * geps( g_i )
            cdpt_vld_idx = cdpt_vld_idx & g_i_vld;
            
            XN_ROW_GVAL_I = XN_ROW_GVAL + g_i - 1;
            g_i_min       = min( X_n( XN_ROW_GVAL_I, : ) );
            if g_i_min >= 0
                continue;
            else                
                w_viol_idx = ( db_cdpt( CDPT_G_I + CDPT_G_EST, : ) < 0 );
                w_viol_div = 1 + w_viol_idx .* db_cdpt( CDPT_G_I + CDPT_G_EST, : ) / g_i_min;
                w_viol     = w_viol ./ w_viol_div;
            end            
        end
        w_vld = ( db_cdpt( CDPT_ROW_F_LAMBDA, : ) ) .* cdpt_vld_idx / fgam; %  - 2 * feps
        
        [ ~, mdpt_idx ] = max( (1-delta)*w_vld + delta*w_unc.*w_bst.*w_viol.*db_cdpt( end-1, : ) + phi.*db_cdpt( end, : ) );
        x_n = db_cdpt( 1:D, mdpt_idx );
        
        % the 'while' snippet below is a last-resort way to avoid multiple samples at the same location
        % but practically should not run at all!
        while ~isempty( x_n ) && min( sum( ( X_n( 1:D, : ) - x_n * ones( 1, iter ) ) .^ 2 ) ) < 1e-9 % if chosen point x_n is too close to an existing sample, choose another
%             disp( 'Duplicate sampling detected. Choosing another exploration point.' );
            db_cdpt( : , mdpt_idx ) = [];
            w_vld( : , mdpt_idx )   = [];
            w_unc( : , mdpt_idx )   = [];
            w_bst( : , mdpt_idx )   = [];
            w_viol( : , mdpt_idx )  = [];
            [ ~, mdpt_idx ] = max( (1-delta)*w_vld + delta*w_unc.*w_bst.*w_viol.*db_cdpt( end-1, : ) + phi.*db_cdpt( end, : ) );        
            x_n = db_cdpt( 1:D, mdpt_idx );
        end
        mode_hist(iter) = MODE_EXPLORE;
        mode_prev       = MODE_EXPLORE;
        disp( [ num2str( iter ), ' Exploration' ] );
    end

    % after i cleared the chosen column in db_cdpt (look at lines 101, 424), why do duplicates still exist?
    if sum( all( X_n( 1:D, : ) == x_n ) ) > 0 % this means that i'm selecting a cdpt even though it's already sampled!
        fprintf('Iteration %d: %d duplicate/s detected.\n', iter, sum( all( X_n( 1:D, : ) == x_n ) ) );
    end

    calc_t_hist( iter ) = toc( calc_time );
    x_hist = X_n( 1:D, : );    
    if ~exist( 'opt_x', 'var' )
        opt_x = NaN( D, 1 );
    end    
    
    %% Debug section        
    if ( fgam == Inf ) || ismember( Inf, ggam )
       iter; 
    end
end

z_hist = X_n( XN_ROW_FVAL, : );
c_hist = [];
for g_i = 1:g_len
    XN_ROW_GVAL_I = XN_ROW_GVAL + g_i - 1;
    c_hist = [ c_hist; X_n( XN_ROW_GVAL_I, : ) ];
end
    
out.opt_x       = opt_x;
out.opt_z       = opt_z;
out.x_hist      = x_hist;
out.z_hist      = z_hist;
out.c_hist      = c_hist;
out.opt_x_hist  = opt_x_hist;
out.opt_z_hist  = opt_z_hist;
out.calc_t_hist = calc_t_hist;
out.mode_hist   = mode_hist;
out.eps_hist    = eps_hist;
out.gam_hist    = gam_hist;
out.tr_hist     = tr_hist;