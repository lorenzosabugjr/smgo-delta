function out = opt_smgo_d(options)
% SMGO-Delta
% - uses global Lipschitz constants for f and g
% - uses SM function estimate of g for fulfillment/violation of the
%   constraint
% - includes multiple constraints

% Arguments to this function (contained in struct 'options'):
%   .objfun  - test function/experiment, accepts a function handle
%   .conlen  - [used if objfun returns the values of the constraint functions together with the objective, in a matrix] 
%                the number of constraints returned with the objective 
%   .confun  - [used if the constraints values are provided by separate constraint function/s]
%                an array of function handles, where array length is the number of constraints
%   .startpt - starting point (N x 1) column matrix, can be [].
%   .bounds  - lower and upper bounds. must be an (N x 2) matrix
%   .maxiter - maximum iterations for a run

% (Optional arguments) tunable optimization parameters (also in struct 'options'):
%   .alpha - expected improvement threshold
%   .beta  - weighting between lower bound and central approximation
%   .B     - number of candidate points generated across segments between candidate points
%   .phi   - age multiplier

% Outputs of this function are contained in the the struct 'out':
%   .opt_x - best (feasible) sampled point
%   .opt_z - objective value at opt_x
%   .opt_c - constraint values at opt_x
%   .opt_x_hist - history of best sampled points throughout the optimization process
%   .opt_z_hist - history of best objective (corresponding to opt_x_hist)
%   .opt_c_hist - history of constraints at the best point (corresponding to opt_x_hist)
%   .x_hist - history of all sampled points
%   .z_hist - history of all sampled objectives
%   .c_hist - history of all sampled constraints
%   .calc_t_hist - history of time elapsed to calculate the next point to sample
%   .mode_hist - history of SMGO-D mode taken (exploitation or exploration)
%   .eps_hist - history of calculated noise bounds
%   .gam_hist - history of calculated Lipschitz constants
%   .tr_hist - history of trust region exponents (higher exponent means smaller trust region size)

%% Acquiring function arguments/options
max_iter  = options.maxiter;
f         = options.objfun;
bounds    = options.bounds;
D         = size(bounds, 1);

% ===========================================
% quasi-random points distribution parameters
% ===========================================
if isfield(options, 'sobolsize')
    sbl_size = options.sobolsize;
else
    sbl_size = 500;
end
if isfield(options, 'sobol')
    if ~options.sobol
        sbl_size = 0;
    end
end
sbl_seq = sobolset(D);
sbl_seq = sbl_seq(1:sbl_size, :)';
sbl_seq(:,2) = [];
sbl_size = sbl_size - 1;
% ===========================================

% ===========================================
% trust region parameters
% ===========================================
tr_size = 0.1;
if isfield(options, 'trustregion')
    if ~options.trustregion
        tr_size = 1.0;
    end
end
tr_coeff = 0.5;
% ===========================================

% ===========================================
% constraint function/s parameters
% ===========================================
if isfield(options, 'confun')
    g = options.confun;
else
    g = {};
end

if isfield(options, 'conlen')
    g_len = options.conlen;
else
    g_len = 0;
end

if ~exist('g_len', 'var')
    g_len = length(g); % number of constraints (g is an array of function handles)
end
% ===========================================

% ===========================================
% starting point
% ===========================================
if isfield(options, 'startpt')
    x0 = options.startpt;
else
    x0 = rand(D, 1);
    disp('Starting from a random start point.');
end
% ===========================================

% ===========================================
% SMGO-D advanced parameters
% ===========================================
if isfield(options, 'alpha')
    alpha = options.alpha;
else
    alpha = 0.005;
end

if isfield(options, 'beta')
    beta = options.beta;
else
    beta = 0.1;
end

if isfield(options, 'delta')
    delta = options.delta;
else
    delta = 0.2;
end

if isfield(options, 'phi')
    phi = options.phi;
else
    phi = 1e-9;
end

if isfield(options, 'B')
    B = options.B;
else
    B = 5;
end
% ===========================================


%% Define initial variable values

diam  = sqrt(D);   % diameter of search space

% Lipschitz constants
fgam_0    = 1e-6;   ggam_0    = 1e-6 * ones(g_len, 1);
fgam      = fgam_0; ggam      = ggam_0;
fgam_prev = fgam_0; ggam_prev = ggam_0;

% Noise estimations
feps_0    = 0;      geps_0    = zeros(g_len, 1);
feps      = feps_0; geps      = geps_0;

X_n = [];

MODE_EXPLOIT = 1;
MODE_EXPLORE = 2;
mode_hist   = NaN(1, max_iter);
mode_prev   = MODE_EXPLOIT;

opt_x       = NaN(D, 1);            opt_x_hist  = NaN(D, max_iter);
opt_z       = Inf;                  opt_z_hist  = NaN(1, max_iter);
opt_c       = NaN(g_len, 1);        opt_c_hist  = NaN(g_len, max_iter);

calc_t_hist = NaN(1, max_iter);
eps_hist    = NaN(1 + g_len, max_iter);
gam_hist    = NaN(1 + g_len, max_iter);

tr_exp   = 0;
tr_exp_0 = 0;
tr_hist  = NaN(1, max_iter); 
mode1_thr = Inf;

% additional rows for X_n
XN_ROW_FVAL = D + 1;
XN_ROW_GVAL = D + 2;

%% Construct the database of segments and candidate points
% Database of Mode psi: from existing points, we calculate the candidate 
%   points (cdpts), and the value of uncertainty. Starting as empty
db_cdpt = [];

% Array structure:
% - rows 1 to D: the cdpoint coordinate
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

INFTY = 1e25;

%% ============================================================
% Main SMGO-?? loop
% =============================================================
% (1) Evaluation of (black-box) objective and constraint functions
% (2) Calculation of noise bounds and Lipschitz constants
% (3) Set membership (SM) model building (the bulk part of the algorithm in terms of memory and calc)
% (3.1) (Iterative) update of SM bounds information on existing candidate points
% (3.2) Generation of new candidate points and calculating their SM bounds
% (4) Exploitation routine (simply select a candidate point, using the updated SM bounds)
% (5) Exploration routine (similarly just a selection process of candidate point)

for iter = 1:max_iter    
    if iter == 1
        x_n = x0;        
    end
       
    if isfield(options, 'confun')
        z_n = f(normd2real(x_n, bounds));
        c_n = NaN(g_len, 1);
        for g_i = 1:g_len
            c_n(g_i) = g{g_i}(normd2real(x_n, bounds));
        end    
    else
        zc = f(normd2real(x_n, bounds)); 
        z_n = zc(1);
        c_n = NaN(g_len, 1);
        for g_i = 1:g_len
            c_n(g_i) = zc(1 + g_i);
        end    
    end
    X_n = [X_n [x_n; z_n; c_n]];          % side-stacking the column vector to collection of samples (fat matrix)
    

    calc_time = tic;
    if iter > 1        
        %% ================================
        %  Calculating noise
        %  ================================

        % collecting which points are within radius 'rad' from each other
        % increment 'rad' when there are no points inside the ball of any other
        for rad = 0.1:0.1:1.0
            ball_rad = rad * diam/2;
            ball_idx = zeros(iter);
            for k = 1:iter
                ball_idx(k, :) = vecnorm(X_n(1:D, :) - repmat(X_n(1:D, k), 1, iter)) <= ball_rad;
                ball_idx(k, k) = 0;
            end
            if sum(ball_idx,'all') > 0
                break;
            end
        end
        ball_idx = logical(ball_idx);
        
        % collecting the noise values
        ns_lkup = NaN(1, iter);
        for j = 1:iter
            if sum(ball_idx(j, :)) > 0
                ns_lkup(j) = max(abs(X_n(XN_ROW_FVAL, ball_idx(j, :)) - X_n(XN_ROW_FVAL, j)));
            end
        end
        feps = sum(ns_lkup(ns_lkup > 0)) / (2 * sum(ball_idx,'all'));
        
        for g_i = 1:g_len
            XN_ROW_GVAL_I = XN_ROW_GVAL + g_i - 1;            
            ns_lkup = NaN(1, iter);
            for j = 1:iter
                if sum(ball_idx(j, :)) > 0
                    ns_lkup(j) = max(abs(X_n(XN_ROW_GVAL_I, ball_idx(j, :)) - X_n(XN_ROW_GVAL_I, j)));
                end            
            end
            geps(g_i) = sum(ns_lkup(ns_lkup > 0)) / (2 * sum(ball_idx,'all'));
        end 
        
        %% ========================================
        %  Generate/update Lipschitz constants
        %  ========================================
        fgam_prev = fgam; 
        fgam_db = zeros(iter, 1);
        for i = 1:iter
            elg_idx = abs(X_n(XN_ROW_FVAL, :) - X_n(XN_ROW_FVAL, i)) > 2 * feps;            
            if sum(elg_idx) > 0
                fgam_db(i) = max((abs(X_n(XN_ROW_FVAL, elg_idx) - X_n(XN_ROW_FVAL, i)) - 2*feps) ./ ...
                                  vecnorm(X_n(1:D, elg_idx) - repmat(X_n(1:D, i), 1, sum(elg_idx))));
            else 
                fgam_db(i) = 1e-6;
            end            
        end
        fgam = max(fgam_db);
        
        ggam_prev = ggam;
        for g_i = 1:g_len
            ggam_db = zeros(iter, 1);
            XN_ROW_GVAL_I = XN_ROW_GVAL + g_i - 1;
            for i = 1:iter
                elg_idx = abs(X_n(XN_ROW_GVAL_I, :) - X_n(XN_ROW_GVAL_I, i)) > 2 * geps(g_i);
                if sum(elg_idx) > 0
                    ggam_db(i) = max((abs(X_n(XN_ROW_GVAL_I, elg_idx) - X_n(XN_ROW_GVAL_I, i)) - 2*geps(g_i)) ./ ...
                                      vecnorm(X_n(1:D, elg_idx) - repmat(X_n(1:D, i), 1, sum(elg_idx))));
                else 
                    ggam_db(i) = 1e-6;
                end
            end
            ggam(g_i) = max(ggam_db);
        end
    end    
    eps_hist(:, iter) = [feps; geps];
    gam_hist(:, iter) = [fgam; ggam];

    %% Updating the current best point
    if sum(c_n >= 0) == g_len && z_n < opt_z % build exploitation table only if all constraints are fulfilled
        opt_z = z_n;
        opt_x = x_n;        
        opt_c = c_n;
        opt_z_new = true;
    else
        opt_z_new = false;
    end    
    
    %% ===============================================
    % Iteratively update the midpoint database db_cdpt
    % ================================================
    % - saves the vertex value of the upper (lower) bound
    % - also saves the hypercone height value
    % - only the height value is scaled
    % - if the new incoming upper (lower) bound is tighter than existing
    %   - replace the vertex value
    %   - replace the height value
    if iter ~= 1
        db_cdpt([CDPT_ROW_F_UB_HGT CDPT_ROW_F_LB_HGT], :) = ...
            (fgam / fgam_prev) .* db_cdpt([CDPT_ROW_F_UB_HGT CDPT_ROW_F_LB_HGT], :);
        
        for g_i = 1:g_len
            CDPT_G_I = CDPT_ROW_G_INFO + 6*(g_i-1);
            db_cdpt(CDPT_G_I + [CDPT_G_UB_HGT CDPT_G_LB_HGT], :) = ...
                (ggam(g_i) / ggam_prev(g_i)) .* db_cdpt(CDPT_G_I + [CDPT_G_UB_HGT CDPT_G_LB_HGT], :);
        end
        cdpts = size(db_cdpt, 2);
        
        % calculating if i should update the f upper (lower) bounds
        new_cone_hgt = fgam * vecnorm(db_cdpt(1:D, :) - repmat(x_n, 1, cdpts));
        new_f_ub     = z_n + feps + new_cone_hgt; 
        new_f_lb     = z_n - feps - new_cone_hgt;
        f_ub_update  = (db_cdpt(CDPT_ROW_F_UB_VTX, :) + db_cdpt(CDPT_ROW_F_UB_HGT, :) + feps) > new_f_ub;
        f_lb_update  = (db_cdpt(CDPT_ROW_F_LB_VTX, :) - db_cdpt(CDPT_ROW_F_LB_HGT, :) - feps) < new_f_lb;
        
        % replacing the f upper (lower) bounds data if needed
        db_cdpt(CDPT_ROW_F_UB_VTX, f_ub_update) = z_n;
        db_cdpt(CDPT_ROW_F_UB_HGT, f_ub_update) = new_cone_hgt(:, f_ub_update);
        db_cdpt(CDPT_ROW_F_LB_VTX, f_lb_update) = z_n;
        db_cdpt(CDPT_ROW_F_LB_HGT, f_lb_update) = new_cone_hgt(:, f_lb_update);
        
        db_cdpt(end-1, :) = min(db_cdpt(end-1, :), vecnorm(db_cdpt(1:D, :) - repmat(x_n, 1, cdpts)) / diam);
        % incrementing age of existing candidate points
        db_cdpt(end, :) = db_cdpt(end, :) + 1;
        
        for g_i = 1:g_len
            CDPT_G_I = CDPT_ROW_G_INFO + 6*(g_i-1);
            % calculating if i should update the g upper (lower) bounds
            new_cone_hgt = ggam(g_i) * vecnorm(db_cdpt(1:D, :) - repmat(x_n, 1, cdpts));
            new_g_ub     = c_n(g_i) + geps(g_i) + new_cone_hgt; 
            new_g_lb     = c_n(g_i) - geps(g_i) -  new_cone_hgt;
            g_ub_update  = new_g_ub < (db_cdpt(CDPT_G_I + CDPT_G_UB_VTX, :) + db_cdpt(CDPT_G_I + CDPT_G_UB_HGT, :)) + geps(g_i);
            g_lb_update  = new_g_lb > (db_cdpt(CDPT_G_I + CDPT_G_LB_VTX, :) - db_cdpt(CDPT_G_I + CDPT_G_LB_HGT, :)) - geps(g_i);

            % replacing the g upper (lower) bounds data if needed
            db_cdpt(CDPT_G_I + CDPT_G_UB_VTX, g_ub_update) = c_n(g_i);
            db_cdpt(CDPT_G_I + CDPT_G_UB_HGT, g_ub_update) = new_cone_hgt(:, g_ub_update);
            db_cdpt(CDPT_G_I + CDPT_G_LB_VTX, g_lb_update) = c_n(g_i);
            db_cdpt(CDPT_G_I + CDPT_G_LB_HGT, g_lb_update) = new_cone_hgt(:, g_lb_update);
        end
    end 
    
    %% introduce additional cdpts to db_cdpt
    % Generate candidate points: calculate end points from x_n in the cardinal directions
    db_cdpt_end1 = [x_n * ones(1, 2*D);
                     zeros(5 + 6*g_len + 2, 2*D)];
    n = size(db_cdpt_end1, 1);
    db_cdpt_end1(1:(2*n+1):end)     = 0;
    db_cdpt_end1((n+1):(2*n+1):end) = 1;

    if iter > 1
        db_cdpt_end1 = [db_cdpt_end1 [X_n(1:D, 1:end-1); zeros(5 + 6*g_len + 2, iter-1)]];
    end

    % draw candidate points along the directions from x_n to the cdpt_ends
    if iter == 1
        db_cdpt_iter = [sbl_seq; zeros(5 + 6*g_len + 2, sbl_size)];
    else
        db_cdpt_iter = [];
    end
    for cdpt_w = (1:(B-1)) / B
        db_cdpt_iter = [db_cdpt_iter (cdpt_w*db_cdpt_end1 + (1-cdpt_w)*([x_n; zeros(5 + 6*g_len + 2, 1)]))];
    end
        
    [~, new_cdpts] = size(db_cdpt_iter);
    for new_cdpt = 1:new_cdpts
        % calculate the upper and lower bounds w.r.t. f
        % - *_ulb_lookup stores the locations of the points x_n, and the
        %   heights of the cones from x_n's
        f_ulb_lookup = [X_n(XN_ROW_FVAL, 1:iter); ...
                         fgam * vecnorm(X_n(1:D, 1:iter) - repmat(db_cdpt_iter(1:D, new_cdpt), 1, iter))];
        [~, ub_idx] = min(f_ulb_lookup(1, :) + f_ulb_lookup(2, :)); % calculating individual upper bound values
        [~, lb_idx] = max(f_ulb_lookup(1, :) - f_ulb_lookup(2, :)); % calculating individual lower bound values        
        db_cdpt_iter(CDPT_ROW_F_UB_VTX:CDPT_ROW_F_LB_HGT, new_cdpt) = [f_ulb_lookup(1:2, ub_idx);
                                                                       f_ulb_lookup(1:2, lb_idx)];
                                                                      
        % calculate the upper and lower bounds w.r.t. g
        for g_i = 1:g_len
            XN_ROW_GVAL_I = XN_ROW_GVAL + g_i - 1;
            CDPT_G_I      = CDPT_ROW_G_INFO + 6*(g_i-1);
            g_ulb_lookup = [X_n(XN_ROW_GVAL_I, 1:iter); ...
                             ggam(g_i) * vecnorm(X_n(1:D, 1:iter) - repmat(db_cdpt_iter(1:D, new_cdpt), 1, iter))];
            [~, ub_idx] = min(g_ulb_lookup(1, :) + g_ulb_lookup(2, :)); % calculating individual upper bound values
            [~, lb_idx] = max(g_ulb_lookup(1, :) - g_ulb_lookup(2, :)); % calculating individual lower bound values        
            db_cdpt_iter(CDPT_G_I + (CDPT_G_UB_VTX:CDPT_G_LB_HGT), new_cdpt) = [g_ulb_lookup(1:2, ub_idx);
                                                                                g_ulb_lookup(1:2, lb_idx)];
        end
        db_cdpt_iter(end-1, new_cdpt) = min(vecnorm(X_n(1:D, 1:end) - repmat(db_cdpt_iter(1:D, new_cdpt), 1, iter))) / diam;
    end
    db_cdpt = [db_cdpt db_cdpt_iter];
    
    db_cdpt(CDPT_ROW_F_LAMBDA, :) = (db_cdpt(CDPT_ROW_F_UB_VTX, :) + db_cdpt(CDPT_ROW_F_UB_HGT, :)) ... % lambda is upper bound ...
                                  - (db_cdpt(CDPT_ROW_F_LB_VTX, :) - db_cdpt(CDPT_ROW_F_LB_HGT, :)) + 2 * feps;    % minus the lower bound    
    for g_i = 1:g_len
        CDPT_G_I = CDPT_ROW_G_INFO + 6*(g_i-1);
        db_cdpt(CDPT_G_I + CDPT_G_LAMBDA, :) = ...
            (db_cdpt(CDPT_G_I + CDPT_G_UB_VTX, :) + db_cdpt(CDPT_G_I + CDPT_G_UB_HGT, :)) ... 
          - (db_cdpt(CDPT_G_I + CDPT_G_LB_VTX, :) - db_cdpt(CDPT_G_I + CDPT_G_LB_HGT, :)) + 2 * geps(g_i);
        db_cdpt(CDPT_G_I + CDPT_G_EST, :)    = ...
            delta * ((db_cdpt(CDPT_G_I + CDPT_G_UB_VTX, :) + db_cdpt(CDPT_G_I + CDPT_G_UB_HGT, :)) ... 
                   + (db_cdpt(CDPT_G_I + CDPT_G_LB_VTX, :) - db_cdpt(CDPT_G_I + CDPT_G_LB_HGT, :))) / 2 + ...
            (1 - delta) * (db_cdpt(CDPT_G_I + CDPT_G_LB_VTX, :) - db_cdpt(CDPT_G_I + CDPT_G_LB_HGT, :));
    end
    
    % generate/update trust region hyperbox
    if iter > 1
        if (mode_prev == MODE_EXPLOIT) && opt_z_new
            if z_n < mode1_thr
                % if the actual sample is actually less than expected improvement, 
                % then enlarge the trust region hyperbox
                if tr_exp > 0
                    tr_exp = tr_exp - 1;
                end
            end
        else
            % shrink the trust region hyperbox
            if tr_exp < 10
                tr_exp = tr_exp + 1;
            else
                tr_exp = tr_exp_0;
            end
        end    
    end
    tr_bounds = ones(D, 1) * [-0.5 0.5] * tr_size * (tr_coeff ^ tr_exp) + opt_x * ones(1, 2);
    tr_bounds = max(min(tr_bounds, 1.0), 0.0);
    tr_hist(iter) = tr_exp;
    
    % create candidate points inside tr_bounds, decided by location of sblset points
    % sblset points are D x num matrix    
    sbl_cdpt = [sbl_seq * tr_size * (tr_coeff ^ tr_exp) + tr_bounds(:, 1) * ones(1, sbl_size); 
                 zeros(2 + 2*g_len, sbl_size)];
    for sbl_i = 1:sbl_size
        sbl_dst = vecnorm(X_n(1:D,:) - repmat(sbl_cdpt(1:D, sbl_i), 1, iter));
        sbl_cdpt(SBL_ROW_FUB, sbl_i) = min(X_n(XN_ROW_FVAL,:) + feps + fgam*sbl_dst);  
        sbl_cdpt(SBL_ROW_FLB, sbl_i) = max(X_n(XN_ROW_FVAL,:) - feps - fgam*sbl_dst);
        for g_i = 1:g_len
            SBL_G_I = SBL_ROW_G_INFO + 2*(g_i-1);
            XN_ROW_GVAL_I = XN_ROW_GVAL + g_i - 1;
            sbl_cdpt(SBL_G_I + SBL_ROW_GUB, sbl_i) = min(X_n(XN_ROW_GVAL_I,:) + geps(g_i) + ggam(g_i)*sbl_dst);  
            sbl_cdpt(SBL_G_I + SBL_ROW_GLB, sbl_i) = max(X_n(XN_ROW_GVAL_I,:) - geps(g_i) - ggam(g_i)*sbl_dst);
        end
    end
    
    
    %% Algorithm (1) exploitation
    mode1_ok = 0;
    
    %% - Processing and choosing exploitation point from CDPT database
    [~ , cdpt_len] = size(db_cdpt);
    cdpt_vld_idx     = true(1, cdpt_len);
    % filtering due to constraints satisfaction
    for g_i = 1:g_len
        CDPT_G_I     = CDPT_ROW_G_INFO + 6*(g_i-1);
        cdpt_vld_idx = cdpt_vld_idx & (db_cdpt(CDPT_G_I + CDPT_G_EST, :) >= 0);
    end    
    % filtering due to trust region
    cdpt_vld_idx = cdpt_vld_idx & prod(db_cdpt(1:D, :) > tr_bounds(:, 1));
    cdpt_vld_idx = cdpt_vld_idx & prod(db_cdpt(1:D, :) < tr_bounds(:, 2));
    
    % TODO: avoid making copies of filtered out matrices
    if sum(cdpt_vld_idx) > 0
        w_min = ((db_cdpt(CDPT_ROW_F_UB_VTX, :) + db_cdpt(CDPT_ROW_F_UB_HGT, :)) + ...
                 (db_cdpt(CDPT_ROW_F_LB_VTX, :) - db_cdpt(CDPT_ROW_F_LB_HGT, :))) / 2;
        w_vld = -w_min + beta*db_cdpt(CDPT_ROW_F_LAMBDA, :);
        [~, cdpt_idx] = max(w_vld + (~cdpt_vld_idx * -INFTY));

        % the 'while' snippet below is a last-resort way to avoid multiple samples 
        % at the same location, but practically should not run at all!
        cdpt_tmp = db_cdpt(1:D, cdpt_idx);
        % if chosen point x_n is too close to an existing sample, choose another
        while ~isempty(cdpt_tmp) && min(vecnorm(X_n(1:D, :) - repmat(cdpt_tmp, 1, iter))) < 1e-12
            fprintf('Iteration %d: Choosing another exploitation point.\n', iter);
            w_vld(:, cdpt_idx)     = [];
            cdpt_vld_idx(cdpt_idx) = [];
            db_cdpt(:, cdpt_idx)   = [];
            [~, cdpt_idx]         = max(w_vld + (~cdpt_vld_idx * -INFTY));
            cdpt_tmp              = db_cdpt(1:D, cdpt_idx);
        end
    else
        cdpt_tmp = [];
    end
    
    %% Processing and choosing exploitation point from Sobol-generated candidate points
    sbl_vld_idx   = true(1, sbl_size);
    % filtering due to constraints satisfaction
    for g_i = 1:g_len
        SBL_G_I = SBL_ROW_G_INFO + 2*(g_i-1);
        sbl_vld_idx = sbl_vld_idx & (delta * (sbl_cdpt(SBL_G_I + SBL_ROW_GUB, :) + sbl_cdpt(SBL_G_I + SBL_ROW_GLB, :)) / 2 + ...
                                      (1 - delta) * sbl_cdpt(SBL_G_I + SBL_ROW_GLB, :) >= 0);
    end
    sbl_vld = sbl_cdpt(:, sbl_vld_idx);
    if ~isempty(sbl_vld)
        [~, cdpt_idx] = max(-(sbl_vld(SBL_ROW_FUB, :) + sbl_vld(SBL_ROW_FLB, :)) / 2 + beta*(sbl_vld(SBL_ROW_FUB, :) - sbl_vld(SBL_ROW_FLB, :)));

        % the 'while' snippet below is a last-resort way to avoid multiple samples
        % at the same location but practically should not run at all!
        sbl_tmp = sbl_vld(1:D, cdpt_idx);
        % if chosen point x_n is too close to an existing sample, choose another
        while ~isempty(sbl_tmp) && min(vecnorm(X_n(1:D, :) - repmat(sbl_tmp, 1, iter))) < 1e-9 
            fprintf('Iteration %d: Choosing another exploitation point.\n', iter);
            sbl_vld(:, cdpt_idx) = [];
            [~, cdpt_idx]        = min(sbl_vld(end, :));        
            sbl_tmp              = sbl_vld(1:D, cdpt_idx);            
        end
    else
        sbl_tmp = [];
    end
    
    % just selecting if the candidate point from db_cdpt or sbl_cdpt should be chosen
    if ~isempty(cdpt_tmp) && ~isempty(sbl_tmp)
        if max((X_n(XN_ROW_FVAL,:) - feps) - fgam*vecnorm(X_n(1:D,:) - repmat(cdpt_tmp, 1, iter))) < ...
           max((X_n(XN_ROW_FVAL,:) - feps) - fgam*vecnorm(X_n(1:D,:) - repmat(sbl_tmp, 1, iter)))
            x_n_tmp = cdpt_tmp;
        else
            x_n_tmp = sbl_tmp;
        end
    elseif ~isempty(cdpt_tmp)
        x_n_tmp = cdpt_tmp;        
    elseif ~isempty(sbl_tmp)
        x_n_tmp = sbl_tmp;
    else
        x_n_tmp = [];
    end
    
    mode1_thr = opt_z - alpha * fgam;
    if ~isempty(x_n_tmp)
        if opt_z < Inf && max((X_n(XN_ROW_FVAL,:) - feps) - fgam*vecnorm(X_n(1:D,:) - repmat(x_n_tmp, 1, iter))) < mode1_thr
            x_n = x_n_tmp;
            mode_hist(iter) = MODE_EXPLOIT;
            mode_prev       = MODE_EXPLOIT;
            mode1_ok        = true;
        end
    end
    
    %% Algorithm (2) exploration
    if ~mode1_ok
        [~ , cdpt_len] = size(db_cdpt);
        cdpt_vld_idx   = true(1, cdpt_len);
        w_bst          = ones(1, cdpt_len);
        w_unc          = zeros(1, cdpt_len); 

        for g_i = 1:g_len
            CDPT_G_I     = CDPT_ROW_G_INFO + 6*(g_i-1);
%             g_i_vld      = ((db_cdpt(CDPT_G_I + CDPT_G_UB_VTX, :) + db_cdpt(CDPT_G_I + CDPT_G_UB_HGT, :)) ... 
%                           + (db_cdpt(CDPT_G_I + CDPT_G_LB_VTX, :) - db_cdpt(CDPT_G_I + CDPT_G_LB_HGT, :))) / 2 >= 0;
            g_i_vld      = (db_cdpt(CDPT_G_I + CDPT_G_EST, :) >= 0);
            w_bst        = w_bst .* (g_i_vld+1);
            w_unc        = w_unc + (db_cdpt(CDPT_G_I + CDPT_G_LAMBDA, :))/ggam(g_i);
            cdpt_vld_idx = cdpt_vld_idx & g_i_vld;
        end
        w_vld = (db_cdpt(CDPT_ROW_F_LAMBDA, :)) .* cdpt_vld_idx / fgam;
        
        [~, cdpt_idx] = max(((1-delta)*w_vld + delta*w_unc.*w_bst/(2^g_len)).*db_cdpt(end-1, :) + phi.*db_cdpt(end, :));
        x_n = db_cdpt(1:D, cdpt_idx);
        
        % the 'while' snippet below is a last-resort way to avoid multiple samples at the same location
        % but practically should not run at all!
        while ~isempty(x_n) && min(vecnorm(X_n(1:D, :) - repmat(x_n, 1, iter))) < 1e-12  % if chosen point x_n is too close to an existing sample, choose another
            fprintf('Iteration %d: Choosing another exploration point.\n', iter);
            db_cdpt(:, cdpt_idx) = [];
            w_vld(:, cdpt_idx)   = [];
            w_unc(:, cdpt_idx)   = [];
            w_bst(:, cdpt_idx)   = [];
            [~, cdpt_idx] = max(((1-delta)*w_vld + delta*w_unc.*w_bst/(2^g_len)).*db_cdpt(end-1, :) + phi.*db_cdpt(end, :));        
            x_n = db_cdpt(1:D, cdpt_idx);
        end
        mode_hist(iter) = MODE_EXPLORE;
        mode_prev       = MODE_EXPLORE;
    end

    % after i cleared the chosen column in db_cdpt (look at lines 101, 424), why do duplicates still exist?
    if sum(all(X_n(1:D, :) == x_n)) > 0 % this means that i'm selecting a cdpt even though it's already sampled!
        fprintf('Iteration %d: %d duplicate/s detected.\n', iter, sum(all(X_n(1:D, :) == x_n)));
    end

    calc_t_hist(iter) = toc(calc_time);
    if ~exist('opt_x', 'var')
        opt_x = NaN(D, 1);
    end    
    opt_x_hist(:, iter) = opt_x;
    opt_z_hist(iter)    = opt_z;
    opt_c_hist(:, iter) = opt_c;
    
end

x_hist = X_n(1:D, :);    
z_hist = X_n(XN_ROW_FVAL, :);
c_hist = X_n(XN_ROW_GVAL:end, :);

out.opt_x = opt_x;        out.x_hist = x_hist;       out.opt_x_hist = opt_x_hist;
out.opt_z = opt_z;        out.z_hist = z_hist;       out.opt_z_hist = opt_z_hist;
out.opt_c = opt_c;        out.c_hist = c_hist;       out.opt_c_hist = opt_c_hist;

out.calc_t_hist = calc_t_hist;
out.mode_hist   = mode_hist;
out.eps_hist    = eps_hist;
out.gam_hist    = gam_hist;
out.tr_hist     = tr_hist;
