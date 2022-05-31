clear all; close all; clc;
addpath('test_fns', 'util');
if ~exist('out_dat', 'dir')
    mkdir('out_dat');
end
addpath('out_dat');

% Declare definitions, including default settings for SMGO-D
defines;

% Possible options in test problems: FN_T1 FN_STYB_10D G04 G05 G06 G08 G09 G10 G12 G23 G24 
test_fn_list  = [ FN_T1 ];

%% Automated run for all test functions, opt algos, and exploit weighting
for test_fn = test_fn_list
    test_init;
    x0_list_normd = rand( opt_dim, max_trial );
    fprintf( "Solving %s using SMGO-D.\n", test_fn_names{test_fn} );

    opt_x_stack       = [];       calc_t_hist_stack = [];
    opt_z_stack       = [];       mode_hist_stack   = [];
    opt_c_stack       = [];       eps_hist_stack    = [];
    x_hist_stack      = [];       gam_hist_stack    = [];
    z_hist_stack      = [];       tr_hist_stack     = [];
    c_hist_stack      = [];
    opt_x_hist_stack  = [];
    opt_z_hist_stack  = [];
    opt_c_hist_stack  = [];

    for trial = 1:max_trial
        options.startpt = x0_list_normd( : , trial );
        tic;
        smgo_d_res = opt_smgo_d( options );
        toc;
        opt_x       = smgo_d_res.opt_x;
        opt_z       = smgo_d_res.opt_z;
        opt_c       = smgo_d_res.opt_c;
        x_hist      = smgo_d_res.x_hist;
        z_hist      = smgo_d_res.z_hist;
        c_hist      = smgo_d_res.c_hist;                    
        opt_x_hist  = smgo_d_res.opt_x_hist;
        opt_z_hist  = smgo_d_res.opt_z_hist;
        opt_c_hist  = smgo_d_res.opt_c_hist;
        calc_t_hist = smgo_d_res.calc_t_hist;
        mode_hist   = smgo_d_res.mode_hist;
        eps_hist    = smgo_d_res.eps_hist;
        gam_hist    = smgo_d_res.gam_hist;
        tr_hist     = smgo_d_res.tr_hist;

        opt_x_stack       = [ opt_x_stack; opt_x' ];
        opt_z_stack       = [ opt_z_stack; opt_z ];  
        opt_c_stack       = [ opt_c_stack; opt_c' ];                                      
        x_hist_stack      = [ x_hist_stack; x_hist ];
        z_hist_stack      = [ z_hist_stack; z_hist ];
        c_hist_stack      = [ c_hist_stack; c_hist ];
        opt_x_hist_stack  = [ opt_x_hist_stack; opt_x_hist ];
        opt_z_hist_stack  = [ opt_z_hist_stack; opt_z_hist ];
        opt_c_hist_stack  = [ opt_c_hist_stack; opt_c_hist ];
        calc_t_hist_stack = [ calc_t_hist_stack; calc_t_hist ];
        mode_hist_stack   = [ mode_hist_stack; mode_hist ];
        eps_hist_stack    = [ eps_hist_stack; eps_hist ];
        gam_hist_stack    = [ gam_hist_stack; gam_hist ];
        tr_hist_stack     = [ tr_hist_stack; tr_hist ];
        fprintf("Trial %4d: %f\n", trial, opt_z)
    end

    % save important data about the optimization runs
    mat_name = sprintf( "out_dat/SMGO_D-%s-%.2f.mat", test_fn_names{test_fn}, beta );
    save( mat_name, 'x0_list_normd', 'opt_x_stack', 'opt_z_stack', ...
                    'x_hist_stack', 'z_hist_stack', 'c_hist_stack', ...
                    'opt_x_hist_stack', 'opt_z_hist_stack', 'calc_t_hist_stack', ...
                    'mode_hist_stack', 'eps_hist_stack', 'gam_hist_stack', 'tr_hist_stack' );
end
