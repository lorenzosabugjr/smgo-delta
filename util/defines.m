test_case = 0;
FN_G04    = 1;
FN_G05    = 2;
FN_G06    = 3;
FN_G08    = 4;
FN_G09    = 5;
FN_G10    = 6;
FN_G12    = 7;
FN_G23    = 8;
FN_G24    = 9;
FN_T1     = 10;
FN_STYB_10D = 11;
test_fn_names = { 'G04', 'G05', 'G06', 'G08', 'G09', 'G10', 'G12', 'G23', 'G24', 'T1', 'STYB_10D' };

max_trial = 10;

% SMGO-D settings
alpha = 0.005; sobol = 1;
beta  = 0.1;   trust = 1;
delta = 0.2; 
options.maxiter     = 500;
options.alpha       = alpha;
options.beta        = beta;
options.delta       = delta;
options.sobol       = sobol;
options.trustregion = trust;