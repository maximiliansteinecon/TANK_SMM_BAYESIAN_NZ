// ====================================================================
// FILE 2: SMM ESTIMATION
// ====================================================================
var c_SP c_SA c y k i g R_rate rK w L L_SP L_SA T T_SP T_SA A eps_G gy gc;
varexo eta_A eta_G;
parameters beta delta alpha gamma_y omega rho_A l_SP l_SA rho_G;

// -------------------------------------------------------------------------
// 1. OBSERVABLES
// -------------------------------------------------------------------------
verbatim;
options_.varobs = {'gy', 'gc', 'R_rate'};
end;

// -------------------------------------------------------------------------
// 2. CALIBRATION
// -------------------------------------------------------------------------
omega    = 0.3;       
alpha    = 1/3;       
rho_A    = 0.9;       
l_SP     = 1/3;       
l_SA     = 1/3;       
beta     = 0.99;      
delta    = 0.025;     
gamma_y  = 1.005;     
rho_G    = 0.0;

// -------------------------------------------------------------------------
// 3. MODEL
// -------------------------------------------------------------------------
model;
    // Househoulds
    c_SP = w*L_SP - T_SP;
    c_SA = (c - omega*c_SP)/(1-omega);
    k = ((1-delta)/gamma_y)*k(-1) + (1/gamma_y)*i(-1);
    1/c_SA = beta*(R_rate(+1)/gamma_y)*(1/c_SA(+1));
    R_rate = rK(+1) + (1-delta);
    
    // Firms
    y = A * k(-1)^alpha * L^(1-alpha);
    w = (1-alpha) * y / L;
    rK = alpha * y / k(-1);
    
    // Government
    g = T;
    g = 0.2 * y * eps_G;
    
    // Aggregation
    L = omega*L_SP + (1-omega)*L_SA;
    c = y - i - g;
    
    // Exogenous Processes
    log(A) = rho_A*log(A(-1)) + eta_A;
    log(eps_G) = rho_G*log(eps_G(-1)) + eta_G;
    
    // Identities
    L_SP = l_SP; L_SA = l_SA;
    T_SP = T; T_SA = T;
    
    // Observables
    gy = log(gamma_y) + log(y) - log(y(-1));
    gc = log(gamma_y) + log(c) - log(c(-1));
end;

steady_state_model;
    A = 1; eps_G = 1; 
    L_SP = l_SP; L_SA = l_SA; L = omega*L_SP + (1-omega)*L_SA;
    R_rate = gamma_y/beta; rK = R_rate - (1-delta);
    y = ((alpha^alpha * L^(1-alpha))/(rK^alpha))^(1/(1-alpha));
    k = (alpha/rK)*y; i = (gamma_y-(1-delta))*k; w = (1-alpha)*y/L;
    g = 0.2*y*eps_G; T = g; T_SP = T; T_SA = T;
    c_SP = w*L_SP - T_SP; c = y - i - g; c_SA = (c - omega*c_SP)/(1-omega);
    gy = log(gamma_y); gc = log(gamma_y);
end;

steady;

// -------------------------------------------------------------------------
// 4. SMM ESTIMATION
// -------------------------------------------------------------------------
shocks;
    var eta_G = 0; 
end;

matched_moments;
    gy; gc; R_rate;
    gy*gy; gc*gc;
end;

estimated_params;
    // Format: param_name, INITIAL_GUESS, lower_bound, upper_bound;
    beta, 0.995, 0.9, 0.9999;
    gamma_y, 1.008, 1.0, 1.02;
    stderr eta_A, 0.01, 0.0001, 0.1;
    
    rho_A, 0.99, 0.01, 0.9999;
    
    omega, 0.5, 0.0, 0.99;
end;

// Using mode_compute=4 (CSMINWEL)
method_of_moments(mom_method = SMM, datafile = 'macro_data_NZ.csv', mode_compute = 5, order = 1);
