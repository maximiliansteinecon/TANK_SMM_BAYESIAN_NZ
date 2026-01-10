// ====================================================================
// FILE 3: BAYESIAN ESTIMATION
// ====================================================================

var c_SP c_SA c y k i g R_rate rK w L L_SP L_SA T T_SP T_SA A eps_G gy gc;
varexo eta_A eta_G;
parameters beta delta alpha gamma_y omega rho_A l_SP l_SA rho_G;

// -------------------------------------------------------------------------
// 1. OBSERVABLES
// -------------------------------------------------------------------------
verbatim;
options_.varobs = {'gy', 'gc'};
options_.varobs_id = [strmatch('gy', M_.endo_names, 'exact'); 
                      strmatch('gc', M_.endo_names, 'exact')];
end;

// -------------------------------------------------------------------------
// 2. CALIBRATION (SMM RESULTS)
// -------------------------------------------------------------------------
beta     = 0.9950;  
gamma_y  = 1.0080;   
rho_A    = 0.9988;   
omega    = 0.9361;    

// Standard Parameters
alpha    = 1/3;       
l_SP     = 1/3;       
l_SA     = 1/3;       
delta    = 0.025;     
rho_G    = 0.5;       

// -------------------------------------------------------------------------
// 3. MODEL
// -------------------------------------------------------------------------
model;
    c_SP = w*L_SP - T_SP;
    c_SA = (c - omega*c_SP)/(1-omega);
    k = ((1-delta)/gamma_y)*k(-1) + (1/gamma_y)*i(-1);
    1/c_SA = beta*(R_rate(+1)/gamma_y)*(1/c_SA(+1));
    R_rate = rK(+1) + (1-delta);
    
    y = A * k(-1)^alpha * L^(1-alpha);
    w = (1-alpha) * y / L;
    rK = alpha * y / k(-1);
    
    g = T;
    g = 0.2 * y * eps_G;
    
    L = omega*L_SP + (1-omega)*L_SA;
    c = y - i - g;
    
    log(A) = rho_A*log(A(-1)) + eta_A;
    log(eps_G) = rho_G*log(eps_G(-1)) + eta_G;
    
    L_SP = l_SP; L_SA = l_SA;
    T_SP = T; T_SA = T;
    
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
// 4. BAYESIAN ESTIMATION
// -------------------------------------------------------------------------
shocks;
    var eta_A = (0.0098)^2; 
    var eta_G = 0.01; 
end;

estimated_params;
    rho_G, beta_pdf, 0.5, 0.2;
    stderr eta_G, inv_gamma_pdf, 0.01, inf;
end;

// Running with Smoother
estimation(datafile='macro_data_NZ.csv', mh_replic=2000, mh_nblocks=2, 
           mh_jscale=0.3, order=1, 
           smoother, forecast=0, tex);

// -------------------------------------------------------------------------
// 5. ANALYSIS: SHOCK DECOMPOSITION
// -------------------------------------------------------------------------
// Plots the contribution of eta_A vs eta_G.
shock_decomposition y c;

// -------------------------------------------------------------------------
// 6. COUNTERFACTUALS
// -------------------------------------------------------------------------
// Generates the IRFs for different omega values.
@#for w_val in ["0.1", "0.3", "0.5"]
    set_param_value('omega', @{w_val});
    stoch_simul(order=1, irf=20) y c;
@#endfor
