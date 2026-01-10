// ====================================================================
// FILE 1 - Two-Agent Household DSGE Model - Baseline calibration .inc 
// Based on Setup and Report by Maximilian Stein
// Model: Hand-to-Mouth (SP) and Wealthy (SA) Households
// Version: 17th December 2025
// ====================================================================

var
    c_SP c_SA c y k i g r rK w L L_SP L_SA T T_SP T_SA A;
    // gy gc
varexo
    eta_A;
parameters
    beta delta alpha gamma_y omega rho_A sigma_A l_SP l_SA;
    // rho_G sigma_G

// ---- Calibration from (set-up) ----
omega    = 0.3;       // Share of hand-to-mouth households
alpha    = 1/3;       // Capital share in production
rho_A    = 0.9;       // TFP persistence
sigma_A  = 0.01;      // TFP shock standard deviation
l_SP     = 1/3;       // inelastic labor, spenders
l_SA     = 1/3;       // inelastic labor, savers

// Calibration - Empirical definitions
beta     = 0.99;      // Quarterly discount factor
delta    = 0.025;     // Quarterly depreciation rate
gamma_y  = 1.005;     // Quarterly growth rate (approx 2% annual)

// -------------------------------------------------------------------------
// Model
// -------------------------------------------------------------------------
model;
    // 1 SP budget
    c_SP = w*L_SP - T_SP;

    // 2 SA budget
    // Modified due to conflict with goods market clearing otherwise
    c_SA = (c - omega*c_SP)/(1-omega);

    // 3 Capital accumulation (detrended)
    k = ((1-delta)/gamma_y)*k(-1) + (1/gamma_y)*i(-1);

    // 4 Euler (bonds)
    1/c_SA = beta*(r(+1)/gamma_y)*(1/c_SA(+1));

    // 5 No-arbitrage
    r = rK(+1) + (1-delta);

    // 6 Production (detrended)
    y = A * k(-1)^alpha * L^(1-alpha);

    // 7 Labor demand
    w = (1-alpha) * y / L;

    // 8 Capital demand
    rK = alpha * y / k(-1);

    // 9 Government balanced budget, G = taxshare*Y, uniform taxes
    g = T;

    // 10 Spending Rule
    g = 0.2 * y;

    // 11 Aggregate Labor
    L = omega*L_SP + (1-omega)*L_SA;

    // 12 Aggregate consumption
    // c = omega*c_SP + (1-omega)*c_SA;

    // 13 Bond Market Clearing
    // b = 0, but as imposed in (2) commented out here

    // 14 TFP process
    log(A) = rho_A*log(A(-1)) + eta_A;

    // 15 Goods market clearing
    c = y - i - g;

    // 16 Inelastic labor supplies
    L_SP = l_SP;
    L_SA = l_SA;

    // 17 Uniform taxes across agents
    T_SP = T;
    T_SA = T;
end;

// -------------------------------------------------------------------------
// Steady state (analytical block)
// Dynare supports steady_state_model to provide steady state equations.
// -------------------------------------------------------------------------
steady_state_model;
    // Labour Components
    A = 1;
    L_SP = l_SP;
    L_SA = l_SA;
    L = omega*L_SP + (1-omega)*L_SA;

    // From Bond Euler in detrended steady state: 1 = beta*(r/gy)
    r = gamma_y/beta;

    // From Capital Euler in detrended steady state: r = rk + 1 - delta
    rK = r - (1-delta);

    // Solve for y,k from:
    // rk = alpha*y/k  and  y = k^alpha * L^(1-alpha)
    y = ( (alpha^alpha * L^(1-alpha)) / (rK^alpha) )^(1/(1-alpha));
    k = (alpha/rK)*y;

    // Investment from detrended capital accumulation steady state:
    // k = ((1-delta)/gy)k + (1/gy)i  =>  i = (gy - (1-delta))*k
    i = (gamma_y-(1-delta))*k;

    // Wages
    w = (1-alpha)*y/L;

    // Government Components
    g = 0.2*y;
    T = g;
    T_SP = T;
    T_SA = T;

    // Consumption Components
    c_SP = w*L_SP - T_SP;
    c = y - i - g;
    c_SA = (c - omega*c_SP)/(1-omega);

  
end;

// -------------------------------------------------------------------------
// Shocks and simulation
// -------------------------------------------------------------------------
steady;
check;
resid(non_zero);
model_diagnostics;

shocks;
    var eta_A = sigma_A^2;
end;

stoch_simul(order=1, irf=40, periods=1000);