#include "lbm_shan_chen.h"
#include <cmath>
#include <algorithm>
#include <cstring>

LBMShanChen::LBMShanChen(int nx, int ny, double tau, double G, 
                         double rho_liquid, double rho_gas)
    : nx_(nx), ny_(ny), tau_(tau), G_(G), 
      rho_liquid_(rho_liquid), rho_gas_(rho_gas) {
    
    n_total_ = nx_ * ny_;
    
    // Allocate memory for all fields
    rho_.resize(n_total_, rho_gas_);
    ux_.resize(n_total_, 0.0);
    uy_.resize(n_total_, 0.0);
    f_.resize(n_total_ * Q, 0.0);
    f_eq_.resize(n_total_ * Q, 0.0);
    Fx_.resize(n_total_, 0.0);
    Fy_.resize(n_total_, 0.0);
    
    // Initialize distribution functions to equilibrium
    computeEquilibrium();
    std::copy(f_eq_.begin(), f_eq_.end(), f_.begin());
}

LBMShanChen::~LBMShanChen() = default;

void LBMShanChen::initializeDroplet(int center_x, int center_y, double radius) {
    double radius_sq = radius * radius;
    
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int idx = index(i, j);
            
            // Distance from center
            double dx = i - center_x;
            double dy = j - center_y;
            double dist_sq = dx * dx + dy * dy;
            
            // Set density: liquid inside droplet, gas outside
            if (dist_sq <= radius_sq) {
                rho_[idx] = rho_liquid_;
            } else {
                rho_[idx] = rho_gas_;
            }
            
            // Initialize velocities to zero
            ux_[idx] = 0.0;
            uy_[idx] = 0.0;
        }
    }
    
    // Recompute equilibrium and distribution functions
    computeEquilibrium();
    std::copy(f_eq_.begin(), f_eq_.end(), f_.begin());
}

void LBMShanChen::step() {
    // Main LBM time step sequence:
    // 1. Compute Shan-Chen interaction force (based on current density)
    computeShanChenForce();
    
    // 2. Compute equilibrium distribution
    // Note: We'll incorporate force into equilibrium by modifying velocity
    computeEquilibrium();
    
    // 3. Apply collision operator (BGK)
    collide();
    
    // 4. Stream distribution functions
    stream();
    
    // 5. Update macroscopic variables (density, velocity with force)
    updateMacroscopic();
}

double LBMShanChen::psi(double rho) const {
    // Shan-Chen pseudo-potential: psi(rho) = 1 - exp(-rho)
    return 1.0 - std::exp(-rho);
}

void LBMShanChen::computeShanChenForce() {
    // Shan-Chen force: F = -G * psi(rho) * sum(w_i * psi(rho_neighbor) * c_i)
    // This force creates phase separation when G < 0
    
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int idx = index(i, j);
            double rho_center = rho_[idx];
            double psi_center = psi(rho_center);
            
            double Fx = 0.0;
            double Fy = 0.0;
            
            // Sum over all velocity directions
            for (int q = 0; q < Q; ++q) {
                // Neighbor position
                int i_nb = i + cx_[q];
                int j_nb = j + cy_[q];
                
                // Periodic boundary conditions
                if (i_nb < 0) i_nb = nx_ - 1;
                if (i_nb >= nx_) i_nb = 0;
                if (j_nb < 0) j_nb = ny_ - 1;
                if (j_nb >= ny_) j_nb = 0;
                
                int idx_nb = index(i_nb, j_nb);
                double rho_nb = rho_[idx_nb];
                double psi_nb = psi(rho_nb);
                
                // Accumulate force contribution
                Fx -= G_ * psi_center * w_[q] * psi_nb * cx_[q];
                Fy -= G_ * psi_center * w_[q] * psi_nb * cy_[q];
            }
            
            Fx_[idx] = Fx;
            Fy_[idx] = Fy;
        }
    }
}

void LBMShanChen::computeEquilibrium() {
    // Equilibrium distribution function:
    // f_i^eq = w_i * rho * [1 + 3*(c_i·u) + 4.5*(c_i·u)^2 - 1.5*u^2]
    // For force implementation, we use velocity with force contribution:
    // u_eq = u + F/(2*rho)
    
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int idx = index(i, j);
            double rho = rho_[idx];
            
            // Velocity with force contribution for equilibrium
            double ux, uy;
            if (rho > 1e-10) {
                ux = ux_[idx] + 0.5 * Fx_[idx] / rho;
                uy = uy_[idx] + 0.5 * Fy_[idx] / rho;
            } else {
                ux = 0.0;
                uy = 0.0;
            }
            
            double u_sq = ux * ux + uy * uy;
            
            for (int q = 0; q < Q; ++q) {
                double c_dot_u = cx_[q] * ux + cy_[q] * uy;
                double f_eq = w_[q] * rho * (1.0 + 3.0 * c_dot_u + 
                                              4.5 * c_dot_u * c_dot_u - 
                                              1.5 * u_sq);
                
                int f_idx = f_index(i, j, q);
                f_eq_[f_idx] = f_eq;
            }
        }
    }
}

void LBMShanChen::collide() {
    // BGK collision operator:
    // f_i' = f_i - (f_i - f_i^eq) / tau
    
    double omega = 1.0 / tau_;  // Collision frequency
    
    for (int idx = 0; idx < n_total_; ++idx) {
        for (int q = 0; q < Q; ++q) {
            int f_idx = idx * Q + q;
            f_[f_idx] = f_[f_idx] - omega * (f_[f_idx] - f_eq_[f_idx]);
        }
    }
}

void LBMShanChen::stream() {
    // Streaming step: f_i(x + c_i, t+1) = f_i'(x, t)
    // We need a temporary buffer to avoid overwriting during streaming
    
    std::vector<double> f_new(n_total_ * Q, 0.0);
    
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            for (int q = 0; q < Q; ++q) {
                // Source position
                int i_src = i - cx_[q];
                int j_src = j - cy_[q];
                
                // Periodic boundary conditions
                if (i_src < 0) i_src = nx_ - 1;
                if (i_src >= nx_) i_src = 0;
                if (j_src < 0) j_src = ny_ - 1;
                if (j_src >= ny_) j_src = 0;
                
                int f_src = f_index(i_src, j_src, q);
                int f_dst = f_index(i, j, q);
                
                f_new[f_dst] = f_[f_src];
            }
        }
    }
    
    // Copy back
    std::copy(f_new.begin(), f_new.end(), f_.begin());
}

void LBMShanChen::updateMacroscopic() {
    // Update density and velocity from distribution functions
    // rho = sum(f_i)
    // u = (sum(f_i * c_i) + F/2) / rho
    // This accounts for the force effect on momentum
    
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int idx = index(i, j);
            
            // Compute density and momentum
            double rho = 0.0;
            double jx = 0.0;
            double jy = 0.0;
            
            for (int q = 0; q < Q; ++q) {
                int f_idx = f_index(i, j, q);
                double f_val = f_[f_idx];
                rho += f_val;
                jx += f_val * cx_[q];
                jy += f_val * cy_[q];
            }
            
            // Update density
            rho_[idx] = rho;
            
            // Update velocity with force contribution
            // u = (j + F/2) / rho
            if (rho > 1e-10) {  // Avoid division by zero
                ux_[idx] = (jx + 0.5 * Fx_[idx]) / rho;
                uy_[idx] = (jy + 0.5 * Fy_[idx]) / rho;
            } else {
                ux_[idx] = 0.0;
                uy_[idx] = 0.0;
            }
        }
    }
}

