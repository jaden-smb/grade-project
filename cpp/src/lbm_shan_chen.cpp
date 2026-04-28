#include "lbm_shan_chen.h"
#include <cmath>
#include <algorithm>
#include <cstring>

LBMShanChen::LBMShanChen(int nx, int ny, double tau, double G,
                         double rho_liquid, double rho_gas)
    : nx_(nx), ny_(ny), tau_(tau), G_(G),
      rho_liquid_(rho_liquid), rho_gas_(rho_gas),
      eos_type_(EOS_ORIGINAL_SC), rho0_(1.0/1.5),
      T_eos_(0.7 * 0.09433), a_cs_(0.5), b_cs_(2.0) {
    
    n_total_ = nx_ * ny_;
    
    rho_.resize(n_total_, rho_gas_);
    ux_.resize(n_total_, 0.0);
    uy_.resize(n_total_, 0.0);
    f_.resize(n_total_ * Q, 0.0);
    f_eq_.resize(n_total_ * Q, 0.0);
    f_tmp_.resize(n_total_ * Q, 0.0);
    Fx_.resize(n_total_, 0.0);
    Fy_.resize(n_total_, 0.0);
    
    computeEquilibrium();
    std::copy(f_eq_.begin(), f_eq_.end(), f_.begin());
}

LBMShanChen::~LBMShanChen() = default;

void LBMShanChen::initializeDroplet(int center_x, int center_y, double radius) {
    constexpr double interface_width = 10.0;
    double rho_avg = 0.5 * (rho_liquid_ + rho_gas_);
    double rho_diff = 0.5 * (rho_liquid_ - rho_gas_);
    
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int idx = index(i, j);
            double dx = i - center_x;
            double dy = j - center_y;
            double r = std::sqrt(dx*dx + dy*dy);
            rho_[idx] = rho_avg + rho_diff * std::tanh((radius - r) / interface_width);
            ux_[idx] = 0.0;
            uy_[idx] = 0.0;
        }
    }
    
    computeEquilibrium();
    std::copy(f_eq_.begin(), f_eq_.end(), f_.begin());
}

void LBMShanChen::step() {
    computeShanChenForce();
    computeEquilibrium();
    collide();
    stream();
    updateMacroscopic();
}

double LBMShanChen::psi(double rho) const {
    if (eos_type_ == EOS_ORIGINAL_SC) {
        return 1.0 - std::exp(-rho / rho0_);
    }

    // EOS_CARNAHAN_STARLING: Yuan-Schaefer (2006) formulation
    // psi = sqrt(2*(p_CS - rho*cs2) / (G*cs2))
    // p_CS = rho*T*f(eta) - a*rho^2,  eta = b*rho/4,
    //        f(eta) = (1+eta+eta^2-eta^3)/(1-eta)^3
    constexpr double cs2 = 1.0 / 3.0;
    double eta = b_cs_ * rho / 4.0;
    if (eta >= 0.99) eta = 0.99;  // stay below packing singularity
    const double d  = 1.0 - eta;
    const double f_eta = (1.0 + eta + eta*eta - eta*eta*eta) / (d*d*d);
    const double p_cs  = rho * T_eos_ * f_eta - a_cs_ * rho * rho;
    const double arg   = 2.0 * (p_cs - rho * cs2) / (G_ * cs2);
    return (arg > 0.0) ? std::sqrt(arg) : 0.0;
}

void LBMShanChen::computeShanChenForce() {
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int idx = index(i, j);
            double psi_center = psi(rho_[idx]);
            double Fx = 0.0, Fy = 0.0;
            
            for (int q = 0; q < Q; ++q) {
                int i_nb = (i + cx_[q] + nx_) % nx_;
                int j_nb = (j + cy_[q] + ny_) % ny_;
                double psi_nb = psi(rho_[index(i_nb, j_nb)]);
                Fx -= G_ * psi_center * w_[q] * psi_nb * cx_[q];
                Fy -= G_ * psi_center * w_[q] * psi_nb * cy_[q];
            }
            
            Fx_[idx] = Fx;
            Fy_[idx] = Fy;
        }
    }
}

void LBMShanChen::computeEquilibrium() {
    const double u_max_sq = 0.04;

    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int idx = index(i, j);
            double rho = rho_[idx];
            double ux, uy;
            if (rho > 1e-10) {
                ux = ux_[idx] + 0.5 * Fx_[idx] / rho;
                uy = uy_[idx] + 0.5 * Fy_[idx] / rho;
            } else {
                ux = uy = 0.0;
            }

            double u_sq = ux*ux + uy*uy;
            if (u_sq > u_max_sq) {
                double scale = std::sqrt(u_max_sq / u_sq);
                ux *= scale;
                uy *= scale;
                u_sq = u_max_sq;
                ++clamp_count_;
            }

            for (int q = 0; q < Q; ++q) {
                double cu = cx_[q]*ux + cy_[q]*uy;
                f_eq_[f_index(i, j, q)] = w_[q] * rho * (1.0 + 3.0*cu + 4.5*cu*cu - 1.5*u_sq);
            }

            double feq_sum = 0.0;
            for (int q = 0; q < Q; ++q)
                feq_sum += f_eq_[f_index(i, j, q)];
            if (std::abs(feq_sum - rho) > 1e-14) {
                double correction = rho / feq_sum;
                for (int q = 0; q < Q; ++q)
                    f_eq_[f_index(i, j, q)] *= correction;
            }
        }
    }
}

void LBMShanChen::collide() {
    double omega = 1.0 / tau_;
    for (int idx = 0; idx < n_total_; ++idx) {
        for (int q = 0; q < Q; ++q) {
            int f_idx = idx * Q + q;
            f_[f_idx] -= omega * (f_[f_idx] - f_eq_[f_idx]);
        }
    }
}

void LBMShanChen::stream() {
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            for (int q = 0; q < Q; ++q) {
                int i_src = (i - cx_[q] + nx_) % nx_;
                int j_src = (j - cy_[q] + ny_) % ny_;
                f_tmp_[f_index(i, j, q)] = f_[f_index(i_src, j_src, q)];
            }
        }
    }
    f_.swap(f_tmp_);
}

void LBMShanChen::updateMacroscopic() {
    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int idx = index(i, j);
            double rho = 0.0, jx = 0.0, jy = 0.0;
            
            for (int q = 0; q < Q; ++q) {
                double fv = f_[f_index(i, j, q)];
                rho += fv;
                jx  += fv * cx_[q];
                jy  += fv * cy_[q];
            }
            
            rho_[idx] = rho;
            if (rho > 1e-10) {
                ux_[idx] = (jx + 0.5 * Fx_[idx]) / rho;
                uy_[idx] = (jy + 0.5 * Fy_[idx]) / rho;
            } else {
                ux_[idx] = uy_[idx] = 0.0;
            }
        }
    }
}
