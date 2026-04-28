#ifndef LBM_SHAN_CHEN_H
#define LBM_SHAN_CHEN_H

#include <vector>
#include <memory>

// EoS selector constants
static constexpr int EOS_ORIGINAL_SC        = 0;  // psi = 1 - exp(-rho/rho0)
static constexpr int EOS_CARNAHAN_STARLING  = 1;  // Yuan-Schaefer with CS EoS

class LBMShanChen {
public:
    LBMShanChen(int nx, int ny, double tau, double G,
                double rho_liquid, double rho_gas);

    ~LBMShanChen();

    void initializeDroplet(int center_x, int center_y, double radius);
    void step();

    const std::vector<double>& getDensity() const { return rho_; }
    const std::vector<double>& getVelocityX() const { return ux_; }
    const std::vector<double>& getVelocityY() const { return uy_; }
    int getNx() const { return nx_; }
    int getNy() const { return ny_; }
    void setG(double G) { G_ = G ; }
    double getG() const { return G_; }
    double getTau() const { return tau_; }
    int getClampCount() const { return clamp_count_; }
    void resetClampCount() { clamp_count_ = 0; }

    // EoS selection (call before running steps; default = EOS_ORIGINAL_SC)
    void setEosType(int eos_type) { eos_type_ = eos_type; }
    int  getEosType() const       { return eos_type_; }

    // Temperature for non-ideal EoS (Carnahan-Starling); default = 0.7 * Tc
    void setTemperature(double T) { T_eos_ = T; }
    double getTemperature() const { return T_eos_; }

    // rho0 for original SC: psi = 1 - exp(-rho/rho0); default = 1/1.5
    void setRho0(double rho0) { rho0_ = rho0; }
    double getRho0() const { return rho0_; }

private:
    int nx_, ny_;
    int n_total_;

    double tau_;
    double G_;
    double rho_liquid_;
    double rho_gas_;

    // EoS parameters (all have safe defaults matching original behaviour)
    int    eos_type_ = EOS_ORIGINAL_SC;
    double rho0_    = 1.0 / 1.5;   // original SC: psi = 1 - exp(-1.5*rho)
    double T_eos_   = 0.7 * 0.09433; // CS temperature (0.7 * Tc, Tc~0.09433 for a=0.5,b=2)
    double a_cs_    = 0.5;          // CS mean-field attraction
    double b_cs_    = 2.0;          // CS hard-sphere repulsion
    
    static const int Q = 9;
    static constexpr double w_[Q] = {
        4.0/9.0,               // rest
        1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,    // cardinal
        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0  // diagonal
    };
    static constexpr int cx_[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    static constexpr int cy_[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    
    std::vector<double> rho_;
    std::vector<double> ux_, uy_;
    std::vector<double> f_, f_eq_, f_tmp_;
    std::vector<double> Fx_, Fy_;
    
    int index(int i, int j) const { return j * nx_ + i; }
    int f_index(int i, int j, int q) const { return (j * nx_ + i) * Q + q; }

    void computeEquilibrium();   // f_eq(u_eq), u_eq = u_phys + F/(2*rho)
    void computeShanChenForce(); // F = -G*psi(rho)*sum(w*psi(nb)*c)
    void collide();              // BGK: f' = f - (f - f_eq)/tau
    void stream();               // f(x+c, t+1) = f'(x, t)
    void updateMacroscopic();    // rho = sum(f), u = (sum(f*c) + F/2)/rho
    double psi(double rho) const; // dispatches based on eos_type_

    int clamp_count_ = 0;
};

#endif
