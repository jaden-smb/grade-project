#ifndef LBM_SHAN_CHEN_H
#define LBM_SHAN_CHEN_H

#include <vector>
#include <memory>

static constexpr int EOS_ORIGINAL_SC        = 0;
static constexpr int EOS_CARNAHAN_STARLING  = 1;

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

    void setEosType(int eos_type) { eos_type_ = eos_type; }
    int  getEosType() const       { return eos_type_; }

    void setTemperature(double T) { T_eos_ = T; }
    double getTemperature() const { return T_eos_; }

    void setRho0(double rho0) { rho0_ = rho0; }
    double getRho0() const { return rho0_; }

private:
    int nx_, ny_;
    int n_total_;

    double tau_;
    double G_;
    double rho_liquid_;
    double rho_gas_;

    int    eos_type_ = EOS_ORIGINAL_SC;
    double rho0_    = 1.0 / 1.5;
    double T_eos_   = 0.7 * 0.09433;
    double a_cs_    = 0.5;
    double b_cs_    = 2.0;
    
    static const int Q = 9;
    static constexpr double w_[Q] = {
        4.0/9.0,
        1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,
        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0
    };
    static constexpr int cx_[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    static constexpr int cy_[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    
    std::vector<double> rho_;
    std::vector<double> ux_, uy_;
    std::vector<double> f_, f_eq_, f_tmp_;
    std::vector<double> Fx_, Fy_;
    
    int index(int i, int j) const { return j * nx_ + i; }
    int f_index(int i, int j, int q) const { return (j * nx_ + i) * Q + q; }

    void computeEquilibrium();
    void computeShanChenForce();
    void collide();
    void stream();
    void updateMacroscopic();
    double psi(double rho) const;

    int clamp_count_ = 0;
};

#endif
