#ifndef LBM_SHAN_CHEN_H
#define LBM_SHAN_CHEN_H

#include <vector>
#include <memory>

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

private:
    int nx_, ny_;
    int n_total_;

    double tau_;
    double G_;
    double rho_liquid_;
    double rho_gas_;
    
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
    double psi(double rho) const; // 1 - exp(-1.5*rho)

    int clamp_count_ = 0;
};

#endif
