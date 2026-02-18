#ifndef LBM_SHAN_CHEN_H
#define LBM_SHAN_CHEN_H

#include <vector>
#include <memory>

// D2Q9 LBM with Shan-Chen pseudo-potential for liquid-gas phase separation.
// Uses BGK collision and psi(rho) = 1 - exp(-rho).
class LBMShanChen {
public:
    // tau: relaxation time (>0.5), G: cohesion (~-4.5 to -7, critical ~-4)
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
    std::vector<double> f_, f_eq_;
    std::vector<double> Fx_, Fy_;
    
    int index(int i, int j) const { return j * nx_ + i; }
    int f_index(int i, int j, int q) const { return (j * nx_ + i) * Q + q; }
    
    void computeEquilibrium();   // f_eq = w*rho*(1 + 3cu + 4.5cu^2 - 1.5u^2)
    void computeShanChenForce(); // F = -G*psi(rho)*sum(w*psi(nb)*c)
    void collide();              // BGK: f' = f - (f - f_eq)/tau
    void stream();               // f(x+c, t+1) = f'(x, t)
    void updateMacroscopic();    // rho = sum(f), u = (sum(f*c) + F/2)/rho
    double psi(double rho) const; // 1 - exp(-rho)
};

#endif // LBM_SHAN_CHEN_H
