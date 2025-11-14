#ifndef LBM_SHAN_CHEN_H
#define LBM_SHAN_CHEN_H

#include <vector>
#include <memory>

/**
 * Two-Phase Lattice Boltzmann Method using Shan-Chen Model
 * 
 * This class implements a D2Q9 LBM simulation with Shan-Chen pseudo-potential
 * force for liquid-gas phase separation. The model uses:
 * - D2Q9 velocity set (9 discrete velocities in 2D)
 * - BGK collision operator
 * - Shan-Chen interaction force with psi(rho) = 1 - exp(-rho)
 * 
 * Mathematical Model:
 * - Distribution function: f_i(x, t) represents particles moving in direction i
 * - Collision: f_i' = f_i - (f_i - f_i^eq) / tau
 * - Streaming: f_i(x + c_i, t+1) = f_i'(x, t)
 * - Shan-Chen force: F = -G * psi(rho) * sum(w_i * psi(rho_neighbor) * c_i)
 * - Density: rho = sum(f_i)
 * - Velocity: u = (sum(f_i * c_i) + F/2) / rho
 */
class LBMShanChen {
public:
    /**
     * Constructor
     * @param nx Grid width
     * @param ny Grid height
     * @param tau Relaxation time (controls viscosity)
     * @param G Cohesion parameter (negative for liquid-gas separation)
     * @param rho_liquid Initial liquid density
     * @param rho_gas Initial gas density
     */
    LBMShanChen(int nx, int ny, double tau, double G, 
                double rho_liquid, double rho_gas);
    
    ~LBMShanChen();
    
    /**
     * Initialize a circular liquid droplet in gas
     * @param center_x X coordinate of droplet center
     * @param center_y Y coordinate of droplet center
     * @param radius Droplet radius in lattice units
     */
    void initializeDroplet(int center_x, int center_y, double radius);
    
    /**
     * Perform one LBM time step
     * This includes: collision, force application, streaming, and density update
     */
    void step();
    
    // Accessors for Python
    const std::vector<double>& getDensity() const { return rho_; }
    const std::vector<double>& getVelocityX() const { return ux_; }
    const std::vector<double>& getVelocityY() const { return uy_; }
    int getNx() const { return nx_; }
    int getNy() const { return ny_; }
    
private:
    // Grid dimensions
    int nx_, ny_;
    int n_total_;  // Total number of lattice sites
    
    // LBM parameters
    double tau_;           // Relaxation time
    double G_;             // Cohesion parameter
    double rho_liquid_;    // Liquid density
    double rho_gas_;       // Gas density
    
    // D2Q9 velocity set (9 directions)
    static const int Q = 9;
    static constexpr double w_[Q] = {
        4.0/9.0,   // i=0: rest particle
        1.0/9.0,   // i=1-4: cardinal directions
        1.0/9.0,
        1.0/9.0,
        1.0/9.0,
        1.0/36.0,  // i=5-8: diagonal directions
        1.0/36.0,
        1.0/36.0,
        1.0/36.0
    };
    static constexpr int cx_[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    static constexpr int cy_[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    
    // Lattice fields
    std::vector<double> rho_;      // Density field
    std::vector<double> ux_;       // X-velocity field
    std::vector<double> uy_;       // Y-velocity field
    std::vector<double> f_;        // Distribution functions (n_total * Q)
    std::vector<double> f_eq_;     // Equilibrium distribution functions
    std::vector<double> Fx_;       // X-component of Shan-Chen force
    std::vector<double> Fy_;       // Y-component of Shan-Chen force
    
    // Helper functions
    int index(int i, int j) const { return j * nx_ + i; }
    int f_index(int i, int j, int q) const { return (j * nx_ + i) * Q + q; }
    
    /**
     * Compute equilibrium distribution function
     * f_i^eq = w_i * rho * [1 + 3*(c_i·u) + 4.5*(c_i·u)^2 - 1.5*u^2]
     */
    void computeEquilibrium();
    
    /**
     * Compute Shan-Chen pseudo-potential force
     * F = -G * psi(rho) * sum(w_i * psi(rho_neighbor) * c_i)
     * where psi(rho) = 1 - exp(-rho)
     */
    void computeShanChenForce();
    
    /**
     * Apply BGK collision operator
     * f_i' = f_i - (f_i - f_i^eq) / tau
     */
    void collide();
    
    /**
     * Stream distribution functions
     * f_i(x + c_i, t+1) = f_i'(x, t)
     */
    void stream();
    
    /**
     * Update density and velocity from distribution functions
     * rho = sum(f_i)
     * u = (sum(f_i * c_i) + F/2) / rho
     */
    void updateMacroscopic();
    
    /**
     * Shan-Chen pseudo-potential function
     * psi(rho) = 1 - exp(-rho)
     */
    double psi(double rho) const;
};

#endif // LBM_SHAN_CHEN_H

