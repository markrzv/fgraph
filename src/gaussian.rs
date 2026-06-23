use std::ops::{Div, Mul};

#[derive(Debug, Clone, Copy)]
pub struct Gaussian {
    tau: f64,
    rho: f64,
}

impl Gaussian {
    pub fn new(tau: f64, rho: f64) -> Self {
        Self { tau, rho }
    }

    pub fn from_mean_and_variance(mean: f64, variance: f64) -> Self {
        let tau = mean / variance;
        let rho = 1.0 / variance;
        Self { tau, rho }
    }
}

impl Gaussian {
    pub fn tau(&self) -> f64 {
        self.tau
    }

    pub fn rho(&self) -> f64 {
        self.rho
    }

    pub fn mean(&self) -> f64 {
        self.tau / self.rho
    }

    pub fn variance(&self) -> f64 {
        1.0 / self.rho
    }
}

impl Gaussian {
    pub fn absdiff(&self, other: &Self) -> f64 {
        let tau_diff = (self.tau - other.tau).abs();
        let rho_diff = (self.rho - other.rho).abs().sqrt();
        tau_diff.max(rho_diff)
    }
}

impl Mul<Gaussian> for Gaussian {
    type Output = Gaussian;

    fn mul(self, other: Gaussian) -> Gaussian {
        Gaussian {
            tau: self.tau + other.tau,
            rho: self.rho + other.rho,
        }
    }
}

impl Div<Gaussian> for Gaussian {
    type Output = Gaussian;

    fn div(self, other: Gaussian) -> Gaussian {
        Gaussian {
            tau: self.tau - other.tau,
            rho: self.rho - other.rho,
        }
    }
}

impl Default for Gaussian {
    fn default() -> Self {
        Gaussian { tau: 0.0, rho: 0.0 }
    }
}
