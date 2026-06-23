use std::f64::consts::{PI, SQRT_2};

use libm::erf;

use crate::{gaussian::Gaussian, variable::Variable};

pub struct GreaterThanFactor {
    x: Variable,
    epsilon: f64,
    message_to_x: Gaussian,
}

impl GreaterThanFactor {
    pub fn new(x: &Variable, epsilon: f64) -> Self {
        Self {
            x: x.clone(),
            epsilon,
            message_to_x: Gaussian::default(),
        }
    }
}

impl GreaterThanFactor {
    pub fn update_message_to_x(&mut self) -> Result<f64, String> {
        let old_marginal = self.x.get();
        let old_message = self.message_to_x;

        let marginal_without_message = old_marginal / old_message;

        if marginal_without_message.rho() <= 0.0 {
            return Err(format!(
                "Invalid Gaussian in update_message_to_x: rho = {} <= 0",
                marginal_without_message.rho()
            ));
        }

        let mean_value = marginal_without_message.mean();
        let variance_value = marginal_without_message.variance();
        let standard_deviation = variance_value.sqrt();

        let normalized_mean = (mean_value - self.epsilon) / standard_deviation;

        let truncated_mean = mean_value + standard_deviation * v_function(normalized_mean);
        let truncated_variance = variance_value * (1.0 - w_function(normalized_mean));

        let new_marginal = Gaussian::from_mean_and_variance(truncated_mean, truncated_variance);
        let new_message = new_marginal / marginal_without_message;

        self.x.set(new_marginal);
        self.message_to_x = new_message;

        Ok(old_marginal.absdiff(&new_marginal))
    }
}

fn phi(x: f64) -> f64 {
    (-0.5 * x * x).exp() / (2.0 * PI).sqrt()
}

fn big_phi(x: f64) -> f64 {
    0.5 * (1.0 + erf(x / SQRT_2))
}

fn v_function(x: f64) -> f64 {
    let denom = big_phi(x);
    if denom < 1e-10 { -x } else { phi(x) / denom }
}

fn w_function(x: f64) -> f64 {
    let v = v_function(x);
    v * (v + x)
}
