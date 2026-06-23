use crate::{gaussian::Gaussian, variable::Variable};

pub struct GaussianFactor {
    x: Variable,
    prior: Gaussian,
    message_to_x: Gaussian,
}

impl GaussianFactor {
    pub fn new(x: &Variable, prior: Gaussian) -> Self {
        Self {
            x: x.clone(),
            prior,
            message_to_x: Gaussian::default(),
        }
    }
}

impl GaussianFactor {
    pub fn update_message_to_x(&mut self) -> Result<f64, String> {
        let previous_message = self.message_to_x;
        let new_message = self.prior;

        let current_value = self.x.get();
        let updated = current_value / previous_message * new_message;

        if updated.rho() <= 0.0 {
            return Err(format!("Invalid Gaussian: rho = {} <= 0", updated.rho()));
        }

        self.x.set(updated);
        self.message_to_x = new_message;

        Ok(previous_message.absdiff(&new_message))
    }
}
