use crate::{gaussian::Gaussian, variable::Variable};

pub struct GaussianMeanFactor {
    x: Variable,
    y: Variable,
    beta_squared: f64,
    message_to_x: Gaussian,
    message_to_y: Gaussian,
}

impl GaussianMeanFactor {
    pub fn new(x: &Variable, y: &Variable, beta_squared: f64) -> Self {
        Self {
            x: x.clone(),
            y: y.clone(),
            beta_squared,
            message_to_x: Gaussian::default(),
            message_to_y: Gaussian::default(),
        }
    }
}

fn mean_factor_message(
    source: Gaussian,
    source_message: Gaussian,
    target: Gaussian,
    target_message: Gaussian,
    beta_squared: f64,
) -> Result<(Gaussian, Gaussian, f64), String> {
    let source_without_message = source / source_message;

    if source_without_message.rho() < 0.0 {
        return Err(format!(
            "Invalid Gaussian in mean_factor_message: rho = {} < 0",
            source_without_message.rho()
        ));
    }

    if source_without_message.rho() == 0.0 {
        return Ok((target, target_message, 0.0));
    }

    let new_message = Gaussian::from_mean_and_variance(
        source_without_message.mean(),
        source_without_message.variance() + beta_squared,
    );

    let new_marginal = target / target_message * new_message;

    if new_marginal.rho() <= 0.0 {
        return Err(format!(
            "Invalid Gaussian in mean_factor_message: resulting marginal rho = {} <= 0",
            new_marginal.rho()
        ));
    }

    let difference = target.absdiff(&new_marginal);

    Ok((new_marginal, new_message, difference))
}

impl GaussianMeanFactor {
    pub fn update_message_to_x(&mut self) -> Result<f64, String> {
        let (new_marginal, new_message, difference) = mean_factor_message(
            self.y.get(),
            self.message_to_y,
            self.x.get(),
            self.message_to_x,
            self.beta_squared,
        )?;

        self.x.set(new_marginal);
        self.message_to_x = new_message;

        Ok(difference)
    }

    pub fn update_message_to_y(&mut self) -> Result<f64, String> {
        let (new_marginal, new_message, difference) = mean_factor_message(
            self.x.get(),
            self.message_to_x,
            self.y.get(),
            self.message_to_y,
            self.beta_squared,
        )?;

        self.y.set(new_marginal);
        self.message_to_y = new_message;

        Ok(difference)
    }
}
