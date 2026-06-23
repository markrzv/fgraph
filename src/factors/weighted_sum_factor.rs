use crate::{gaussian::Gaussian, variable::Variable};

pub struct WeightedSumFactor {
    x: Variable,
    y: Variable,
    z: Variable,
    x_coefficient: f64,
    y_coefficient: f64,
    message_to_x: Gaussian,
    message_to_y: Gaussian,
    message_to_z: Gaussian,
}

impl WeightedSumFactor {
    pub fn new(
        x: &Variable,
        y: &Variable,
        z: &Variable,
        x_coefficient: f64,
        y_coefficient: f64,
    ) -> Self {
        Self {
            x: x.clone(),
            y: y.clone(),
            z: z.clone(),
            x_coefficient,
            y_coefficient,
            message_to_x: Gaussian::default(),
            message_to_y: Gaussian::default(),
            message_to_z: Gaussian::default(),
        }
    }
}

fn weighted_sum_message(
    source_one: Gaussian,
    source_one_message: Gaussian,
    source_one_coefficient: f64,
    source_two: Gaussian,
    source_two_message: Gaussian,
    source_two_coefficient: f64,
    target_marginal: Gaussian,
    target_message: Gaussian,
) -> Result<(Gaussian, Gaussian, f64), String> {
    let source_one_without_message = source_one / source_one_message;
    let source_two_without_message = source_two / source_two_message;

    if source_one_without_message.rho() < 0.0 || source_two_without_message.rho() < 0.0 {
        return Err("Invalid Gaussian in weighted_sum_message: negative rho".to_string());
    }

    if source_one_without_message.rho() == 0.0 || source_two_without_message.rho() == 0.0 {
        return Ok((target_marginal, target_message, 0.0));
    }

    let new_mean = source_one_coefficient * source_one_without_message.mean()
        + source_two_coefficient * source_two_without_message.mean();

    let new_variance = source_one_coefficient.powi(2) * source_one_without_message.variance()
        + source_two_coefficient.powi(2) * source_two_without_message.variance();

    let new_message = Gaussian::from_mean_and_variance(new_mean, new_variance);

    let new_marginal = (target_marginal / target_message) * new_message;

    if new_marginal.rho() <= 0.0 {
        return Err(format!(
            "Invalid Gaussian in weighted_sum_message: resulting marginal rho = {} <= 0",
            new_marginal.rho()
        ));
    }

    let difference = target_marginal.absdiff(&new_marginal);

    Ok((new_marginal, new_message, difference))
}

impl WeightedSumFactor {
    pub fn update_message_to_x(&mut self) -> Result<f64, String> {
        let (new_marginal, new_message, difference) = weighted_sum_message(
            self.z.get(),
            self.message_to_z,
            1.0 / self.x_coefficient,
            self.y.get(),
            self.message_to_y,
            -self.y_coefficient / self.x_coefficient,
            self.x.get(),
            self.message_to_x,
        )?;

        self.x.set(new_marginal);
        self.message_to_x = new_message;

        Ok(difference)
    }

    pub fn update_message_to_y(&mut self) -> Result<f64, String> {
        let (new_marginal, new_message, difference) = weighted_sum_message(
            self.z.get(),
            self.message_to_z,
            1.0 / self.y_coefficient,
            self.x.get(),
            self.message_to_x,
            -self.x_coefficient / self.y_coefficient,
            self.y.get(),
            self.message_to_y,
        )?;

        self.y.set(new_marginal);
        self.message_to_y = new_message;

        Ok(difference)
    }

    pub fn update_message_to_z(&mut self) -> Result<f64, String> {
        let (new_marginal, new_message, difference) = weighted_sum_message(
            self.x.get(),
            self.message_to_x,
            self.x_coefficient,
            self.y.get(),
            self.message_to_y,
            self.y_coefficient,
            self.z.get(),
            self.message_to_z,
        )?;

        self.z.set(new_marginal);
        self.message_to_z = new_message;

        Ok(difference)
    }
}
