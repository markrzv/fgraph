use fgraph::{
    factors::{GaussianFactor, GaussianMeanFactor, GreaterThanFactor, WeightedSumFactor},
    gaussian::Gaussian,
    variable::Variable,
};

fn approx(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-6
}

#[test]
fn gaussian_factor_update() {
    let s = Variable::default();
    assert_eq!(s.get().tau(), 0.0);
    assert_eq!(s.get().rho(), 0.0);

    let mut factor1 = GaussianFactor::new(&s, Gaussian::from_mean_and_variance(2.0, 42.0));
    let mut factor2 = GaussianFactor::new(&s, Gaussian::from_mean_and_variance(1.0, 1.0));

    assert!(approx(
        factor1.update_message_to_x().unwrap(),
        0.1543033499620919
    ));
    assert!(approx(factor1.update_message_to_x().unwrap(), 0.0));
    assert!(approx(s.get().mean(), 2.0));
    assert!(approx(s.get().variance(), 42.0));

    assert!(approx(factor2.update_message_to_x().unwrap(), 1.0));
    assert!(approx(s.get().mean(), 1.0232558139534884));
    assert!(approx(s.get().variance(), 0.9767441860465117));
}

#[test]
fn gaussian_mean_factor_update() {
    let s1 = Variable::default();
    let s2 = Variable::default();

    let mut factor1 = GaussianFactor::new(&s1, Gaussian::from_mean_and_variance(3.0, 1.0));
    let mut factor2 = GaussianMeanFactor::new(&s2, &s1, 0.5);

    assert!(approx(factor1.update_message_to_x().unwrap(), 3.0));
    assert!(approx(factor2.update_message_to_x().unwrap(), 2.0));

    assert!(approx(s1.get().mean(), 3.0));
    assert!(approx(s2.get().mean(), 3.0));
    assert!(approx(s1.get().variance(), 1.0));
    assert!(approx(s2.get().variance(), 1.5));

    assert!(approx(factor2.update_message_to_y().unwrap(), 0.0));
}

#[test]
fn weighted_sum_factor_update() {
    let s1 = Variable::default();
    let s2 = Variable::default();
    let s3 = Variable::default();

    let mut f1 = GaussianFactor::new(&s1, Gaussian::from_mean_and_variance(1.0, 1.0));
    let mut f2 = GaussianFactor::new(&s2, Gaussian::from_mean_and_variance(2.0, 4.0));
    let mut f3 = GaussianFactor::new(&s3, Gaussian::from_mean_and_variance(2.0, 0.5));
    let mut weighted_sum_factor = WeightedSumFactor::new(&s1, &s2, &s3, 0.5, 0.5);

    assert!(approx(f1.update_message_to_x().unwrap(), 1.0));
    assert!(approx(f2.update_message_to_x().unwrap(), 0.5));

    assert!(approx(
        weighted_sum_factor.update_message_to_z().unwrap(),
        1.2
    ));
    assert!(approx(s3.get().mean(), 1.5));
    assert!(approx(s3.get().variance(), 1.25));

    assert!(approx(f3.update_message_to_x().unwrap(), 4.0));
    assert!(approx(s3.get().mean(), 1.8571428571428574));
    assert!(approx(s3.get().variance(), 0.35714285714285715));

    assert!(approx(s1.get().mean(), 1.0));
    assert!(approx(s1.get().variance(), 1.0));
    assert!(approx(
        weighted_sum_factor.update_message_to_x().unwrap(),
        0.40824829046386313
    ));
    assert!(approx(s1.get().mean(), 1.142857142857143));
    assert!(approx(s1.get().variance(), 0.8571428571428571));

    assert!(approx(s2.get().mean(), 2.0));
    assert!(approx(s2.get().variance(), 4.0));
    assert!(approx(
        weighted_sum_factor.update_message_to_y().unwrap(),
        1.0000000000000002
    ));
    assert!(approx(s2.get().mean(), 2.571428571428572));
    assert!(approx(s2.get().variance(), 1.7142857142857144));
}

#[test]
fn greater_than_factor_update() {
    let s = Variable::default();

    let mut f = GaussianFactor::new(&s, Gaussian::from_mean_and_variance(1.0, 1.0));
    let mut greater_than_factor = GreaterThanFactor::new(&s, 0.0);

    assert!(approx(f.update_message_to_x().unwrap(), 1.0));
    assert!(approx(s.get().mean(), 1.0));
    assert!(approx(s.get().variance(), 1.0));

    assert!(approx(
        greater_than_factor.update_message_to_x().unwrap(),
        1.0448277182202785
    ));
    assert!(approx(s.get().mean(), 1.2875999709391783));
    assert!(approx(s.get().variance(), 0.6296862857766055));
}
