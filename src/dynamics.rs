/// Represents the demographic breakdown of a population at a specific point in time.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PopulationState {
    pub juveniles: u64,
    pub adults: u64,
}

impl PopulationState {
    pub fn new(juveniles: u64, adults: u64) -> Self {
        Self { juveniles, adults }
    }

    pub fn total(&self) -> u64 {
        self.juveniles + self.adults
    }
}

/// Represents a generalized discrete-time population projection model.
/// Used to simulate population dynamics, epidemiological growth, and demographic shifts.
#[derive(Debug, Clone)]
pub struct AgeStructuredModel {
    /// The number of offspring produced per reproductive individual per time step.
    pub fecundity: u64,
    /// The probability of an individual surviving to the next time step (0.0 to 1.0).
    /// Set to 1.0 for immortal models; lower values simulate mortality.
    pub survival_rate: f64,
}

impl AgeStructuredModel {
    pub fn new(fecundity: u64, survival_rate: f64) -> Self {
        Self {
            fecundity,
            survival_rate: survival_rate.clamp(0.0, 1.0),
        }
    }

    /// Projects the population state forward by a given number of elapsed time intervals,
    /// starting from an initial population distribution.
    pub fn project(
        &self,
        initial_state: PopulationState,
        elapsed_intervals: u32,
    ) -> PopulationState {
        let mut juveniles = initial_state.juveniles as f64;
        let mut adults = initial_state.adults as f64;

        // Shadow with the right format
        let fecundity = self.fecundity as f64;

        // If elapsed_intervals is 0, this loop does not execute,
        // and initial_state is returned perfectly intact.
        for _ in 1..=elapsed_intervals {
            let new_offspring = adults * fecundity;
            let surviving_juveniles = juveniles * self.survival_rate;
            let surviving_adults = adults * self.survival_rate;

            // Shift cohorts forward
            adults = surviving_adults + surviving_juveniles;
            juveniles = new_offspring;
        }

        PopulationState::new(juveniles.round() as u64, adults.round() as u64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_binary_fission_with_immortality() {
        // Standard binary fission / ideal growth (fecundity=1, survival=1.0)
        let model = AgeStructuredModel::new(1, 1.0);
        let initial = PopulationState::new(10, 0);
        let step_3 = model.project(initial, 3);
        assert_eq!(step_3.adults, 20);
        assert_eq!(step_3.juveniles, 10);
    }

    #[test]
    fn test_population_decay_with_mortality() {
        // No reproduction, 50% survival rate
        let model = AgeStructuredModel::new(0, 0.5);
        // Start with a substantial real-world population
        let initial = PopulationState::new(1000, 0);

        // Step 1: 1000 juveniles -> 50% survive to become 500 adults (0 juveniles born)
        let step_1 = model.project(initial, 1);
        assert_eq!(step_1.juveniles, 0);
        assert_eq!(step_1.adults, 500);

        // Step 2: 500 adults -> 50% survive to become 250 adults
        let step_2 = model.project(initial, 2);
        assert_eq!(step_2.total(), 250);
    }
}
