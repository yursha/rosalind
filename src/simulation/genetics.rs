/// Represents a population distribution of Mendelian genotypes for a single factor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AllelePopulation {
    pub homozygous_dominant: u64,
    pub heterozygous: u64,
    pub homozygous_recessive: u64,
}

impl AllelePopulation {
    pub fn new(homozygous_dominant: u64, heterozygous: u64, homozygous_recessive: u64) -> Self {
        Self {
            homozygous_dominant,
            heterozygous,
            homozygous_recessive,
        }
    }

    /// Calculates the probability that two randomly selected organisms
    /// will produce offspring possessing at least one dominant allele.
    pub fn dominant_phenotype_probability(&self) -> f64 {
        let dom_homozygotes = self.homozygous_dominant as f64;
        let heterozygotes = self.heterozygous as f64;
        let rec_homozygotes = self.homozygous_recessive as f64;

        let total_population = dom_homozygotes + heterozygotes + rec_homozygotes;
        if total_population < 2.0 {
            return 0.0; // Avoid division by zero or mating a single organism
        }

        let num_ordered_pairs = total_population * (total_population - 1.0);

        // Calculate the probability of picking two parents that produce a homozygous recessive (aa) child:

        // 1. Both parents are homozygous recessive (aa + aa -> 100% aa offspring)
        let p_rec_from_rec_rec = (rec_homozygotes * (rec_homozygotes - 1.0)) / num_ordered_pairs;

        // 2. One parent is heterozygous, one is homozygous recessive (Aa + aa -> 50% aa offspring)
        // Note: this captures events "Aa + aa" and "aa + Aa" so we multiply it by 2
        let p_rec_from_het_rec = 2.0 * 0.5 * (heterozygotes * rec_homozygotes) / num_ordered_pairs;

        // 3. Both parents are heterozygous (Aa + Aa -> 25% aa offspring)
        let p_rec_from_het_het = 0.25 * (heterozygotes * (heterozygotes - 1.0)) / num_ordered_pairs;

        let p_rec = p_rec_from_rec_rec + p_rec_from_het_rec + p_rec_from_het_het;

        // Complement rule: P(at least one dominant) = 1 - P(both recessive)
        1.0 - p_rec
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper assertion function to handle floating-point precision tolerances
    fn assert_near(actual: f64, expected: f64) {
        let epsilon = 1e-5;
        assert!(
            (actual - expected).abs() < epsilon,
            "Expected {}, got {}",
            expected,
            actual
        );
    }

    #[test]
    fn test_equal_distribution() {
        let population = AllelePopulation::new(2, 2, 2);
        let probability = population.dominant_phenotype_probability();

        assert_near(probability, 0.78333);
    }

    #[test]
    fn test_insufficient_population() {
        // Population of 1 cannot mate
        let population = AllelePopulation::new(1, 0, 0);
        assert_near(population.dominant_phenotype_probability(), 0.0);

        // Empty population
        let empty_population = AllelePopulation::new(0, 0, 0);
        assert_near(empty_population.dominant_phenotype_probability(), 0.0);
    }

    #[test]
    fn test_pure_homozygous_dominant() {
        // If everyone is AA, probability of dominant offspring must be 1.0
        let population = AllelePopulation::new(10, 0, 0);
        assert_near(population.dominant_phenotype_probability(), 1.0);
    }

    #[test]
    fn test_pure_homozygous_recessive() {
        // If everyone is aa, probability of dominant offspring must be 0.0
        let population = AllelePopulation::new(0, 0, 10);
        assert_near(population.dominant_phenotype_probability(), 0.0);
    }

    #[test]
    fn test_dominant_and_recessive_only() {
        // AA (100) and aa (100).
        let population = AllelePopulation::new(100, 0, 100);
        assert_near(population.dominant_phenotype_probability(), 0.75126);
    }
}
