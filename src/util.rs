use std::fmt;
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RnaBase {
    A, C, G, U,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DnaBase {
    A, C, G, T,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AminoAcid {
    A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V,
}

#[derive(Debug, PartialEq, Eq)]
pub struct InvalidSymbolError(pub char);

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DnaSequence(pub Vec<DnaBase>);

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RnaSequence(pub Vec<RnaBase>);

impl fmt::Display for DnaBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let c = match self {
            DnaBase::A => 'A', DnaBase::C => 'C', DnaBase::G => 'G', DnaBase::T => 'T',
        };
        write!(f, "{}", c)
    }
}

impl fmt::Display for RnaBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let c = match self {
            RnaBase::A => 'A', RnaBase::C => 'C', RnaBase::G => 'G', RnaBase::U => 'U',
        };
        write!(f, "{}", c)
    }
}

impl TryFrom<char> for DnaBase {
    type Error = InvalidSymbolError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c.to_ascii_uppercase() {
            'A' => Ok(DnaBase::A),
            'C' => Ok(DnaBase::C),
            'G' => Ok(DnaBase::G),
            'T' => Ok(DnaBase::T),
            _ => Err(InvalidSymbolError(c)),
        }
    }
}

impl TryFrom<char> for RnaBase {
    type Error = InvalidSymbolError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c.to_ascii_uppercase() {
            'A' => Ok(RnaBase::A),
            'C' => Ok(RnaBase::C),
            'G' => Ok(RnaBase::G),
            'U' => Ok(RnaBase::U),
            _ => Err(InvalidSymbolError(c)),
        }
    }
}

impl FromStr for DnaSequence {
    type Err = InvalidSymbolError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let bases = s
            .chars()
            .map(|c| match c.to_ascii_uppercase() {
                'A' => Ok(DnaBase::A),
                'C' => Ok(DnaBase::C),
                'G' => Ok(DnaBase::G),
                'T' => Ok(DnaBase::T),
                _ => Err(InvalidSymbolError(c)),
            })
            .collect::<Result<Vec<DnaBase>, _>>()?;
        Ok(DnaSequence(bases))
    }
}

impl fmt::Display for DnaSequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for base in &self.0 {
            let c = match base {
                DnaBase::A => 'A', DnaBase::C => 'C', DnaBase::G => 'G', DnaBase::T => 'T',
            };
            write!(f, "{}", c)?;
        }
        Ok(())
    }
}

impl FromStr for RnaSequence {
    type Err = InvalidSymbolError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let bases = s
            .chars()
            .map(|c| match c.to_ascii_uppercase() {
                'A' => Ok(RnaBase::A),
                'C' => Ok(RnaBase::C),
                'G' => Ok(RnaBase::G),
                'U' => Ok(RnaBase::U),
                _ => Err(InvalidSymbolError(c)),
            })
            .collect::<Result<Vec<RnaBase>, _>>()?;
        Ok(RnaSequence(bases))
    }
}

impl fmt::Display for RnaSequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for base in &self.0 {
            let c = match base {
                RnaBase::A => 'A', RnaBase::C => 'C', RnaBase::G => 'G', RnaBase::U => 'U',
            };
            write!(f, "{}", c)?;
        }
        Ok(())
    }
}

fn transcribe_dna_to_rna(dna_sequence: &[DnaBase]) -> Vec<RnaBase> {
    dna_sequence
        .iter()
        .map(|base| match base {
            DnaBase::A => RnaBase::A,
            DnaBase::C => RnaBase::C,
            DnaBase::G => RnaBase::G,
            DnaBase::T => RnaBase::U,
        })
    .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dna_parsing_valid() {
        let input = "AGCTagct";
        let expected = vec![
            DnaBase::A, DnaBase::G, DnaBase::C, DnaBase::T,
            DnaBase::A, DnaBase::G, DnaBase::C, DnaBase::T,
        ];

        let parsed: Vec<DnaBase> = input
            .chars()
            .map(DnaBase::try_from)
            .collect::<Result<_, _>>()
            .expect("Valid DNA characters should parse cleanly");

        assert_eq!(parsed, expected);
    }

    #[test]
    fn test_dna_parsing_invalid() {
        let input = "GATXG";
        let result: Result<Vec<DnaBase>, _> = input.chars().map(DnaBase::try_from).collect();

        assert_eq!(result, Err(InvalidSymbolError('X')));
    }

    #[test]
    fn test_dna_to_rna_transcription() {
        let dna = vec![DnaBase::G, DnaBase::A, DnaBase::T, DnaBase::C];
        let expected_rna = vec![RnaBase::G, RnaBase::A, RnaBase::U, RnaBase::C];

        assert_eq!(transcribe_dna_to_rna(&dna), expected_rna);
    }
}

#[cfg(test)]
mod serialization_tests {
    use super::*;

    #[test]
    fn test_dna_sequence_serialization_round_trip() {
        let original_raw = "GATGGAACTTGACTACGTAAATT";

        // 1. Deserialize (String -> Struct)
        let sequence: DnaSequence = original_raw.parse()
            .expect("Valid DNA string should deserialize seamlessly");

        // Verify internal structural representation
        assert_eq!(sequence.0[0], DnaBase::G);
        assert_eq!(sequence.0[2], DnaBase::T);

        // 2. Serialize (Struct -> String)
        let serialized_output = sequence.to_string();

        // Assert perfect round-trip parity
        assert_eq!(serialized_output, original_raw);
    }

    #[test]
    fn test_rna_sequence_serialization_round_trip() {
        let original_raw = "GAUGGAACUUGACUACGUAAAUU";

        // 1. Deserialize
        let sequence: RnaSequence = original_raw.parse()
            .expect("Valid RNA string should deserialize seamlessly");

        // 2. Serialize
        let serialized_output = sequence.to_string();

        assert_eq!(serialized_output, original_raw);
    }

    #[test]
    fn test_deserialization_failure_on_invalid_text() {
        // 'X' is invalid in DNA
        let corrupt_dna = "GATGXAACTT";
        let result: Result<DnaSequence, _> = corrupt_dna.parse();
        assert_eq!(result, Err(InvalidSymbolError('X')));

        // 'T' is invalid in RNA
        let corrupt_rna = "GAUGUAACTT";
        let rna_result: Result<RnaSequence, _> = corrupt_rna.parse();
        assert_eq!(rna_result, Err(InvalidSymbolError('T')));
    }
}
