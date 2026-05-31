use std::io::BufRead;

/// Represents a single parsed record from a FASTA file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    /// The unique sequence identifier (the first string token following the '>' symbol).
    pub id: String,
    /// Optional free-text description metadata appearing after the identifier.
    pub description: Option<String>,
    /// The raw biological sequence characters, stripped of line breaks and structural whitespace.
    pub sequence: String,
}

/// A zero-copy overhead, memory-efficient streaming Multi-FASTA parser.
pub struct FastaReader<R: BufRead> {
    reader: R,
    scratch_buf: String,
    next_header: Option<String>,
}

impl<R: BufRead> FastaReader<R> {
    /// Instantiates a new streaming parser from any input source implementing `BufRead`.
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            scratch_buf: String::new(),
            next_header: None,
        }
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = Result<FastaRecord, std::io::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        // Resolve the header line for this iteration block
        let header = match self.next_header.take() {
            Some(cached_header) => cached_header,
            None => loop {
                self.scratch_buf.clear();
                match self.reader.read_line(&mut self.scratch_buf) {
                    Ok(0) => return None, // Clean EOF encountered before a new record starts
                    Ok(_) => {
                        let trimmed = self.scratch_buf.trim();
                        if trimmed.is_empty() {
                            continue; // Skip safely over sporadic leading or separating blank lines
                        }
                        if trimmed.starts_with('>') {
                            break trimmed.to_string();
                        }
                        return Some(Err(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            "Malformed FASTA record: Found sequence payload characters prior to an initial sequence header line ('>')",
                        )));
                    }
                    Err(err) => return Some(Err(err)),
                }
            },
        };

        // Parse identity key and free-text description from the header metadata
        let header_payload = &header[1..]; // Strip out structural leading sequence marker '>'
        let mut tokens = header_payload.splitn(2, |c: char| c.is_whitespace());
        let id = tokens.next().unwrap_or("").to_string();
        let description = tokens
            .next()
            .map(|desc| desc.trim().to_string())
            .filter(|desc| !desc.is_empty());

        //  Stream and accumulate downstream line data into the sequence buffer
        let mut sequence = String::new();
        loop {
            self.scratch_buf.clear();
            match self.reader.read_line(&mut self.scratch_buf) {
                Ok(0) => break, // Reaching EOF here terminates this final sequence payload safely
                Ok(_) => {
                    let trimmed = self.scratch_buf.trim();
                    if trimmed.is_empty() {
                        continue; // Skip isolated empty internal rows
                    }
                    if trimmed.starts_with('>') {
                        // Lookahead protection: Cache this line for the subsequent next() evaluation pass
                        self.next_header = Some(trimmed.to_string());
                        break;
                    }
                    sequence.push_str(trimmed);
                }
                Err(err) => return Some(Err(err)),
            }
        }

        Some(Ok(FastaRecord {
            id,
            description,
            sequence,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    #[test]
    fn test_valid_standard_multi_fasta() {
        let fasta_payload = "
>ID_1234 Extended descriptive metadata field
GATTACA
GATTACA
>ID_5678
AAAAA
        ";

        let mut parser = FastaReader::new(BufReader::new(fasta_payload.as_bytes()));

        let rec_1 = parser.next().unwrap().unwrap();
        assert_eq!(rec_1.id, "ID_1234");
        assert_eq!(
            rec_1.description,
            Some("Extended descriptive metadata field".to_string())
        );
        assert_eq!(rec_1.sequence, "GATTACAGATTACA");

        let rec_2 = parser.next().unwrap().unwrap();
        assert_eq!(rec_2.id, "ID_5678");
        assert_eq!(rec_2.description, None);
        assert_eq!(rec_2.sequence, "AAAAA");

        assert!(parser.next().is_none());
    }

    #[test]
    fn test_invalid_fasta_sequence_first() {
        let broken_payload = "ATGCTAGCTAG\n>Rosalind_1111\nATGC";
        let mut parser = FastaReader::new(BufReader::new(broken_payload.as_bytes()));

        let result = parser.next().unwrap();
        assert!(result.is_err());
        if let Err(e) = result {
            assert_eq!(e.kind(), std::io::ErrorKind::InvalidData);
        }
    }

    #[test]
    fn test_empty_input_returns_none() {
        let empty_payload = "   \n   \n";
        let mut parser = FastaReader::new(BufReader::new(empty_payload.as_bytes()));
        assert!(parser.next().is_none());
    }
}
