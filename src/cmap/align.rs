use anyhow::{Context, Result};
use memchr::memchr;
use serde::{Deserialize, Serialize};
use std::str::from_utf8;

/// A structure to represent an alignment record.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentRecord {
    pub size: u32,
    pub dt: u32,
    pub dq: u32,
    pub is_last: bool
}

impl AlignmentRecord {
    /// Parse a byte array into a vector of alignment records.
    ///
    /// # Arguments
    /// * `align` - A byte array
    ///
    /// # Returns
    /// * Vec<AlignmentRecord>
    ///
    /// # Example
    ///
    /// ```
    /// use chaintools as chain;
    ///
    /// let line = b"9\t1\t0\n";
    /// let data = chain::AlignmentRecord::parse(line);
    /// println!("{:?}", data);
    ///
    /// > [AlignmentRecord { size: 9, dt: 1, dq: 0, is_last: false}];
    /// ```
    pub fn parse(align: &[u8]) -> Vec<AlignmentRecord> {
        Self::parse_byte(align).expect("ERROR: Failed to parse alignment record")
    }

    /// Parse a byte array from a single alignment line into an alignment record.    
    ///
    /// # Arguments
    /// * `align` - A byte array
    ///
    /// # Returns
    /// * Option<AlignmentRecord>
    ///
    /// # Example
    ///
    /// ```
    /// use chaintools as chain;
    ///
    /// let line = b"9\t1\t0\n";
    /// let data = chain::AlignmentRecord::from(line);
    /// println!("{:?}", data);
    ///
    /// > Some(AlignmentRecord { size: 9, dt: 1, dq: 0, is_last: false });
    /// ```
    pub fn from(align: &[u8]) -> Option<Self> {
        Self::parse_byte(align)
            .expect("Failed to parse alignment record")
            .first()
            .cloned()
    }

    /// Create a new alignment record.
    ///
    /// # Arguments
    /// * `size` - The size of the ungapped alignment
    /// * `dt` -  The difference between the end of this block and the beginning of the next block (reference/target sequence)
    /// * `dq` - The difference between the end of this block and the beginning of the next block (query sequence)
    ///
    /// # Returns
    /// * AlignmentRecord
    ///
    /// # Example
    ///
    /// ```
    /// use chaintools as chain;
    ///
    /// let data = chain::AlignmentRecord::new(9, 1, 0);
    /// println!("{:?}", data);
    ///
    /// > AlignmentRecord { size: 9, dt: 1, dq: 0 };
    /// ```
    pub fn new(size: u32, dt: u32, dq: u32, is_last: bool) -> Self {
        Self { size, dt, dq, is_last}
    }

    /// Dummy getter for the size of the ungapped alignment.
    ///
    /// # Arguments
    /// * `self` - An alignment record
    ///
    /// # Returns
    /// * u32
    ///
    /// # Example
    ///
    /// ```
    /// use chaintools as chain;
    ///
    /// let data = chain::AlignmentRecord::new(9, 1, 0, false);
    /// let line = data.size();
    /// println!("{:?}", line);
    ///
    /// > 9
    /// ```
    pub fn size(&self) -> u32 {
        self.size
    }

    /// Dummy getter for dt
    ///
    /// # Arguments
    /// * `self` - An alignment record
    ///
    /// # Returns
    /// * u32
    ///
    /// # Example
    ///
    /// ```
    /// use chaintools as chain;
    ///
    /// let data = chain::AlignmentRecord::new(9, 1, 0, false);
    /// let line = data.dt();
    /// println!("{:?}", line);
    ///
    /// > 1
    /// ```
    pub fn dt(&self) -> u32 {
        self.dt
    }

    /// Dummy getter for dq
    ///
    /// # Arguments
    /// * `self` - An alignment record
    ///
    /// # Returns
    /// * u32
    ///
    /// # Example
    ///
    /// ```
    /// use chaintools as chain;
    ///
    /// let data = chain::AlignmentRecord::new(9, 1, 0, false);
    /// let line = data.dq();
    /// println!("{:?}", line);
    ///
    /// > 0
    /// ```
    pub fn dq(&self) -> u32 {
        self.dq
    }

    /// Get the alignment record as a string.
    ///
    /// # Arguments
    /// * `self` - An alignment record
    ///
    /// # Returns
    /// * String
    ///
    /// # Example
    ///
    /// ```
    /// use chaintools as chain;
    ///
    /// let data = chain::AlignmentRecord::new(9, 1, 0, false);
    /// let line = data.to_string();
    /// println!("{:?}", line);
    ///
    /// > "9\t1\t0"
    /// ```
    pub fn to_string(&self) -> String {
        format!("{}\t{}\t{}", self.size, self.dt, self.dq)
    }

    /// Parse a byte array into a vector of alignment records.
    ///
    /// # Arguments
    /// * `align` - A byte array
    ///
    /// # Returns
    /// * Result<Vec<AlignmentRecord>, Box<dyn std::error::Error>>
    ///
    /// # Example
    /// ```
    /// use chaintools as chain;
    ///
    /// let line = b"9\t1\t0\n";
    /// let data = chain::AlignmentRecord::parse_byte(line);
    /// println!("{:?}", data);
    ///
    /// > Ok([AlignmentRecord { size: 9, dt: 1, dq: 0, is_last: false}]);
    /// ```
    fn parse_byte(align: &[u8]) -> Result<Vec<AlignmentRecord>> {
        let mut acc = vec![];
        let mut align = &align[..];

        loop {
            let Some(sep) = memchr(b'\t', align) else {
                if align.len() < 2 {
                    break;
                } // YM: Does this mean that the last block does not get recorded?

                let end = memchr(b'\n', &align[..]).with_context(|| {
                    format!(
                        "Failed to find separator in: {:?}. Bad formatted line!",
                        String::from_utf8_lossy(align)
                    )
                })?;

                let size = from_utf8(&align[..end])
                    .context("Failed to parse size")
                    .and_then(|s| s.parse::<u32>().context("Failed to parse size"))?;

                acc.push(AlignmentRecord {
                    size: size,
                    dt: 0,
                    dq: 0,
                    is_last: true
                });
                break;
            };

            let end = memchr(b'\n', &align[sep..]).unwrap();
            let mid = memchr(b'\t', &align[sep + 1..]).unwrap();

            let size = from_utf8(&align[..sep])?.parse::<u32>().with_context(|| {
                format!("Failed to parse size: {:?}", String::from_utf8_lossy(align))
            })?;

            let dt = from_utf8(&align[sep + 1..sep + mid + 1])?
                .parse::<u32>()
                .with_context(|| {
                    format!(
                        "Failed to parse dt from slice: {:?}",
                        String::from_utf8_lossy(&align[sep + 1..])
                    )
                })?;
            let dq = from_utf8(&align[sep + mid + 2..sep + end])?
                .parse::<u32>()
                .with_context(|| {
                    format!(
                        "Failed to parse dq: {:?}",
                        String::from_utf8_lossy(&align[sep + mid + 2..])
                    )
                })?;

            acc.push(AlignmentRecord {
                size: size,
                dt: dt,
                dq: dq,
                is_last: false
            });

            align = &align[sep + end + 1..];
        }

        Ok(acc)
    }
}
