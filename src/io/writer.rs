use anyhow::Result;
use bincode;
use serde::Serialize;
use std::{fmt::Debug, fs::File, io::prelude::*, path::Path};

/// Write chaintools object to a file
pub struct Writer;

impl Writer {
    /// Encode and write chaintools object to a file
    ///
    /// # Arguments
    /// * `data` - A chaintools object
    ///
    /// # Returns
    /// * Result<()>
    ///
    /// # Example
    ///
    /// ```
    /// use chaintools as chain;
    ///
    ///
    /// let line = b"chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1\n9\t1\t0\n\n";
    /// let data = chain::Reader::from_bytes(line)?;
    /// chain::Writer::write_index("/path/to/output.bin", data).unwrap();
    /// ```
    pub fn to_bin<K, T>(data: K, path: T) -> Result<()>
    where
        K: Serialize + Clone + Debug + Sync + Send,
        T: AsRef<Path> + Debug,
    {
        let encoded: Vec<u8> = bincode::serialize(&data)?;
        let mut file = File::create(path)?;
        file.write_all(&encoded).expect("Failed to write to file");
        Ok(())
    }

    /// Encode, compress and write chaintools object to a file
    ///
    /// # Arguments
    /// * `data` - A chaintools object
    /// * `path` - A path to the output file
    ///
    /// # Returns
    /// * Result<()>
    ///
    /// # Example
    /// ```
    /// use chaintools as chain;
    ///
    /// let line = b"chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1\n9\t1\t0\n\n";
    /// let data = chain::Reader::from_bytes(line)?;
    /// chain::Writer::write_index_gz("/path/to/output.bin.gz", data).unwrap();
    /// ```
    ///
    pub fn to_bin_gz<K, T>(data: K, path: T) -> Result<()>
    where
        K: Serialize + Clone + Debug + Sync + Send,
        T: AsRef<Path> + Debug,
    {
        let encoded: Vec<u8> = bincode::serialize(&data)?;
        let mut file =
            flate2::write::GzEncoder::new(File::create(path)?, flate2::Compression::default());
        file.write_all(&encoded).expect("Failed to write to file");
        Ok(())
    }
}
