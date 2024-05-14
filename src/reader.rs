use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use fxhash::FxHashMap;
use memchr::memchr;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::{fmt::Debug, fs::File, io::Read, path::Path};

use crate::chain::Chain;
use crate::chainmap::ChainMap;

/// A reader for chain files.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Reader;

impl Reader {
    /// Create a new reader.
    ///
    /// # Arguments
    ///
    /// * `file` - A path to a chain file.
    ///
    /// # Returns
    ///
    /// A `Result` containing a `FxHashMap` of `Chain` objects.
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use chainder as chain;
    ///
    /// let data = chain::Reader::from_file("/path/to/chainfile")?;
    /// ```
    pub fn from_file<T>(file: T) -> Result<ChainMap>
    where
        T: AsRef<Path> + Debug,
    {
        let data = Self::open(file)?;
        Self::parse(&data)
    }

    /// Parser for chain files.
    ///
    /// # Arguments
    /// * `data` - A reference to a byte slice.
    ///
    /// # Returns
    /// A `Result` containing a `FxHashMap` of `Chain` objects.
    ///
    /// # Example
    ///
    /// ```
    /// use chainder as chain;
    ///
    /// let data = chain::Reader::parse(&data)?;
    /// ```
    pub fn parse(data: &[u8]) -> Result<ChainMap> {
        let mut vacc: Vec<(&[u8], &[u8])> = Vec::new();
        let mut data = &data[..];
        loop {
            let sep = memchr(b'\n', &data).with_context(|| {
                format!(
                    "Failed to find separator in: {:?}. Bad formatted line!",
                    String::from_utf8_lossy(data)
                )
            })?;
            let Some(end) = memchr(b'c', &data[sep..]) else {
                let header = &data[..sep];
                let block = &data[sep + 1..];
                vacc.push((header, block));
                break;
            };
            let header = &data[..sep];
            let block = &data[sep + 1..sep + end - 1];
            vacc.push((header, block));
            data = &data[sep + end..];
        }

        let chainfile = vacc
            .par_iter()
            .filter_map(|(header, block)| Chain::from(header, block).ok())
            .fold(
                || FxHashMap::default(),
                |mut acc, chain| {
                    acc.insert(chain.0, chain.1);
                    acc
                },
            )
            .reduce(
                || FxHashMap::default(),
                |mut acc, map| {
                    acc.extend(map);
                    acc
                },
            );

        Ok(ChainMap { map: chainfile })
    }

    /// Create a new reader from a byte slice.
    ///
    /// # Arguments
    ///
    /// * `data` - A reference to a byte slice.
    ///
    /// # Returns
    ///
    /// A `Result` containing a `FxHashMap` of `Chain` objects.
    ///
    /// # Example
    ///
    /// ```
    /// use chainder as chain;
    ///
    /// let line = b"chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1\n9\t1\t0\n\n";
    /// let data = chain::Reader::from_bytes(line)?;
    ///
    /// assert_eq!(data.len(), 1);
    /// ```
    pub fn from_bytes(data: &[u8]) -> Result<ChainMap> {
        Reader::parse(data)
    }

    /// Create a new reader from a binary file.
    ///
    /// # Arguments
    /// * `bin` - A path to a binary file.
    ///
    /// # Returns
    /// A `Result` containing a `FxHashMap` of `Chain` objects.
    ///
    /// # Example
    ///
    /// ```
    /// use chainder as chain;
    ///
    /// let data = chain::Reader::from_bin("/path/to/binfile")?;
    /// ```
    pub fn from_bin<T>(bin: T) -> Result<ChainMap>
    where
        T: AsRef<Path> + Debug,
    {
        let data = Self::open(bin)?;
        let decoded: ChainMap =
            bincode::deserialize(&data).with_context(|| "Deserialization failed")?;
        Ok(decoded)
    }

    /// Load a single chain from a binary file.
    ///
    /// # Arguments
    /// * `id` - A chain id.
    /// * `bin` - A path to a binary file.
    ///
    /// # Returns
    /// A `Result` containing a `Chain` object.
    ///
    /// # Example
    ///
    /// ```
    /// use chainder as chain;
    ///
    /// let data = chain::Reader::load_chain(&123, "/path/to/binfile")?;
    /// ```
    pub fn load_chain<T>(id: u32, bin: T) -> Result<Chain>
    where
        T: AsRef<Path> + Debug,
    {
        let data = Self::open(bin)?;
        let decoded: ChainMap =
            bincode::deserialize(&data).with_context(|| "Deserialization failed")?;
        Ok(decoded.get(&id).expect("Failed to get chain").clone())
    }

    /// Private opener for files.
    ///
    /// # Arguments
    /// * `file` - A path to a file.
    ///
    /// # Returns
    /// A `Result` containing a `Vec` of bytes.
    ///
    /// # Example
    ///
    /// ```
    /// use chainder as chain;
    ///
    /// let data = chain::Reader::open("/path/to/file")?;
    /// ```
    fn open<T>(file: T) -> Result<Vec<u8>>
    where
        T: AsRef<Path> + Debug,
    {
        let stat = std::fs::metadata(&file)
            .with_context(|| format!("Failed to get metadata for {:?}", file))?;
        let mut data = vec![];
        data.reserve(stat.len() as usize + 1);

        let mut f = File::open(&file).with_context(|| format!("Failed to open file {:?}", file))?;

        if file.as_ref().extension().unwrap() == "gz" {
            let mut decoder = MultiGzDecoder::new(f);
            decoder
                .read_to_end(&mut data)
                .with_context(|| format!("Failed to read file {:?}", file))?;
        } else {
            f.read_to_end(&mut data)
                .with_context(|| format!("Failed to read file {:?}", file))?;
        }

        Ok(data)
    }
}
