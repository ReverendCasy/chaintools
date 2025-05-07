use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use fxhash::FxHashMap;
use memchr::{memchr, memrchr};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
// use std::borrow::Cow;
// use std::ops::RangeBounds;
use std::{fmt::{Debug, Display}, fs::File, io::{BufRead, BufReader, Read, Seek, SeekFrom}, path::Path};

use crate::cmap::chain::Chain;
use crate::cmap::map::ChainMap;
// use crate::io::indexer::BinaryIndex;

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
    /// use chaintools as chain;
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
    /// use chaintools as chain;
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
    /// use chaintools as chain;
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
    /// use chaintools as chain;
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
    /// use chaintools as chain;
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
    /// use chaintools as chain;
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

    /// Extract selected chains from the binary vector
    /// 
    /// # Arguments
    /// * `file` - A path to the chain file 
    /// 
    /// * `chains` - a vector of string literals containing the chain IDs to extract
    pub fn extract<U>(file: U, chains: Vec<&str>) -> Result<ChainMap> 
    where 
        U: AsRef<Path> + Debug,
        // T:  IntoIterator + RangeBounds<str>,
        // T::IntoIter: ExactSizeIterator
    {
        let mut chainmap: FxHashMap<u32, Chain> = FxHashMap::default();
        if chains.len() == 0 {
            return Ok(ChainMap { map: chainmap });
        }
        let data: Vec<u8>  = Self::open(file)?;
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
                let id_field: usize = memrchr(b' ', header).with_context(|| {
                    format!(
                        "Improperly formatted header line: {:?} !",
                        String::from_utf8_lossy(header)
                    )
                }
                )?;
                let id_str = std::str::from_utf8(&header[id_field+1..])?;
                if chains.contains(&id_str) {
                    let block = &data[sep + 1..];
                    let (chain_id, chain_obj) = Chain::from(header, block).unwrap();
                    chainmap.insert(chain_id, chain_obj);
                } 
                break;
            };
            let header = &data[..sep];
            let id_field: usize = memrchr(b' ', header).with_context(|| {
                format!(
                    "Improperly formatted header line: {:?} !",
                    String::from_utf8_lossy(header)
                )
            }
            )?;
            let id_str: &str = std::str::from_utf8(&header[id_field+1..])?;
            // add the chain if its ID is in the requested chains vector
            if chains.contains(&id_str) {
                let block = &data[sep + 1..sep + end - 1];
                let (chain_id, chain_obj) = Chain::from(header, block).unwrap();
                chainmap.insert(chain_id, chain_obj);
                // break once all the chains were extracted
                // TODO: Check for duplicate elements in the chain ID vector??
                if chainmap.len() == chains.len() {break};
            } 
            data = &data[sep + end..];
        }
        Ok(ChainMap { map: chainmap })
    }

    // Private function for index reading
    //
    // # Arguments
    // 
    // * `file` - an index file
    fn read_index<U>(file: U, chains: &Vec<u64>, all: bool) -> Result<FxHashMap<u64, (u64, u64)>>
    where
        U: AsRef<Path> + Debug + Display
    {
        let mut index: FxHashMap<u64, (u64, u64)> = FxHashMap::default();
        let f = File::open(file)?;
        let buf = BufReader::new(f).lines();
        for line in buf.map_while(Result::ok) {
            let line_data: Vec<u64> = line.split('\t')
                .map(|x|
                    x.parse::<u64>().expect("Invalid numeric value found in the index file")
                )
                .collect::<Vec<u64>>();
            if chains.contains(&line_data[0]) || all {
                index.insert(line_data[0], (line_data[1], line_data[2]));
            }
            if index.len() == chains.len() && !all {break}
        }
        Ok(index)
    }

    /// Extract selected chains from an indexed file
    /// 
    /// # Arguments
    /// 
    /// 
    /// 
    pub fn extract_ix<U>(file: U, chains: Vec<u64>, all: bool) -> Result<ChainMap>
    where
        U: AsRef<Path> + Debug + Display
    {
        let mut chainmap: FxHashMap<u32, Chain> = FxHashMap::default();
        if chains.len() == 0 {return Ok(ChainMap{ map: chainmap })}
        let index_file: String = format!("{}.ix", &file);
        // let index: FxHashMap<u64, (usize, usize)> = BinaryIndex::read_index(index_file)?;
        let index: FxHashMap<u64, (u64, u64)> = Self::read_index(index_file, &chains, all)?;

        let mut f = File::open(file)?;
        for chain_id in chains {
            // get chain coordinates, panic if the chain Id is missing from the index file
            let (start, end) = match index.get(&chain_id) {
                Some(x) => {(x.0 as u64, x.1 as u64)},
                None => {
                    panic!("Chain {} was not found in the index file", chain_id)
                }
            };
            // extract the chain bytes
            f.seek(SeekFrom::Start(start))?;
            let mut chain_string = vec![0; (end-start) as usize];
            f.read(&mut chain_string[..])?;
            // split chain into body and header
            let header_pos = memchr(b'\n', &chain_string).with_context(|| {
                format!(
                    "Failed to find newline separator in: {:?}",
                    String::from_utf8_lossy(&chain_string)
                )
            })?;
            let header: &[u8] = &chain_string[..header_pos];
            let block: &[u8] = &chain_string[header_pos+1..];
            let (id, chain) = Chain::from(header, block)?;
            chainmap.insert(id, chain);
        }
        Ok(ChainMap{ map: chainmap } )
    }
}
