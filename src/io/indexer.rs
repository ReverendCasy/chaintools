use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use fxhash::FxHashMap;
use memchr::{memchr, memrchr};
use std::{fmt::{Debug, Display}, fs::File,  io::{Read, Write}, path::Path};

// use crate::io::writer::Writer;

/// Binary indexer for chain file

#[derive(Debug)]
pub struct BinaryIndex;

impl BinaryIndex{

    /// Private opener for files. Doubles the Reader::open function
    ///
    /// # Arguments
    /// * `file` - A path to a file.
    ///
    /// # Returns
    /// A `Result` containing a `Vec` of bytes.
    ///
    /// # Example
    ///
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

    /// Creates a binary index of a file
    /// 
    /// # Arguments
    /// 
    /// * `file` - A path to a chain file
    /// 
     /// # Returns
    /// A `Result` object containing a unit type
    /// 
    /// # Example
    ///
    /// ```
    /// use chaintools as chain;
    ///
    /// let data = chain::Indexer::index("/path/to/file")?;
    /// // the index will be saved to /path/to/file.ix
    /// ```
    pub fn index<U>(file: U) -> Result<()>
    where
        U: AsRef<Path> + Debug + Display,
    {
        // open file as binary
        let data: Vec<u8> = Self::open(&file)?;
        let mut data: &[u8] = &data[..];

        // create the necessary variables for position tracking
        let out_file: String = format!("{}.ix", &file);
        let mut out_handle = File::create(&out_file[..])?;
        // let mut results: Vec<(u64, usize, usize)> = Vec::new();
        // let mut results: FxHashMap<u64, (usize, usize)> = FxHashMap::default();
        let mut start_byte: usize = 0;
        let mut end_byte: usize;
        let mut header_end: usize = 0;
        let mut offset: usize = 0;
        // let mut curr_byte: u64 = 0;
        let mut chain: u64 = 0;
        let mut chain_encountered: bool = false;
        // here goes the indexing code

        // iterate over file
        loop {
            let Some(chain_start) = memchr(b'c', &data) else {
                // record last byte to be the end byte for the last chain
                end_byte = start_byte + header_end + 1+ data.len();
                // results.insert(chain, (start_byte, end_byte));
                writeln!(out_handle, "{}\t{}\t{}", chain, start_byte, end_byte)?; // TODO: Better error handling?
                // our job is done here
                break
            };
            // once the `chain` keyword is encountered, record the start byte
            // if another chain record has been already encountered, save the start and end bytes for it
            if chain_encountered {
                end_byte = chain_start + offset;//start_byte + chain_start + header_end;//j  + header_end;
                // results.insert(chain, (start_byte, end_byte));
                writeln!(out_handle, "{}\t{}\t{}", chain, start_byte, end_byte)?; // TODO: Better error handling?
                // results.push((chain, start_byte, end_byte));
            }
            data = &data[chain_start..];
            header_end = memchr(b'\n', &data).with_context(|| {
                
                format!(
                    "Failed to find separator in: {:?}. Bad formatted line!",
                    String::from_utf8_lossy(data)
                )
            })?;
            let header = &data[..header_end];
            let id_field: usize = memrchr(b' ', header).with_context(|| {
                format!(
                    "Improperly formatted header line: {:?}!",
                    String::from_utf8_lossy(header)
                )
            }
            )?;
            chain = std::str::from_utf8(&header[id_field+1..])
                .unwrap()
                .to_string()
                .parse::<u64>()?;
            start_byte = chain_start + offset; //+ header_end;
            chain_encountered = true;
            data = &data[header_end+1..];
            offset += chain_start + header_end + 1;
        }

        // serialize the results and exit
        // Writer::to_bin(results, out_file)
        // f = File::write(out_file);


        Ok(())
    }

    /// Reads chain file index
    /// 
    /// # Arguments
    /// 
    /// * `index_file` - A path to an index file
    /// 
    /// # Returns
    /// 
    /// A {chain_id: (first_byte, last_byte)} FxHashMap
    /// 
    pub fn read_index<U>(index_file: U) -> Result<FxHashMap<u64, (usize, usize)>>
    where
        U: AsRef<Path> + Debug + Display
    {
        let data = Self::open(index_file)?;
        let index: FxHashMap<u64, (usize, usize)> = bincode::deserialize(&data).with_context(|| 
            "Deserialization failed"
        )?;
        Ok(index)
    }
}
