use anyhow::{Context, Result};
use memchr::memchr;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::str::from_utf8;

use crate::align::AlignmentRecord;

/// A discrete representation of a genomic chain.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Chain {
    pub score: i64,
    pub refs: ChainHead,
    pub query: ChainHead,
    pub alignment: Vec<AlignmentRecord>,
    pub id: u32,
}

impl Chain {
    /// Create a new chain object from a chain block (header, alignment).
    ///
    /// # Arguments
    /// * `head` - A byte array containing the header of the chain block
    /// * `block` - A byte array containing the alignment of the chain block
    ///
    /// # Returns
    /// * Option<(u32, Self)>
    /// * None if the block is empty
    ///
    /// # Example
    ///
    /// ```
    /// use chainder as chain;
    ///
    /// let head = b"chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1\n";
    /// let block = b"9\t1\t0\n";
    /// let data = chain::Chain::from(head, block);
    /// println!("{:?}", data);
    ///
    /// > Some((1, Chain { score: 4900, refs: ChainHead { chr: "chrY", size: 58368225, strand: '+', start: 25985403, end: 25985638 },
    /// query: ChainHead { chr: "chr5", size: 151006098, strand: '-', start: 43257292, end: 43257528 },
    /// alignment: [AlignmentRecord { size: 9, dt: 1, dq: 0 }], id: 1 }));
    /// ```
    pub fn from(head: &[u8], block: &[u8]) -> Result<(u32, Self)> {
        let (score, refs, query, id) = Self::head(head)?;
        let alignment = AlignmentRecord::parse(block);

        Ok((
            id,
            Self {
                score,
                refs,
                query,
                alignment,
                id,
            },
        ))
    }

    /// Get the chain object as a string.
    ///
    /// # Arguments
    /// * `self` - A chain object
    ///
    /// # Returns
    /// * String
    ///
    /// # Example
    ///
    /// ```
    /// use chainder as chain;
    ///
    /// let data = chain::Reader::from_file("/path/to/chainfile")?;
    /// let line = data.get(&12).unwrap().to_string();
    /// println!("{:?}", line);
    ///
    /// > "chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 12\n9\t1\t0\n\n";
    /// ```
    pub fn to_string(&self) -> String {
        format!("{}\n{}", self.header(), self.alignment())
    }

    /// Get the header of a chain object as a string.
    ///
    /// # Arguments
    /// * `self` - A chain object
    ///
    /// # Returns
    /// * String
    ///
    /// # Example
    ///
    /// ```
    /// use chainder as chain;
    ///
    /// let data = chain::Reader::from_file("/path/to/chainfile")?;
    /// let line = data.get(&12).unwrap().header();
    /// println!("{:?}", line);
    ///
    /// > "chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 12";
    /// ```
    pub fn header(&self) -> String {
        format!(
            "chain {} {} {} {}",
            self.score,
            self.refs.to_string(),
            self.query.to_string(),
            self.id
        )
    }

    /// Get the header of a chain object as a vector of strings.
    /// The vector contains the score, refs, query, and id values.
    ///
    /// # Arguments
    /// * `self` - A chain object
    ///
    /// # Returns
    /// * Vec<String>
    ///
    /// # Example
    /// ```
    /// use chainder as chain;
    ///
    /// let data = chain::Reader::from_file("/path/to/chainfile")?;
    /// let line = data.get(&12).unwrap().header_vec();
    /// println!("{:?}", line);
    ///
    /// > ["4900", "chrY", "58368225", "+", "25985403", "25985638", "chr5", "151006098", "-", "43257292", "43257528", "12"];
    /// ```
    pub fn header_vec(&self) -> Vec<String> {
        let mut hv = vec![String::from("chain"), self.score.to_string()];
        hv.extend(self.refs.to_vec());
        hv.extend(self.query.to_vec());
        hv.push(self.id.to_string());
        hv
    }

    /// Get the alignment of a chain object as a string.
    ///
    /// # Arguments
    /// * `self` - A chain object
    ///
    /// # Returns
    /// * String
    ///
    /// # Example
    ///
    /// ```
    /// use chainder as chain;
    ///
    /// let data = chain::Reader::from_file("/path/to/chainfile")?;
    /// let line = data.get(&12).unwrap().alignment();
    /// println!("{:?}", line);
    ///
    /// > "9\t1\t0\n\n";
    /// ```
    pub fn alignment(&self) -> String {
        let len = self.alignment.len();
        self.alignment
            .iter()
            .enumerate()
            .map(|(i, a)| {
                if i == len - 1 {
                    format!("{}\n\n", a.size.to_string())
                } else {
                    format!("{}\n", a.to_string())
                }
            })
            .collect::<String>()
    }

    /// Get the alignment of a chain object as a vector of vectors.
    /// The inner vector contains the size, dt, and dq values.
    ///
    /// # Arguments
    /// * `self` - A chain object
    ///
    /// # Returns
    /// * Vec<Vec<u32>>
    ///
    /// # Example
    ///
    /// ```
    /// use chainder as chain;
    ///
    /// let data = chain::Reader::from_file("/path/to/chainfile")?;
    /// let line = data.get(&12).unwrap().alignment_vec();
    /// println!("{:?}", line);
    ///
    /// > [[9, 1, 0]];
    /// ```
    pub fn alignment_vec(&self) -> Vec<Vec<u32>> {
        self.alignment
            .iter()
            .map(|a| vec![a.size, a.dt, a.dq])
            .collect()
    }

    /// Get the alignment of a chain object as a vector of vectors.
    /// The inner vector contains the size, dt, and dq values.
    ///
    /// # Arguments
    /// * `self` - A chain object
    ///
    /// # Returns
    /// * Vec<Vec<String>>
    ///
    /// # Example
    ///
    /// ```
    /// use chainder as chain;
    ///
    /// let data = chain::Reader::from_file("/path/to/chainfile")?;
    /// let line = data.get(&12).unwrap().alignment_vec_str();
    /// println!("{:?}", line);
    ///
    /// > [["9", "1", "0"]];
    /// ```
    pub fn alignment_vec_str(&self) -> Vec<Vec<String>> {
        self.alignment
            .iter()
            .map(|a| vec![a.size.to_string(), a.dt.to_string(), a.dq.to_string()])
            .collect()
    }

    /// Get the chain object as a byte array.
    ///
    /// # Arguments
    /// * `self` - A chain object
    ///
    /// # Returns
    /// * Vec<u8>
    ///
    /// # Example
    ///
    /// ```
    /// use chainder as chain;
    ///
    /// let data = chain::Reader::from_file("/path/to/chainfile")?;
    /// let line = data.get(&12).unwrap().to_bytes();
    /// println!("{:?}", line);
    ///
    /// > [99, 104, 97, 105, 110, 32, 52, 57, 48, 48, 32, 99, 104, 114, 89, 32, 53, 56, 51,
    /// 54, 56, 50, 50, 53, 32, 43, 32, 50, 53, 57, 56, 53, 52, 48, 51, 32, 50, 53, 57, 56,
    /// 53, 54, 51, 56, 32, 99, 104, 114, 53, 32, 49, 53, 49, 48, 48, 54, 48, 57, 56, 32,
    /// 45, 32, 52, 51, 49, 54, 52, 48, 51, 50, 57, 50, 32, 52, 51, 49, 54, 52, 48, 51, 53,
    /// 50, 56, 32, 49, 50, 10, 57, 9, 49, 9, 48, 10];
    /// ```
    pub fn to_bytes(&self) -> Vec<u8> {
        self.to_string().into_bytes()
    }

    /// Get the chain object as a vector of vectors.
    /// The inner vector contains the header and alignment values.
    /// The header vector contains the score, refs, query, and id values.
    /// The alignment vector contains the size, dt, and dq values.
    ///
    /// # Arguments
    /// * `self` - A chain object
    ///
    /// # Returns
    /// * Vec<Vec<String>>
    ///
    /// # Example
    /// ```
    /// use chainder as chain;
    ///
    /// let data = chain::Reader::from_file("/path/to/chainfile")?;
    /// let line = data.get(&12).unwrap().to_vec();
    /// println!("{:?}", line);
    ///
    /// > [["chain", "4900", "chrY", "58368225", "+", "25985403", "25985638",
    /// "chr5", "151006098", "-", "43257292", "43257528", "12"], ["9", "1", "0"]];
    /// ```
    pub fn to_vec(&self) -> Vec<Vec<String>> {
        let mut v = self.alignment_vec_str();
        v.insert(0, self.header_vec());
        v
    }

    /// Process a chain header into a chainer compatible pre-processing format.
    ///
    /// # Arguments
    /// * `header` - A byte array containing the header of the chain block
    ///
    /// # Returns
    /// * (i64, ChainHead, ChainHead, u32)
    ///
    /// # Example
    ///
    /// ```
    /// use chainder as chain;
    ///
    /// let head = b"chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1\n";
    /// let data = chain::Chain::head(head);
    /// println!("{:?}", data);
    ///
    /// > (4900, ChainHead { chr: "chrY", size: 58368225, strand: '+', start: 25985403, end: 25985638 },
    /// ChainHead { chr: "chr5", size: 151006098, strand: '-', start: 43257292, end: 43257528 }, 1);
    /// ```
    pub fn head(header: &[u8]) -> Result<(i64, ChainHead, ChainHead, u32)> {
        let mut acc = vec![];
        let mut header = &header[..];
        loop {
            let Some(sep) = memchr(b' ', header) else {
                acc.push(header);
                break;
            };
            acc.push(&header[..sep]);
            header = &header[sep + 1..];
        }

        let refs = ChainHead::from(&acc[2..7])?;
        let query = ChainHead::from(&acc[7..12])?;

        let score = from_utf8(&acc[1])
            .unwrap()
            .parse::<i64>()
            .with_context(|| {
                format!(
                    "Failed to parse score in: {:?}. Bad formatted line!",
                    String::from_utf8_lossy(header)
                )
            })?;

        let id = from_utf8(&acc.last().unwrap())
            .unwrap()
            .parse::<u32>()
            .with_context(|| {
                format!(
                    "Failed to parse id in: {:?}. Bad formatted line!",
                    String::from_utf8_lossy(header)
                )
            })?;

        Ok((score, refs, query, id))
    }
}

/// A ref/query chain head object.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChainHead {
    pub chr: String,
    pub size: u64,
    pub strand: char,
    pub start: u64,
    pub end: u64,
}

impl ChainHead {
    /// Create a new chain head object from a byte array representing a ref or query header.
    ///
    /// # Arguments
    /// * `header` - A byte array containing the header of the chain block
    ///
    /// # Returns
    /// * Result<ChainHead>
    /// * An error if the header is not well formatted
    ///
    /// # Example
    ///
    /// ```
    /// use chainder as chain;
    ///
    /// let head = b"chrY 58368225 + 25985403 25985638";
    /// let data = chain::ChainHead::from(head);
    /// println!("{:?}", data);
    ///
    /// > Ok(ChainHead { chr: "chrY", size: 58368225, strand: '+', start: 25985403, end: 25985638 });
    pub fn from(header: &[&[u8]]) -> Result<Self> {
        Ok(Self {
            chr: from_utf8(header[0])
                .context("Failed to decode chr data as UTF-8")?
                .to_string(),

            size: from_utf8(header[1])?.parse::<u64>().with_context(|| {
                format!(
                    "Failed to parse size in: {:?}. Bad formatted line!",
                    String::from_utf8_lossy(header[1])
                )
            })?,

            strand: from_utf8(header[2])?.chars().next().with_context(|| {
                format!(
                    "Failed to parse strand in: {:?}. Bad formatted line!",
                    String::from_utf8_lossy(header[2])
                )
            })?,

            start: from_utf8(header[3])?.parse::<u64>().with_context(|| {
                format!(
                    "Failed to parse start in: {:?}. Bad formatted line!",
                    String::from_utf8_lossy(header[3])
                )
            })?,

            end: from_utf8(header[4])?.parse::<u64>().with_context(|| {
                format!(
                    "Failed to parse end in: {:?}. Bad formatted line!",
                    String::from_utf8_lossy(header[4])
                )
            })?,
        })
    }

    /// Get the chain head object as a string.
    ///
    /// # Arguments
    /// * `self` - A chain head object
    ///
    /// # Returns
    /// * String
    ///
    /// # Example
    ///
    /// ```
    /// use chainder as chain;
    ///
    /// let head = b"chrY 58368225 + 25985403 25985638";
    /// let data = chain::ChainHead::from(head).unwrap().to_string();
    /// println!("{:?}", data);
    ///
    /// > "chrY 58368225 + 25985403 25985638";
    pub fn to_string(&self) -> String {
        format!(
            "{} {} {} {} {}",
            self.chr, self.size, self.strand, self.start, self.end
        )
    }

    /// Get the chain head object as a vector of strings.
    ///
    /// # Arguments
    /// * `self` - A chain head object
    ///
    /// # Returns
    /// * Vec<String>
    ///
    /// # Example
    /// ```
    /// use chainder as chain;
    ///
    /// let head = b"chrY 58368225 + 25985403 25985638";
    /// let data = chain::ChainHead::from(head).unwrap().to_vec();
    /// println!("{:?}", data);
    ///
    /// > ["chrY", "58368225", "+", "25985403", "25985638"];
    pub fn to_vec(&self) -> Vec<String> {
        vec![
            self.chr.to_string(),
            self.size.to_string(),
            self.strand.to_string(),
            self.start.to_string(),
            self.end.to_string(),
        ]
    }
}
