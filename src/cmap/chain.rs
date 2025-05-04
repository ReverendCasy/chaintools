use anyhow::{Context, Result};
use fxhash::FxHashMap;
use memchr::memchr;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
// use std::mem;
use std::str::from_utf8;

use crate::cmap::align::AlignmentRecord;

/// YM - Current issues
/// 
/// 1) Are u64 values really needed for size/coordinate specifiers or is it a bit of an overkill?

/// A discrete representation of a genomic chain.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Chain {
    pub score: u64,
    pub refs: ChainHead,
    pub query: ChainHead,
    pub alignment: Vec<AlignmentRecord>,
    pub id: u32,
}

/// [YM] A structure to represent alignment records as tuples of genomic coordinates
#[derive(Clone, Debug)]
pub enum ChainBlock {
    OneSided { id: String, start: u64, end: u64 },
    DoubleSided { id: String, r_start: u64, r_end: u64, q_start: u64, q_end: u64 }
}
 

/// [YM] An enum specifying for which assemblies block coordinates should be extracted
pub enum BlockSide {
    Ref,
    Query,
    Both
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
    /// use chaintools as chain;
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
    /// use chaintools as chain;
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
    /// use chaintools as chain;
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
    /// use chaintools as chain;
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
    /// use chaintools as chain;
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
    /// use chaintools as chain;
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
    /// use chaintools as chain;
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
    /// use chaintools as chain;
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
    /// use chaintools as chain;
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
    /// use chaintools as chain;
    ///
    /// let head = b"chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1\n";
    /// let data = chain::Chain::head(head);
    /// println!("{:?}", data);
    ///
    /// > (4900, ChainHead { chr: "chrY", size: 58368225, strand: '+', start: 25985403, end: 25985638 },
    /// ChainHead { chr: "chr5", size: 151006098, strand: '-', start: 43257292, end: 43257528 }, 1);
    /// ```
    pub fn head(header: &[u8]) -> Result<(u64, ChainHead, ChainHead, u32)> {
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
            .parse::<u64>()
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

    /// [YM]
    /// Convert the chain into a vector of chain blocks
    /// 
    /// # Arguments
    /// * `self` - A Chain object
    /// 
    /// # Returns
    /// 
    /// * Vector<ChainBlock>
    /// 
    /// # Example
    /// 
    /// ```
    /// use chaintools as chain;
    /// 
    /// let data = chain::Reader::from_file("/path/to/chainfile")?;
    /// ```
    /// 
    /// TODO: 
    /// Implement as a mpsc channel that yields the blocks instead of returning the whole vector
    pub fn to_blocks(&self, side: BlockSide, report_gaps: bool) -> Vec<ChainBlock> {
        let mut r_start: u64 = self.refs.start;
        let q_strand: bool = self.query.strand == '+';
        let mut q_start: u64 = match q_strand {
            true => self.query.start,
            false => self.query.size - self.query.start
        };
        let mut block_num: u32 = 1;
        let mut blocks: Vec<ChainBlock> = Vec::new();
        let mut q_block_start: u64;
        let mut q_block_end: u64;
        // iterate over alignment records
        for b in &self.alignment {
            // the reported data structure depends on the 'side' argument value
            let mut r_block_end: u64 = r_start + (b.size as u64);
            if q_strand {
                q_block_start = q_start;
                q_block_end = q_block_start + (b.size as u64);
            } else {
                q_block_start = q_start - (b.size as u64);
                q_block_end = q_start;
            }
            let block: ChainBlock = match side {
                // report only the reference coordinates
                BlockSide::Ref => {
                    ChainBlock::OneSided{id: block_num.to_string(), start: r_start, end: r_block_end}
                },
                BlockSide::Query => {
                    ChainBlock::OneSided{id: block_num.to_string(), start: q_block_start, end: q_block_end}
                },
                BlockSide::Both => {
                    ChainBlock::DoubleSided { id: block_num.to_string(), r_start: r_start, r_end: r_block_end, q_start: q_block_start, q_end: q_block_end }
                }
            };
            blocks.push(block);
            // r_start += (b.size + b.dt) as u64;
            // q_start = if q_strand {q_start + b.size as u64} else {q_start - b.size as u64};
            // r_start += b.dt as u64;
            // q_start = if q_strand {q_start + (b.dq + b.size) as u64} else {q_start - (b.dq + b.size) as u64};
            r_start += b.size as u64;
            q_start = if q_strand {q_start + b.size as u64} else {q_start - b.size as u64};
            // if chain gap blocks were requested and the first coding block has been passed,
            // add a gap object
            if report_gaps && !(b.dq == 0 && b.dq == 0) {
                let gap_name: String = format!("{}_{}", block_num, block_num + 1);
                r_block_end = r_start + b.dt as u64;
                if q_strand {
                    q_block_start = q_start;
                    q_block_end = q_block_start + (b.dq as u64);
                } else {
                    q_block_start = q_start - (b.dq as u64);
                    q_block_end = q_start;
                }
                let block: ChainBlock = match side {
                    BlockSide::Ref => {
                        ChainBlock::OneSided { id: gap_name, start: r_start, end: r_block_end }
                    },
                    BlockSide::Query => {
                        ChainBlock::OneSided { id: gap_name, start: q_block_start, end: q_block_end }
                    },
                    BlockSide::Both => {
                        ChainBlock::DoubleSided { 
                            id: gap_name, r_start: r_start, r_end: r_block_end, q_start: q_block_start, q_end: q_block_end 
                        }
                    }
                };
                blocks.push(block);
            }
            r_start += b.dt as u64;
            q_start = if q_strand {q_start + b.dq as u64} else {q_start - b.dq as u64};
            block_num += 1;
        }
        blocks
    }

    /// [YM] + NOT FINISHED
    /// Maps coordinates from reference to query
    /// 
    /// # Arguments
    /// 
    /// `intervals` - A collection of objects having "start" and "end" coordinates; using tuples for now
    /// TODO: Define valid types 
    /// 
    /// `abs_threshold` - An absolute value by which an unaligned coordinated can be extrapolated
    /// 
    /// `rel_threshold` - A multiplier of an interval's length specifying the relative threshold of extrapolation
    pub fn map_through(
        &self, 
        intervals: &mut Vec<(&str, u64, u64, &str)>,
        abs_threshold: u64,
        rel_threshold: f64
) -> FxHashMap<&str, (u64, u64)> {
        let output: FxHashMap<&str, (u64, u64)> = FxHashMap::default();
        intervals.sort_by(
        |a, b| if a.1 == b.1 {
            a.2.cmp(&b.2)
        } else {
            a.1.cmp(&b.1)
        }
        );
        let mut min_start: u64 = intervals[0].1;
        let mut max_end: u64 = intervals[intervals.len()].2; // will this panic??
        // create a smart iteration index; iteration will always start from this interval
        let mut curr: usize = 0;
        // record the current interval's end coordinate; this will ensure that the iterator will never
        // skip the nested intervals
        let mut curr_end: u64 = intervals[0].2;
        // create a hash map of relative length threshold; for long interval lists 
        // retrieving those from an array might be faster than calculating them every time anew
        let mut rel_sizes: FxHashMap<&str, u64> = FxHashMap::default();

        // define whether alignment is codirected between reference in query
        // for now we assume that chains always represent the positive strand in the reference sequence
        // this means, 'codirectionality' depends on the query strand alone
        let codirected: bool = &self.query.strand == &'+';

        // initialize the variables standig for block coordinates
        // (see TODO tho)
        // 
        let mut r_start: u64 = self.refs.start;
        let q_strand: bool = self.query.strand == '+';
        let mut q_start: u64 = match q_strand {
            true => self.query.start,
            false => self.query.size - self.query.start
        };
        let mut q_block_start: u64;
        let mut q_block_end: u64;

        // finally, initialize the projected coordinate variables
        let mut start_p: u64 = 0;
        let mut end_p: u64 = 0;

        // all set
        // now, iterate over alignment records
        // TODO: Implement to_blocks() as yielder to avoid code repetition
        for b in &self.alignment {
            // break if the iterator has passed beyond the last interval
            if r_start > max_end {break};
            let mut r_block_end: u64 = r_start + (b.size as u64);
            // skip the block preceding the first interval's start in the reference
            if r_block_end < min_start {continue};
            // define the query coordinates
            if q_strand {
                q_block_start = q_start;
                q_block_end = q_block_start + (b.size as u64);
            } else {
                q_block_start = q_start - (b.size as u64);
                q_block_end = q_start;
            }

            // now, we have a chain block with defined boundaries in both reference and query
            // iteratve over the intervals, chek whether any of their coordinates can be projected 
            // through this block
            for (i, inter) in intervals[curr..].iter().enumerate() {
                // check whether the start coorindate is within the interval
                if (r_start <= inter.1) && (inter.1 <= r_block_end) {
                    let offset: u64 = inter.1 - r_start;
                    if codirected{
                        start_p = q_block_start + offset;
                        // assign to a storage variable
                    } else {
                        end_p = q_block_end - offset;
                        // assign to a storage variable
                    }
                }
                // then, check the end coordinate
                if (r_start <= inter.2) && (inter.2 <= r_block_end) {
                    let offset: u64 = r_block_end - inter.2;
                    if codirected {
                        end_p = q_block_end - offset;
                        // assign to a storage variable
                    } else {
                        start_p = q_block_start + offset;
                        // assign to a storage variable
                    }
                }
                // if an interval has bots its coordinates properly mapped,
                // update the current transcript pointer
            }

            // update the block coordinates
            r_start += b.size as u64;
            q_start = if q_strand {q_start + b.size as u64} else {q_start - b.size as u64};
            // now, process the alignment gap
            // somewhat less trivial than mapping through an aligned block 
            r_block_end = r_start + b.dt as u64;
            if q_strand {
                q_block_start = q_start;
                q_block_end = q_block_start + (b.dq as u64);
            } else {
                q_block_start = q_start - (b.dq as u64);
                q_block_end = q_start;
            }
            // now, iterate through the remaining intervals again
            for (i, inter) in intervals[curr..].iter().enumerate() {
                // start coordinate is within the alignment gap
                if (r_start <= inter.1) && (inter.1 <= r_block_end) {
                    // get the alignment offset
                    let offset: u64 = inter.1 - r_start;
                    // get the relative threshold size
                    let rel_thresh: &u64 = rel_sizes
                        .entry(inter.3)
                        .or_insert((inter.2 - inter.1) * rel_threshold as u64);

                    // check if the offset is within the stated extrapolation limits 
                    if offset > abs_threshold && offset > *rel_thresh {
                        // coordinate is too far to be extrapolated; crop to the chain block's start
                        if codirected {
                            start_p = q_block_start;
                            // assign to a storage variable
                        } else {
                            end_p = q_block_end;
                            // assign to a storage variable
                        }
                    } else {
                        // extrapolated sequence's length does not exceed the stated thresholds
                        if codirected {
                            start_p = q_block_start + offset;
                            // assign to a storage variable
                        } else {
                            end_p = q_block_end - offset;
                            // assign to a storage variable
                        }
                    }
                }

                // and the same for end coordinate
                if (r_start <= inter.2) && (inter.2 <= r_block_end) {
                    // get the alignment offset
                    let offset: u64 = r_block_end - inter.2;
                    // get the relative threshold size
                    let rel_thresh: &u64 = rel_sizes
                        .entry(inter.3)
                        .or_insert((inter.2 - inter.1) * rel_threshold as u64);
                    
                    // check if the offset is within the stated extrapolation limits 
                    if offset > abs_threshold && offset > *rel_thresh {
                        // coordinate is too far to be extrapolated; crop to the chain block's start
                        if codirected {
                            end_p = q_block_end;
                            // assign to a storage variable
                        } else {
                            start_p = q_block_start;
                            // assign to a storage variable
                        }
                    } else {
                        // extrapolated sequence's length does not exceed the stated thresholds
                        if codirected {
                            end_p = q_block_end - offset;
                            // assign to a storage variable
                        } else {
                            start_p = q_block_start + offset;
                            // assign to a storage variable
                        }
                    }
                    // end coordinate has been mapped accordingly
                    // current interval pointer can be updated
                }
            }

            // proceed to the next line
            r_start += b.dt as u64;
            q_start = if q_strand {q_start + b.dq as u64} else {q_start - b.dq as u64};
        }
        output
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
    /// use chaintools as chain;
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
    /// use chaintools as chain;
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
    /// use chaintools as chain;
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
