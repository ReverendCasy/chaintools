use anyhow::{Context, Result};
use cubiculum::merge::merge::intersection;
use cubiculum::structs::structs::{BedEntry, Coordinates, Interval, Named};
use fxhash::FxHashMap;
use std::cmp::{max, min, Ord};
use std::fmt::Debug;
use yield_return::LocalIter;

use crate::cmap::chain::{BlockSide, ChainBlock, DoubleSidedBlock, OneSidedBlock};

impl crate::cmap::chain::Chain {
    /// [YM] Given a vector of cubiculum Interval-like objects, returns a vector
    /// of items overlapping the chain's span
    /// 
    /// # Arguments
    /// `intervals`: A vector of interval objects 
    /// `is_ref`: boolean value indicating whether reference coordinates should be used for the chain;
    /// using query coordinates otherwise
    /// 
    /// # Returns
    /// A Vector of Interval objects whose coordinates overlap the chain span 
    /// 
    pub fn intersect_to_vector<T>(&self, intervals: &Vec<T>, to_ref: bool) -> Vec<T>
    where 
            T: Coordinates + Named + Clone + Debug
    {
        let mut output: Vec<T> = Vec::<T>::new();
        let start: &u64 = if to_ref {&self.refs.start} else {
            if self.query.strand == '+' {&self.query.start} else {&(self.query.size - self.query.end)}
        };
        let end: &u64 = if to_ref {&self.refs.end} else {
            if self.query.strand == '+' {&self.query.end} else {&(self.query.size - self.query.start)}
        };
        for i in intervals {
            let inter_start = match i.start() {
                Some(x) => {x},
                None => continue
            };
            if inter_start >= end {break}
            let inter_end = match i.end() {
                Some(x) => {x},
                None => continue
            };
            if inter_end <= start {continue};
            output.push(i.clone());
        }
        output

    }

    /// [YM] Am implementation of intersect_to_vector() designed for BED8+ BedEntry objects
    /// 
    /// # Arguments
    /// `intervals`: A vector of BedEntry objects. All items are expected to be in BED8+ format 
    /// and have defined thickStart and thickEnd coordinates
    /// `is_ref`: boolean value indicating whether reference coordinates should be used for the chain;
    /// using query coordinates otherwise
    /// 
    /// # Returns
    /// A Vector of BedEntry objects whose coding ('thick') intervals overlap the chain span 
    /// 
    pub fn intersect_to_cds_vector(&self, intervals: &Vec<BedEntry>, to_ref: bool) -> Vec<BedEntry>
    {
        let mut output: Vec<BedEntry> = Vec::<BedEntry>::new();
        let start: u64 = if to_ref {self.refs.start} else {
            if self.query.strand == '+' {self.query.start} else {self.query.size - self.query.end}
        };
        let end: u64 = if to_ref {self.refs.end} else {
            if self.query.strand == '+' {self.query.end} else {self.query.size - self.query.start}
        };
        for i in intervals {
            if i.format() < 8 {continue}
            let inter_start = match i.thick_start() {
                Some(x) => {x},
                None => continue
            };
            if inter_start >= end {
                if self.id == 546368 {
                    println!(
                        "Breaking at {}; chain coords: {}-{}, inter coords: {}-{}", 
                        i.name().unwrap(), start, end, inter_start, i.thick_end().unwrap()
                    )
                }
                break
            }
            let inter_end = match i.thick_end() {
                Some(x) => {x},
                None => continue
            };
            if inter_end <= start {
                if self.id == 546368 {
                    println!(
                        "Continuing after {}; chain coords: {}-{}, inter coords: {}-{}", 
                        i.name().unwrap(), start, end, inter_start, inter_end
                    )
                }
                continue
            };
            output.push(i.clone());
        }
        output

    }

    /// [YM - A first attempt at block yielder]
    /// A 'yielding' version of to_blocks() from cmap::chain; 
    /// returns blocks as a generator-like object of Futures;
    /// not intended to be public, rather aimed at optimizing other functions in the module
    /// 
    /// NOTE: Currently solutions based on this method require more time 
    /// than 'real-time' parsers; this function will be likely modified in future 
    /// to create a real concurrent data channel
    /// 
    /// # Arguments
    /// 
    /// * `side`: A BlockSide enum field, specifying coordinates in which genome 
    /// (reference, query, or both) should be reported
    /// 
    /// * `report_gaps`: a boolean value indicating whether alignment chain gaps should be 
    /// also reported as blocks; by default, only aligned blocks are reported
    /// 
    /// #  Returns
    /// 
    /// A LocalIter<> object of ChainBlock structs
    /// 
    /// ```
    /// use chaintools as chain;
    /// use chain::cmap::chain::BlockSide::*;
    /// 
    /// let head = b"chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1\n";
    /// let block = b"9\t1\t0\n";
    /// let data = chain::Chain::from(head, block);
    /// for x in chain.yield_blocks(Ref, false) {
    ///     println!("Block={:#?}", block);
    /// }
    /// ```
    fn yield_blocks<'a>(
        &'a self, side: BlockSide, 
        report_gaps: bool,
    ) -> yield_return::LocalIter<'a, Box<dyn ChainBlock>>
    {
        let generator = LocalIter::new(|mut x| async move {
            let mut r_start: u64 = self.refs.start;
            let q_strand: bool = self.query.strand == '+';
            let mut q_start: u64 = match q_strand {
                true => self.query.start,
                false => self.query.size - self.query.start
            };
            let mut block_num: u32 = 1;
            // let mut blocks: Vec<ChainBlock> = Vec::new();
            let mut q_block_start: u64;
            let mut q_block_end: u64;
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
                let block: Box<dyn ChainBlock> = match side {
                    // report only the reference coordinates
                    BlockSide::Ref => {
                        Box::new(
                            OneSidedBlock::new(block_num.to_string(), true, r_start, r_block_end)
                        )
                    },
                    BlockSide::Query => {
                        Box::new(
                            OneSidedBlock::new(block_num.to_string(),false, q_block_start, q_block_end)
                        )
                    },
                    BlockSide::Both => {
                        Box::new(
                            DoubleSidedBlock::new(block_num.to_string(), r_start, r_block_end, q_block_start, q_block_end)
                        )
                    }
                };
                // x.ret(block).await;
                x.ret(block).await;
                // blocks.push(block);
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
                    let block: Box<dyn ChainBlock> = match side {
                        BlockSide::Ref => {
                            Box::new(
                                OneSidedBlock::new(gap_name, true, r_start, r_block_end)
                            )
                        },
                        BlockSide::Query => {
                            Box::new(
                                OneSidedBlock::new(gap_name, false, q_block_start, q_block_end)
                            )
                        },
                        BlockSide::Both => {
                            Box::new(
                                DoubleSidedBlock::new( 
                                    gap_name, r_start, r_block_end, q_block_start, q_block_end 
                                )
                            )
                        }
                    };
                    // blocks.push(block);
                    // x.ret(block).await;
                    x.ret(block).await;
                }
                r_start += b.dt as u64;
                q_start = if q_strand {q_start + b.dq as u64} else {q_start - b.dq as u64};
                block_num += 1;
            }

            });
            generator
        }

        pub fn alignment_cov<'a, T>(&self, intervals: &'a mut Vec<T>,) -> Result<FxHashMap<&'a str, u64>> 
        where 
            T: Coordinates + Named + Debug
        {
            // the same routine as above
            // first, sort the input vector
            let mut output: FxHashMap<&str, u64> = FxHashMap::default();
            intervals.sort_by(
                |a, b| if a.start().unwrap() == b.start().unwrap() {
                    a.end().unwrap().cmp(&b.end().unwrap())
                } else {
                    a.start().unwrap().cmp(&b.start().unwrap())
                }
                );
            // define the total span for the input intervals
            let min_start: u64 = *intervals[0].start().with_context(||
                {"Cannot assess coverage for intervals with undefined coordinates"}
            )?;
            let max_end: u64 = *intervals[intervals.len() - 1].end().with_context(||
                {"Cannot assess coverage for intervals with undefined coordinates"}
            )?; // will this panic??
            // create a smart iteration index; iteration will always start from this interval
            let curr: usize = 0;
            // record the current interval's end coordinate; this will ensure that the iterator will never
            // skip the nested intervals
            let curr_end: u64 = *intervals[0].end().with_context(||
                {"Cannot assess coverage for intervals with undefined coordinates"}
            )?;
    
            // create a smart iteration index; iteration will always start from this interval
            let mut curr: usize = 0;
            // record the current interval's end coordinate; this will ensure that the iterator will never
            // skip the nested intervals
            let mut curr_end: u64 = *intervals[0].end().with_context(||
                {"Cannot assess coverage for intervals with undefined coordinates"}
            )?;
    
            // now go
            for (h, b) in self.yield_blocks(BlockSide::Both, false).enumerate() {
                let b_r_start = b.r_start().unwrap();
                let b_r_end = b.r_end().unwrap();

                // continue if the first interval has not yet been reached, break if the last one has been passed
                if b_r_end < min_start {continue};
                if b_r_start > max_end {break};

                for (mut i, inter) in intervals[curr..].iter().enumerate() {
                    i += curr;
                    let inter_start: u64 = *inter.start().with_context(||
                        {format!("Interval {} has an undefined start coordinate which cannot be mapped", i)}
                    )?;
                    let inter_end: u64 = *inter.end().with_context(||
                        {format!("Interval {} has an undefined end coordinate which cannot be mapped", i)}
                    )?;
                    let name: &str = inter.name().with_context(||
                        {"Interval is not named"}
                    )?;
    
                    if !output.contains_key(&inter.name().unwrap()) {
                        output.insert(
                            inter.name().unwrap(),
                            0
                        );
                    }
    
                    // chain block is upstream to the current interval;
                    // since other are guaranteed to start at least in the same position,
                    // the current loop can be safely exited
                    if b_r_end < inter_start {
                        // potentially this is the farthest the intervals have ever reached 
                        // in terms of the  end coordinate; unless this boundary is exceeded, 
                        // the iteration start point will not be updated
                        if inter_end >= curr_end {
                            // curr = i;
                            curr_end = inter_end;
                        }
                        break
                    }
    
                    // chain block is downstream to the current interval;
                    // nothing to do here, proceed to the next interval;
                    if b_r_start > inter_end {
                        // if this interval is not a boundary of the current overlap group,
                        // current transcript pointer can be safely updated;
                        // the next iteration will start downstream to this interval or a nested interval group
                        if inter_end < curr_end {
                            curr += 1;
                        }
                        continue
                    };
    
                    // println!("Assessing intersection for {}", name);
    
                    match intersection(inter_start, inter_end, b_r_start, b_r_end) {
                        Some(x) => {
                            *output.get_mut(name).unwrap() += x;
                        },
                        None => {}
                    }
                }
            }
            Ok(output)
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
    /// 
    /// # Returns
    /// 
    /// Result<&str, Interval> where each interval contains projected coordinates for each input interval 
    pub fn map_through<'a, T>(
        &'a self, 
        // intervals: &mut Vec<(&str, u64, u64, &str)>,
        intervals: &'a mut Vec<T>,
        abs_threshold: u64,
        rel_threshold: f64
    ) -> Result<FxHashMap<&'a str, Interval>> //Result<FxHashMap<&str, (u64, u64)>> 
    where 
        T: Coordinates + Named + Debug
    {
        // let output: FxHashMap<&str, (u64, u64)> = FxHashMap::default();
        let mut output: FxHashMap<&str, Interval> = FxHashMap::default();

        intervals.sort_by(
        |a, b| if a.start().unwrap() == b.start().unwrap() {
            a.end().unwrap().cmp(&b.end().unwrap())
        } else {
            a.start().unwrap().cmp(&b.start().unwrap())
        }
        );
        // define the total span of input intervals:
        // blocks before `min_start` will be ignored; 
        // once `max_end` is passed, iteration over chain stop 
        let mut min_start: u64 = *intervals[0].start().with_context(||
            {"Cannot map intervals with undefined coordinates"}
        )?;
        let max_end: u64 = *intervals[intervals.len() - 1].end().with_context(||
            {"Cannot map intervals with undefined coordinates"}
        )?;
        // create a smart iteration index; iteration will always start from this interval
        let mut curr: usize = 0;
        // record the current interval's end coordinate; this will ensure that the iterator will never
        // skip the nested intervals
        let mut curr_end: u64 = *intervals[0].end().with_context(||
            {"Cannot map intervals with undefined coordinates"}
        )?;

        // create a hash map of relative length threshold; for long interval lists 
        // retrieving those from an array might be faster than calculating them every time anew
        let mut rel_sizes: FxHashMap<&str, u64> = FxHashMap::default();

        // define whether alignment is codirected between reference in query
        // for now we assume that chains always represent the positive strand in the reference sequence
        // this means, 'codirectionality' depends on the query strand alone
        let codirected: bool = &self.query.strand == &'+';

        // initialize the variables standing for block coordinates
        // (see TODO tho)
        // 
        let mut r_start: u64 = self.refs.start;
        let r_end: u64 = self.refs.end;
        let q_strand: bool = self.query.strand == '+';
        let mut q_start: u64 = match q_strand {
            true => self.query.start,
            false => self.query.size - self.query.start
        };

        // finally, initialize the projected coordinate variables
        let mut start_p: u64;
        let mut end_p: u64;

        // all set
        // now, iterate over alignment records
        for (h, b) in self.yield_blocks(BlockSide::Both, true).enumerate() {
            let b_r_start = b.r_start().unwrap();
            let b_r_end = b.r_end().unwrap();
            let b_q_start = b.q_start().unwrap();
            let b_q_end = b.q_end().unwrap();
            let is_gap: bool = b.is_gap();
            // break if the iterator has passed beyond the last interval
            if b_r_start > max_end {break};
            // skip the block preceding the first interval's start in the reference
            if b_r_end < min_start {
                continue
            };

            // check if this is the last block
            let is_last_block: bool = (b_r_start == b_r_end) && (b_q_start == b_q_end);

            // now, we have a chain block with defined boundaries in both reference and query;
            // iterate over the intervals, check whether any of their coordinates can be projected 
            // through this block
            for (mut i, inter) in intervals[curr..].iter().enumerate() {
                i += curr;
                let inter_start: u64 = *inter.start().with_context(||
                    {format!("Interval {} has an undefined start coordinate which cannot be mapped", i)}
                )?;
                let inter_end: u64 = *inter.end().with_context(||
                    {format!("Interval {} has an undefined end coordinate which cannot be mapped", i)}
                )?;

                // add a results block to the the output hash map
                if !output.contains_key(&inter.name().unwrap()) {
                    output.insert(
                        inter.name().unwrap(),
                        Interval::new()
                    );
                    output.
                        entry(&inter.name().unwrap())
                        .and_modify(
                            |x| {
                                x.update_name(inter.name().unwrap().to_string()); // TODO: Will borrow the value!
                                x.update_chrom(self.query.chr.clone()); // TODO: Bad choice altogether
                            }
                        );
                }

                // chain block is upstream to the current interval;
                // since other are guaranteed to start at least in the same position,
                // the current loop can be safely exited
                if b_r_end < inter_start {
                    // potentially this is the farthest the intervals have ever reached 
                    // in terms of the  end coordinate; unless this boundary is exceeded, 
                    // the iteration start point will not be updated
                    if inter_end >= curr_end {
                        // curr = i;
                        curr_end = inter_end;
                    }
                    break
                }

                // chain block is downstream to the current interval;
                // nothing to do here, proceed to the next interval;
                if b_r_start > inter_end {
                    // if this interval is not a boundary of the current overlap group,
                    // current transcript pointer can be safely updated;
                    // the next iteration will start downstream to this interval or a nested interval group
                    if inter_end < curr_end {
                        curr += 1;
                    }
                    continue
                };

                // at this point, it is ascertained that at least on coordinate of the block
                // lies within the the chain block, which makes it potentially mappable;
                // the exact behavior, however, varies depending on whether mapping is performed 
                // through an aligned chain block or am unaligned chain gap

                if !is_gap {
                    for (mut i, inter) in intervals[curr..].iter().enumerate() {
                        i += curr;
                        // check whether the start coordinate is within the block
                        if (b_r_start <= inter_start) && (inter_start <= b_r_end) {
                            //  start coordinate can be mapped
                            let offset: u64 = inter_start - b_r_start;
                            if codirected{
                                start_p = b_q_start + offset;
                                // assign to a storage variable
                                output
                                    .entry(&inter.name().unwrap())
                                    .and_modify(
                                        |x| {
                                            x.update_start(start_p)
                                        }
                                    );
                            } else {
                                end_p = b_q_end - offset;
                                // assign to a storage variable
                                output
                                    .entry(&inter.name().unwrap())
                                    .and_modify(
                                        |x| {
                                            x.update_end(end_p)
                                        }
                                    );
                            }
                            // a special case for the last block; if interval end lies outside of the chain,
                            // try extrapolating the coordinate unless it is too far from the chain 
                            if is_last_block && inter_end >  r_end {
                                // get the alignment offset
                                let offset: u64 = inter_end - b_r_end;
                                // get the relative threshold size
                                let rel_thresh: &u64 = rel_sizes
                                    .entry(
                                        inter.name().unwrap_or("a") // TODO: Find a way to create long-lived string literal IDs or update name() in cubiculum
                                    )
                                    .or_insert((inter.length().unwrap() as f64 * rel_threshold) as u64);
                                
                                // check if the offset is within the stated extrapolation limits 
                                if offset > abs_threshold && offset > *rel_thresh {
                                    // coordinate is too far to be extrapolated; crop to the chain block's start
                                    if codirected {
                                        end_p = b_q_end;
                                        // assign to a storage variable
                                        output
                                            .entry(&inter.name().unwrap())
                                            .and_modify(
                                                |x| {
                                                    x.update_end(end_p)
                                                }
                                            );
                                    } else {
                                        start_p = b_q_start;
                                        // assign to a storage variable
                                        output
                                            .entry(&inter.name().unwrap())
                                            .and_modify(
                                                |x| {
                                                    x.update_start(start_p)
                                                }
                                            );
                                    }
                                } else {
                                    // extrapolated sequence's length does not exceed the stated thresholds
                                    if codirected {
                                        end_p = b_q_end + offset;
                                        // assign to a storage variable
                                        output
                                            .entry(&inter.name().unwrap())
                                            .and_modify(
                                                |x| {
                                                    x.update_end(end_p)
                                                }
                                            );
                                    } else {
                                        start_p = b_q_start - offset;
                                        // assign to a storage variable
                                        output
                                            .entry(&inter.name().unwrap())
                                            .and_modify(
                                                |x| {
                                                    x.update_start(start_p)
                                                }
                                            );
                                    }
                                }
                            }
                        }
                        // then, check the end coordinate
                        if (b_r_start <= inter_end) && (inter_end <= b_r_end) {
                            let offset: u64 = b_r_end - inter_end;
                            if codirected {
                                end_p = b_q_end - offset;
                                // assign to a storage variable
                                output
                                    .entry(&inter.name().unwrap())
                                    .and_modify(
                                        |x| {
                                            x.update_end(end_p)
                                        }
                                    );
                            } else {
                                start_p = b_q_start + offset;
                                // assign to a storage variable
                                output
                                    .entry(&inter.name().unwrap())
                                    .and_modify(
                                        |x| {
                                            x.update_start(start_p)
                                        }
                                    );
                            }
                        }
        
                        // a special case for the first block which extends beyond the chain start
                        if h == 0 && inter_start < r_start {
                            let offset: u64 = r_start - inter_start;
                            // get the relative threshold size
                            let rel_thresh: &u64 = rel_sizes
                                .entry(
                                    inter.name().unwrap_or("a") // TODO: Find a way to create long-lived string literal IDs or update name() in cubiculum
                                )
                                .or_insert((inter.length().unwrap() as f64 * rel_threshold) as u64);
        
                            // check if the offset is within the stated extrapolation limits 
                            if offset > abs_threshold && offset > *rel_thresh {
                                // coordinate is too far to be extrapolated; crop to the chain block's start
                                if codirected {
                                    start_p = b_q_start;
                                    // assign to a storage variable
                                    output
                                        .entry(&inter.name().unwrap())
                                        .and_modify(
                                            |x| {
                                                x.update_start(start_p)
                                            }
                                        );
                                } else {
                                    end_p = b_q_end;
                                    // assign to a storage variable
                                    output
                                        .entry(&inter.name().unwrap())
                                        .and_modify(
                                            |x| {
                                                x.update_end(end_p)
                                            }
                                        );
                                }
                            } else {
                                // extrapolated sequence's length does not exceed the stated thresholds
                                if codirected {
                                    start_p = b_q_start - offset;
                                    // assign to a storage variable
                                    output
                                        .entry(&inter.name().unwrap())
                                        .and_modify(
                                            |x| {
                                                x.update_start(start_p)
                                            }
                                        );
                                } else {
                                    end_p = b_q_end + offset;
                                    // assign to a storage variable
                                    output
                                        .entry(&inter.name().unwrap())
                                        .and_modify(
                                            |x| {
                                                x.update_end(end_p)
                                            }
                                        );
                                }
                            }
                        }
                    }
                } else {
                    for (mut i, inter) in intervals[curr..].iter().enumerate() {
                        i += curr;
                        // let block_id: String = i.to_string();
                        let inter_start: u64 = *inter.start().with_context(||
                            {format!("Interval {} has an undefined start coordinate which cannot be mapped", i)}
                        )?;
                        let inter_end: u64 = *inter.end().with_context(||
                            {format!("Interval {} has an undefined end coordinate which cannot be mapped", i)}
                        )?;
                        // add a results block to the the output hash map
                        if !output.contains_key(&inter.name().unwrap()) {
                            output.insert(
                                inter.name().unwrap(),
                                Interval::new()
                            );
                            output.
                                entry(&inter.name().unwrap())
                                .and_modify(
                                    |x| {
                                        x.update_name(inter.name().unwrap().to_string()); // TODO: Will borrow the value!
                                        x.update_chrom(self.query.chr.clone()); // TODO: Bad choice altogether
                                    }
                                );
                        }
        
                        // start coordinate is within the alignment gap
                        if (r_start <= inter_start) && (inter_start <= b_r_end) {
                            // get the alignment offset
                            let offset: u64 = b_r_end - inter_start;//inter_start - r_start;
                            // get the relative threshold size
                            let rel_thresh: &u64 = rel_sizes
                                .entry(
                                    inter.name().unwrap_or("a") // TODO: Find a way to create long-lived string literal IDs or update name() in cubiculum
                                )
                                .or_insert((inter.length().unwrap() as f64 * rel_threshold) as u64);
        
                            // check if the offset is within the stated extrapolation limits 
                            if offset > abs_threshold && offset > *rel_thresh {
                                // coordinate is too far to be extrapolated; crop to the chain block's start
                                if codirected {
                                    start_p = b_q_start;
                                    // assign to a storage variable
                                    output
                                        .entry(&inter.name().unwrap())
                                        .and_modify(
                                            |x| {
                                                x.update_start(start_p)
                                            }
                                        );
                                } else {
                                    end_p = b_q_end;
                                    // assign to a storage variable
                                    output
                                        .entry(&inter.name().unwrap())
                                        .and_modify(
                                            |x| {
                                                x.update_end(end_p)
                                            }
                                        );
                                }
                            } else {
                                // extrapolated sequence's length does not exceed the stated thresholds
                                if codirected {
                                    start_p = b_q_start - offset;
                                    // assign to a storage variable
                                    output
                                        .entry(&inter.name().unwrap())
                                        .and_modify(
                                            |x| {
                                                x.update_start(start_p)
                                            }
                                        );
                                } else {
                                    end_p = b_q_end + offset;
                                    // assign to a storage variable
                                    output
                                        .entry(&inter.name().unwrap())
                                        .and_modify(
                                            |x| {
                                                x.update_end(end_p)
                                            }
                                        );
                                }
                            }
                        }
        
                        // and the same for end coordinate
                        if (r_start <= inter_end) && (inter_end <= b_r_end) {
                            // get the alignment offset
                            let offset: u64 = inter_end - r_start;//r_block_end - inter_end;
                            // get the relative threshold size
                            let rel_thresh: &u64 = rel_sizes
                                .entry(
                                    inter.name().unwrap_or("a") // TODO: Find a way to create long-lived string literal IDs or update name() in cubiculum
                                )
                                .or_insert((inter.length().unwrap() as f64 * rel_threshold) as u64);
                            
                            // check if the offset is within the stated extrapolation limits 
                            if offset > abs_threshold && offset > *rel_thresh {
                                // coordinate is too far to be extrapolated; crop to the chain block's start
                                if codirected {
                                    end_p = b_q_end;
                                    // assign to a storage variable
                                    output
                                        .entry(&inter.name().unwrap())
                                        .and_modify(
                                            |x| {
                                                x.update_end(end_p)
                                            }
                                        );
                                } else {
                                    start_p = b_q_start;
                                    // assign to a storage variable
                                    output
                                        .entry(&inter.name().unwrap())
                                        .and_modify(
                                            |x| {
                                                x.update_start(start_p)
                                            }
                                        );
                                }
                            } else {
                                // extrapolated sequence's length does not exceed the stated thresholds
                                if codirected {
                                    end_p = b_q_start + offset;
                                    // assign to a storage variable
                                    output
                                        .entry(&inter.name().unwrap())
                                        .and_modify(
                                            |x| {
                                                x.update_end(end_p)
                                            }
                                        );
                                } else {
                                    start_p = b_q_end - offset;
                                    // assign to a storage variable
                                    output
                                        .entry(&inter.name().unwrap())
                                        .and_modify(
                                            |x| {
                                                x.update_start(start_p)
                                            }
                                        );
                                }
                            }
                        }
                    }
                }
            }

            // nothing to look past the last chain block; exit the outer for-loop
            if is_last_block {break};

            // if all the transcripts have been inspected, break the outer loop
            if curr >= intervals.len() {break};
            // update the absolute start of all the transcripts intervals
            min_start = *intervals[curr].start().with_context(||
                {format!("Interval {} has an undefined start coordinate which cannot be mapped", curr)}
            )?;
        }
        Ok(output)
    }

    /// NOTE: Functions with a trailing underscore implement 'real-time' chain body parsing;
    /// they are currently faster than their `yield_blocks()`-based counterparts,
    /// which we aim to change in futute
    /// 

    /// [YM]
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
    /// 
    /// `ignore_undefined` - Boolean flag indicating whether projections for intervals fully enclosed in aligned chain gaps should be left undefined
    /// 
    /// # Returns
    /// Result<&str, Interval> where each interval contains projected coordinates for each input interval 
    pub fn map_through_<'a, T>(
        &'a self, 
        // intervals: &mut Vec<(&str, u64, u64, &str)>,
        intervals: &'a mut Vec<T>,
        abs_threshold: u64,
        rel_threshold: f64,
        ignore_undefined: bool
    ) -> Result<FxHashMap<&'a str, Interval>> //Result<FxHashMap<&str, (u64, u64)>> 
    where 
        T: Coordinates + Named + Debug
    {
        // let output: FxHashMap<&str, (u64, u64)> = FxHashMap::default();
        let mut output: FxHashMap<&str, Interval> = FxHashMap::default();

        intervals.sort_by(
        |a, b| if a.start().unwrap() == b.start().unwrap() {
            a.end().unwrap().cmp(&b.end().unwrap())
        } else {
            a.start().unwrap().cmp(&b.start().unwrap())
        }
        );
        // define the total span of input intervals:
        // blocks before `min_start` will be ignored; 
        // once `max_end` is passed, iteration over chain stop 
        let mut min_start: u64 = *intervals[0].start().with_context(||
            {"Cannot map intervals with undefined coordinates"}
        )?;
        // note, however,  that the elements are sorted by the start coordinate alone,
        // so the last element must not necessarily end farthest
        let mut max_end: u64 = *intervals[intervals.len() - 1].end().with_context(||
            {"Cannot map intervals with undefined coordinates"}
        )?;
        // create a smart iteration index; iteration will always start from this interval
        let mut curr: usize = 0;
        // record the current interval's end coordinate; this will ensure that the iterator will never
        // skip the nested intervals
        let mut curr_end: u64 = *intervals[0].end().with_context(||
            {"Cannot map intervals with undefined coordinates"}
        )?;

        // create a hash map of relative length threshold; for long interval lists 
        // retrieving those from an array might be faster than calculating them every time anew
        let mut rel_sizes: FxHashMap<&str, u64> = FxHashMap::default();

        // define whether alignment is codirected between reference in query
        // for now we assume that chains always represent the positive strand in the reference sequence
        // this means, 'codirectionality' depends on the query strand alone
        let codirected: bool = &self.query.strand == &'+';

        // initialize the variables standing for block coordinates
        // (see TODO tho)
        // 
        let mut r_start: u64 = self.refs.start;
        let r_end: u64 = self.refs.end;
        let q_strand: bool = self.query.strand == '+';
        let mut q_start: u64 = match q_strand {
            true => self.query.start,
            false => self.query.size - self.query.start
        };
        let query_end  = match q_strand {
            true => self.query.end,
            false => self.query.size - self.query.start
        };
        let mut q_block_start: u64;
        let mut q_block_end: u64;

        // finally, initialize the projected coordinate variables
        let mut start_p: u64;
        let mut end_p: u64;

        // all set
        // now, iterate over alignment records
        // TODO: Implement to_blocks() as yielder to avoid code repetition
        'outer: for (h, b) in self.alignment.iter().enumerate() {
            // break if the iterator has passed beyond the last interval
            if r_start > max_end {
                // println!("All blocks passed; r_start={}, r_block_end={}", r_start, r_start + (b.size as u64));
                break
            };
            let mut r_block_end: u64 = r_start + (b.size as u64);
            // skip the block preceding the first interval's start in the reference
            if r_block_end + (b.dt as u64) < min_start {
                r_start += (b.size + b.dt) as u64;
                q_start = if q_strand {q_start + (b.size + b.dq) as u64} else {q_start - (b.size + b.dq) as u64};
                continue
            };

            // define the query coordinates
            if q_strand {
                q_block_start = q_start;
                q_block_end = q_block_start + (b.size as u64);
            } else {
                q_block_start = q_start - (b.size as u64);
                q_block_end = q_start;
            }

            // check if this is the last block
            let is_last_block: bool = (b.dt == 0) && (b.dq == 0);

            // now, we have a chain block with defined boundaries in both reference and query;
            // iterate over the intervals, check whether any of their coordinates can be projected 
            // through this block
            for (mut i, inter) in intervals[curr..].iter().enumerate() {
                i += curr;
                let inter_start: u64 = *inter.start().with_context(||
                    {format!("Interval {} has an undefined start coordinate which cannot be mapped", i)}
                )?;
                let inter_end: u64 = *inter.end().with_context(||
                    {format!("Interval {} has an undefined end coordinate which cannot be mapped", i)}
                )?;
                let inter_name = inter.name().with_context(||
                    {format!("Interval {} has an undefined name value; cannot assign projected coordinates", i)}
                )?;

                // add a results block to the the output hash map
                if !output.contains_key(&inter_name) {
                    output.insert(
                        inter.name().unwrap(),
                        Interval::new()
                    );
                    output.
                        entry(&inter_name)
                        .and_modify(
                            |x| {
                                x.update_name(inter_name.to_string()); // TODO: Will borrow the value!
                                x.update_chrom(self.query.chr.clone()); // TODO: Bad choice altogether
                            }
                        );
                }


                // chain block is upstream to the current interval;
                // since other are guaranteed to start at least in the same position,
                // the current loop can be safely exited
                if r_block_end < inter_start {
                    // println!("Bbbreakpoint!"); X
                    // the pointer can be updated here, but only if the next block is guaranteed to lie further 
                    // downstream to the previous interval;
                    // since the chain block are sorted and do not overlap, the easiest way to prove it
                    // is to check whether the current block's end does not end within the current interval group 
                    if r_block_end > curr_end {
                        // println!("Breakpoint pointer update: i={}, curr={}, b={}, r_start={}, r_block_end={}, inter_start={}, inter_end={}, inter_name={}", i, curr, h, r_start, r_block_end, inter_start, inter_end, inter_name); X
                        curr = i;
                        if curr >= intervals.len() {
                            // println!("All intervals covered; r_start={}, r_block_end={}", r_start, r_block_end); X
                            break 'outer
                        };
                    }
                    // potentially this is the farthest the intervals have ever reached 
                    // in terms of the  end coordinate; unless this boundary is exceeded, 
                    // the iteration start point will not be updated
                    if inter_end >= curr_end {
                        // curr = i;
                        curr_end = inter_end;
                    }
                    break
                }

                // chain block is downstream to the current interval;
                // nothing to do here, proceed to the next interval;
                if r_start > inter_end {
                    // println!("Continue point; i={}, curr={}, inter_name={}", i, curr, inter_name); X
                    // increase the pointer if this is the current leading interval
                    if inter_end == curr_end {
                        curr += 1;
                        if curr >= intervals.len() {
                            // println!("All intervals covered; r_start={}, r_block_end={}", r_start, r_block_end); X
                            break 'outer
                        };
                    }
                    continue
                };

                // first, process the marginal cases
                // a special case for the first block which extends beyond the chain start
                if h == 0 && inter_start < r_start && inter_end >= r_start {
                    let offset: u64 = r_start - inter_start;
                    // get the relative threshold size
                    let rel_thresh: &u64 = rel_sizes
                        .entry(
                            inter_name // TODO: Find a way to create long-lived string literal IDs or update name() in cubiculum
                        )
                        .or_insert((inter.length().unwrap() as f64 * rel_threshold) as u64);

                    // check if the offset is within the stated extrapolation limits 
                    if offset > abs_threshold && offset > *rel_thresh {
                        // coordinate is too far to be extrapolated; crop to the chain block's start
                        if codirected {
                            start_p = q_block_start;
                            // assign to a storage variable
                            output
                                .entry(&inter_name)
                                .and_modify(
                                    |x| {
                                        x.update_start(start_p)
                                    }
                                );
                        } else {
                            end_p = min(q_block_end, query_end);
                            // assign to a storage variable
                            output
                                .entry(&inter_name)
                                .and_modify(
                                    |x| {
                                        x.update_end(end_p)
                                    }
                                );
                        }
                    } else {
                        // extrapolated sequence's length does not exceed the stated thresholds
                        if codirected {
                            start_p = q_block_start.checked_sub(offset).unwrap_or(0);
                            // assign to a storage variable
                            output
                                .entry(&inter_name)
                                .and_modify(
                                    |x| {
                                        x.update_start(start_p)
                                    }
                                );
                        } else {
                            end_p = min(q_block_end + offset, query_end);
                            // assign to a storage variable
                            output
                                .entry(&inter_name)
                                .and_modify(
                                    |x| {
                                        x.update_end(end_p)
                                    }
                                );
                        }
                    }
                }

                // a special case for the last block; if interval end lies outside of the chain,
                // try extrapolating the coordinate unless it is too far from the chain 
                if is_last_block && inter_end >= r_end && inter_start < r_end {
                    // get the alignment offset
                    let offset: u64 = inter_end - r_block_end;
                    // get the relative threshold size
                    let rel_thresh: &u64 = rel_sizes
                        .entry(
                            inter_name // TODO: Find a way to create long-lived string literal IDs or update name() in cubiculum
                        )
                        .or_insert((inter.length().unwrap() as f64 * rel_threshold) as u64);
                    
                    // check if the offset is within the stated extrapolation limits 
                    if offset > abs_threshold && offset > *rel_thresh {
                        // coordinate is too far to be extrapolated; crop to the chain block's start
                        if codirected {
                            end_p = min(q_block_end, query_end);
                            // assign to a storage variable
                            output
                                .entry(&inter_name)
                                .and_modify(
                                    |x| {
                                        x.update_end(end_p)
                                    }
                                );
                        } else {
                            start_p = q_block_start;
                            // assign to a storage variable
                            output
                                .entry(&inter_name)
                                .and_modify(
                                    |x| {
                                        x.update_start(start_p)
                                    }
                                );
                        }
                    } else {
                        // extrapolated sequence's length does not exceed the stated thresholds
                        if codirected {
                            end_p = min(q_block_end + offset, query_end);
                            // assign to a storage variable
                            output
                                .entry(&inter_name)
                                .and_modify(
                                    |x| {
                                        x.update_end(end_p)
                                    }
                                );
                        } else {
                            start_p = q_block_start.checked_sub(offset).unwrap_or(0);
                            // assign to a storage variable
                            output
                                .entry(&inter_name)
                                .and_modify(
                                    |x| {
                                        x.update_start(start_p)
                                    }
                                );
                        }
                    }
                }

                // check whether the start coordinate is within the block
                if (r_start <= inter_start) && (inter_start < r_block_end) {
                    // println!("BLOCK: inter_start={}, r_start={}, r_block_end={}, q_block_start={}, q_block_end={}, i={}, inter_name={}", inter_start, r_start, r_block_end, q_block_start, q_block_end, i, inter_name); X
                    //  start coordinate can be mapped
                    let offset: u64 = inter_start - r_start;
                    if codirected{
                        start_p = q_block_start + offset;
                        // assign to a storage variable
                        output
                            .entry(&inter_name)
                            .and_modify(
                                |x| {
                                    x.update_start(start_p)
                                }
                            );
                    } else {
                        end_p = min(q_block_end - offset, query_end);
                        // assign to a storage variable
                        output
                            .entry(&inter_name)
                            .and_modify(
                                |x| {
                                    x.update_end(end_p)
                                }
                            );
                    }
                }
                // then, check the end coordinate
                if (r_start <= inter_end) && (inter_end < r_block_end) {
                    // println!("BLOCK: inter_end={}, r_start={}, r_block_end={}, q_block_start={}, q_block_end={}, codirected={}, i={}, inter_name={}", inter_end, r_start, r_block_end, q_block_start, q_block_end, codirected, i, inter_name); X
                    let offset: u64 = r_block_end - inter_end;
                    if codirected {
                        end_p = min(q_block_end.checked_sub(offset).unwrap_or(0), query_end);
                        // assign to a storage variable
                        output
                            .entry(&inter_name)
                            .and_modify(
                                |x| {
                                    x.update_end(end_p)
                                }
                            );
                    } else {
                        start_p = q_block_start + offset;
                        // assign to a storage variable
                        output
                            .entry(&inter_name)
                            .and_modify(
                                |x| {
                                    x.update_start(start_p)
                                }
                            );
                    }
                }

                curr_end = max(curr_end, inter_end);
                max_end = max(curr_end, max_end);
            }


            // nothing to look past the last chain block; exit the outer for-loop
            if is_last_block {break};

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

            // ignore blocks standing for full deletions in the reference
            if b.dt == 0 {
                r_start += b.dt as u64;
                q_start = if q_strand {q_start + b.dq as u64} else {q_start - b.dq as u64};
                continue
            }

            // current interval can be potentially exceeded at this point; exit if so
            if curr >= intervals.len() {
                // println!("All intervals covered; r_start={}, r_block_end={}", r_start, r_block_end); X
                break 'outer
            };


            // now, iterate through the remaining intervals again,
            // this time for the chain gap
            for (mut i, inter) in intervals[curr..].iter().enumerate() {
                i += curr;
                // keep track on coordinates enclosed in this gap for this interval
                let mut coords_in_gap: u8 = 0;
                let inter_start: u64 = *inter.start().with_context(||
                    {format!("Interval {} has an undefined start coordinate which cannot be mapped", i)}
                )?;
                let inter_end: u64 = *inter.end().with_context(||
                    {format!("Interval {} has an undefined end coordinate which cannot be mapped", i)}
                )?;
                let inter_name = inter.name().with_context(||
                    {format!("Interval {} has an undefined name value; cannot assign projected coordinates", i)}
                )?;
                // add a results block to the the output hash map
                if !output.contains_key(&inter.name().unwrap()) {
                    output.insert(
                        inter.name().unwrap(),
                        Interval::new()
                    );
                    output.
                        entry(&inter.name().unwrap())
                        .and_modify(
                            |x| {
                                x.update_name(inter.name().unwrap().to_string()); // TODO: Will borrow the value!
                                x.update_chrom(self.query.chr.clone()); // TODO: Bad choice altogether
                            }
                        );
                }

                // again, break if interval iterator has passed the current block 
                if r_block_end < inter_start {
                    // println!("Bbbreakpoint!");
                    // the pointer can be updated here, but only if the next block is guaranteed to lie further 
                    // downstream to the previous interval;
                    // since the chain block are sorted and do not overlap, the easiest way to prove it
                    // is to check whether the current block's end does not end within the current interval group 
                    if r_block_end > curr_end {
                        // println!("Breakpoint pointer update: i={}, curr={}, b={}, r_start={}, r_block_end={}, inter_start={}, inter_end={}, inter_name={}", i, curr, h, r_start, r_block_end, inter_start, inter_end, inter_name); X
                        curr = i;
                        if curr >= intervals.len() {
                            // println!("All intervals covered; r_start={}, r_block_end={}", r_start, r_block_end); X
                            break 'outer
                        };
                    }
                    // potentially this is the farthest the intervals have ever reached 
                    // in terms of the  end coordinate; unless this boundary is exceeded, 
                    // the iteration start point will not be updated
                    if inter_end >= curr_end {
                        // curr = i;
                        curr_end = inter_end;
                    }
                    break
                }

                // continue if the interval iterator has not yet reached the block
                if r_start > inter_end {
                    // println!("Continue point; i={}, curr={}, inter_name={}", i, curr, inter_name); X
                    // increase the pointer if this is the current leading interval
                    if inter_end == curr_end {
                        curr += 1;
                        if curr >= intervals.len() {
                            // println!("All intervals covered; r_start={}, r_block_end={}", r_start, r_block_end); X
                            break 'outer
                        };
                    }
                    continue
                };

                // start coordinate is within the alignment gap
                if (r_start <= inter_start) && (inter_start < r_block_end) {
                    // println!("GAP: inter_start={}, r_start={}, r_block_end={}, q_block_start={}, q_block_end={}, i={}, inter_name={}", inter_start, r_start, r_block_end, q_block_start, q_block_end, i, inter_name); X
                    coords_in_gap += 1;
                    // get the alignment offset
                    let offset: u64 = r_block_end - inter_start;//inter_start - r_start;
                    // get the relative threshold size
                    let rel_thresh: &u64 = rel_sizes
                        .entry(
                            inter.name().unwrap_or("a") // TODO: Find a way to create long-lived string literal IDs or update name() in cubiculum
                        )
                        .or_insert((inter.length().unwrap() as f64 * rel_threshold) as u64);

                    // check if the offset is within the stated extrapolation limits 
                    if offset > abs_threshold && offset > *rel_thresh {
                        // coordinate is too far to be extrapolated; crop to the chain block's start
                        if codirected {
                            start_p = q_block_start;
                            // assign to a storage variable
                            output
                                .entry(&inter.name().unwrap())
                                .and_modify(
                                    |x| {
                                        x.update_start(start_p)
                                    }
                                );
                        } else {
                            end_p = min(q_block_end, query_end);
                            // assign to a storage variable
                            output
                                .entry(&inter.name().unwrap())
                                .and_modify(
                                    |x| {
                                        x.update_end(end_p)
                                    }
                                );
                        }
                    } else {
                        // extrapolated sequence's length does not exceed the stated thresholds
                        if codirected {
                            start_p = q_block_end.checked_sub(offset).unwrap_or(0);
                            // assign to a storage variable
                            output
                                .entry(&inter.name().unwrap())
                                .and_modify(
                                    |x| {
                                        x.update_start(start_p)
                                    }
                                );
                        } else {
                            end_p = min(query_end, q_block_start + offset);
                            // assign to a storage variable
                            output
                                .entry(&inter.name().unwrap())
                                .and_modify(
                                    |x| {
                                        x.update_end(end_p)
                                    }
                                );
                        }
                    }
                }

                // and the same for end coordinate
                if (r_start <= inter_end) && (inter_end < r_block_end) {
                    // println!("GAP: inter_end={}, r_start={}, r_block_end={}, q_block_start={}, q_block_end={}, i={}, inter_name={}", inter_end, r_start, r_block_end, q_block_start, q_block_end, i, inter_name); X
                    coords_in_gap += 1;
                    if coords_in_gap == 2 && ignore_undefined {
                        output
                            .entry(&inter.name().unwrap())
                            .and_modify(
                                |x| {
                                    x.reset_start();
                                    x.reset_end()
                                }
                            );
                        // println!("RESETTING COORDS: {:#?}", output.get(&inter.name().unwrap()).unwrap());
                        continue
                    }
                    // get the alignment offset
                    let offset: u64 = inter_end - r_start;//r_block_end - inter_end;
                    // get the relative threshold size
                    let rel_thresh: &u64 = rel_sizes
                        .entry(
                            inter.name().unwrap_or("a") // TODO: Find a way to create long-lived string literal IDs or update name() in cubiculum
                        )
                        .or_insert((inter.length().unwrap() as f64 * rel_threshold) as u64);
                    
                    // check if the offset is within the stated extrapolation limits 
                    if offset > abs_threshold && offset > *rel_thresh {
                        // coordinate is too far to be extrapolated; crop to the chain block's start
                        if codirected {
                            end_p = min(q_block_end, query_end);
                            // assign to a storage variable
                            output
                                .entry(&inter.name().unwrap())
                                .and_modify(
                                    |x| {
                                        x.update_end(end_p)
                                    }
                                );
                        } else {
                            start_p = q_block_start;
                            // assign to a storage variable
                            output
                                .entry(&inter.name().unwrap())
                                .and_modify(
                                    |x| {
                                        x.update_start(start_p)
                                    }
                                );
                        }
                    } else {
                        // extrapolated sequence's length does not exceed the stated thresholds
                        if codirected {
                            end_p = min(q_block_start + offset, query_end);
                            // assign to a storage variable
                            output
                                .entry(&inter.name().unwrap())
                                .and_modify(
                                    |x| {
                                        x.update_end(end_p)
                                    }
                                );
                        } else {
                            start_p = q_block_end.checked_sub(offset).unwrap_or(0);
                            // assign to a storage variable
                            output
                                .entry(&inter.name().unwrap())
                                .and_modify(
                                    |x| {
                                        x.update_start(start_p)
                                    }
                                );
                        }
                    }
                }
                curr_end = max(curr_end, inter_end);
                max_end = max(curr_end, max_end);
            }

            // if all the transcripts have been inspected, break the outer loop
            if curr >= intervals.len() {
                // println!("All intervals covered; r_start={}, r_block_end={}", r_start, r_block_end); X
                break
            };
            // update the absolute start of all the transcripts intervals
            min_start = *intervals[curr].start().with_context(||
                {format!("Interval {} has an undefined start coordinate which cannot be mapped", curr)}
            )?;

            // proceed to the next line
            r_start += b.dt as u64;
            q_start = if q_strand {q_start + b.dq as u64} else {q_start - b.dq as u64};
        }
        Ok(output)
    }


    // [YM] + NOT FINISHED
    /// Maps coordinates from reference to query
    /// 
    /// # Arguments
    /// 
    /// `intervals` - A collection of objects having "start" and "end" coordinates; using tuples for nows
    /// 
    /// # Returns
    /// 
    /// Result<&str, u64> where key is each interval's name 
    /// and value is the number of bases covered by aligned blocks
    /// 
    /// 
    /// 
    pub fn alignment_cov_<'a, T>(&self, intervals: &'a mut Vec<T>,) -> Result<FxHashMap<&'a str, u64>> 
    where 
        T: Coordinates + Named + Debug
    {
        // the same routine as above
        // first, sort the input vector
        let mut output: FxHashMap<&str, u64> = FxHashMap::default();
        intervals.sort_by(
            |a, b| if a.start().unwrap() == b.start().unwrap() {
                a.end().unwrap().cmp(&b.end().unwrap())
            } else {
                a.start().unwrap().cmp(&b.start().unwrap())
            }
        );
        // define the total span for the input intervals
        let mut min_start: u64 = *intervals[0].start().with_context(||
            {"Cannot assess coverage for intervals with undefined coordinates"}
        )?;
        // note, however,  that the elements are sorted by the start coordinate alone,
        // so the last element must not necessarily end farthest
        let mut max_end: u64 = *intervals[intervals.len() - 1].end().with_context(||
            {"Cannot assess coverage for intervals with undefined coordinates"}
        )?;
        // create a smart iteration index; iteration will always start from this interval
        let mut curr: usize = 0;
        // record the current interval's end coordinate; this will ensure that the iterator will never
        // skip the nested intervals
        let mut curr_end: u64 = *intervals[0].end().with_context(||
            {"Cannot assess coverage for intervals with undefined coordinates"}
        )?;

        // create a smart iteration index; iteration will always start from this interval
        let mut curr: usize = 0;
        // record the current interval's end coordinate; this will ensure that the iterator will never
        // skip the nested intervals
        let mut curr_end: u64 = *intervals[0].end().with_context(||
            {"Cannot assess coverage for intervals with undefined coordinates"}
        )?;

        // initialize the variables standing for block coordinates
        // in this case, only the ref coordinates matter
        let mut r_start: u64 = self.refs.start;
        let mut r_block_end: u64 = 0;
        let r_end: u64 = self.refs.end;

        // now go
        for (h, b) in self.alignment.iter().enumerate() {
            r_block_end = r_start + b.size as u64;
            // continue if the first interval has not yet been reached
            if r_block_end < min_start {
                // don't forget to update the next block's start point
                r_start += (b.size + b.dt) as u64;
                continue
            };
            // break the block loop if the last interval has been passed
            if r_start > max_end {
                break
            };
            for (mut i, inter) in intervals[curr..].iter().enumerate() {
                i += curr;
                let inter_start: u64 = *inter.start().with_context(||
                    {format!("Interval {} has an undefined start coordinate which cannot be mapped", i)}
                )?;
                let inter_end: u64 = *inter.end().with_context(||
                    {format!("Interval {} has an undefined end coordinate which cannot be mapped", i)}
                )?;
                let name: &str = inter.name().with_context(||
                    {"Interval is not named"}
                )?;

                if !output.contains_key(&name) {
                    output.insert(
                        name,
                        0
                    );
                }

                // chain block is upstream to the current interval;
                // since other are guaranteed to start at least in the same position,
                // the current loop can be safely exited
                if r_block_end < inter_start {
                    // the pointer can be updated here, but only if the next block is guaranteed to lie further 
                    // downstream to the previous interval;
                    // since the chain block are sorted and do not overlap, the easiest way to prove it
                    // is to check whether the current block's end does not end within the current interval group 
                    if r_block_end >= curr_end {
                        curr = i
                    }
                    // potentially this is the farthest the intervals have ever reached 
                    // in terms of the  end coordinate; unless this boundary is exceeded, 
                    // the iteration start point will not be updated
                    if inter_end >= curr_end {
                        // curr = i;
                        curr_end = inter_end;
                    }
                    break
                }

                // chain block is downstream to the current interval;
                // nothing to do here, proceed to the next interval;
                if r_start > inter_end {
                    // if inter_end == curr_end {
                    //     curr += 1;
                    // }
                    continue
                };

                // current interval and current block intersect by at least 1 bp;
                // record their intersection
                if let Some(x) = intersection(inter_start, inter_end, r_start, r_block_end) {
                    // *output.get_mut(name).unwrap() += x;
                    output
                        .entry(name)
                        .and_modify(|y| *y += x)
                        .or_insert(0);
                }
                curr_end = max(curr_end, inter_end);
                max_end = max(curr_end, max_end);
            }
            // if the last interval has been passed after the inner for-loop, break the outer one
            if curr >= intervals.len() {println!("Last interval reached; r_start={}, r_block_end={}", r_start, r_block_end); break}
            // otherwise, update the next block's start point
            r_start += (b.size + b.dt) as u64;
            // and the first interval's start point
            min_start  = *intervals[curr].start().unwrap();
        }
        Ok(output)
    }
}