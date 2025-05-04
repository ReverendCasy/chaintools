//!
//! 'chaintools' is a library for serializing and working with genomic chain files in Rust.
//!
//! This crate provides two main entries:
//!
//! 1. `Reader` and `Writer` for reading and writing encoded chain files.
//! 2.  `ChainMap` and inner structs for performing basic operations (filter, sort, etc) on chain files.
//!
//! ## Reading chain files
//!
//! The entrypoint for reading genomic chain files in the 'chaintools' Rust library is through
//! the Reader struct, offering methods like [`Reader::from_file`] to create a reader from a file path
//! and return a [`Result`] with a [`ChainMap`], [`Reader::from_bytes`] to create a reader from a byte slice,
//! [`Reader::from_bin`] to create a reader from a binary file path, and [`Reader::load_chain`] to load a
//! specific chain from a binary file using a chain ID. Internally, the Reader struct handles file
//! operations, parsing byte data efficiently, and utilizes Rust's standard library for file I/O
//! and parallel processing for parsing and filtering chain data, providing a comprehensive interface
//! for working with genomic chain files in various formats within Rust.
//!
//! # Examples
//! ```rust
//! use chaintools::cmap::ChainMap;
//! use chaintools::Reader;
//!
//! let reader: ChainMap = Reader::from_file("path/to/chainfile.chain")?;
//! ```
//!
//! ## Operations over chain files
//!
//! The ChainMap struct provides a variety of methods for performing operations on chain files,
//! including filtering, sorting, inserting, removing, and iterating over chains, making it very
//! flexible to work with. Operations like [`ChainMap::filter`] and [`ChainMap::iter`] allow users
//! to filter and iterate over chains in a chain file, while [`ChainMap::get`] and [`ChainMap::insert`]
//! allow users to get and insert chains into a chain file. The ChainMap struct also provides specific
//! filtering methods like [`ChainMap::filter_id`] and [`ChainMap::filter_score`] to filter chains by ID
//! and score, respectively. This set of methods is inteded to grow over time as more operations are
//! needed for working with genomic chain files in Rust (e.g., merging, splitting, etc).
//!
//! On top of these operations, the Chain struct provides a simple interface for working with individual
//! chains, allowing users to access chain data like ID, score, and blocks, as well as perform basic string
//! operations mainly, such as accessing chain information or converting chain data to a string. A
//! Chain is composed by two inner structs, ChainHead and AlignmentRecord, which represent the chain
//! header and alignment block, respectively, and provide methods for accessing and manipulating chain
//! data. The basic picture of the choosen way to represent a .chain file in this crate is
//! basically a json-like structure, where each chain is a set of inner structs, and each inner struct
//! is a set of fields, which can be accessed and manipulated by the user.
//!
//! # Examples
//! ```rust
//! use chaintools::cmap::ChainMap;
//! use chaintools::Reader;
//! use chaintools::cmap::Chain;
//!
//! let reader: ChainMap = Reader::from_file("path/to/chainfile.chain")?;
//! let chain: Chain = reader.get(10);
//! ```
//!
//! ## Writing chain files
//!
//! The Writer struct provides a simple interface for writing genomic chain files in Rust, offering
//! a discrete set of methods in this release, such as [`Writer::to_bin`] to write a encoded chain
//! file and [`Writer::to_bin_gz`] to write an encoded and gz compressed chain file in cases were
//! the size of the file is a concern.
//!
//! # Examples
//! ```rust
//! use chaintools::cmap::ChainMap;
//! use chaintools::Reader;
//! use chaintools::Writer;
//!
//! let chains: ChainMap = Reader::from_file("path/to/chainfile.chain")?;
//! let filt_chains = chains.filter(|c| c.score > 1000);
//! Writer::to_bin(&filt_chains, "path/to/chainfile.chain");
//! ```

#![warn(rust_2021_compatibility)]
#![warn(rust_2018_idioms)]

pub mod cmap;
pub mod io;

pub use crate::io::*;
pub use crate::cmap::*;
