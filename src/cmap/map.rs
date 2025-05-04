use fxhash::FxHashMap;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::cmap::chain::Chain;

/// A map of chains
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChainMap {
    pub map: FxHashMap<u32, Chain>,
}

impl ChainMap {
    /// Create a new ChainMap
    ///
    /// # Returns
    /// * ChainMap
    ///
    /// # Example
    /// ```
    /// use chaintools::cmap::ChainMap;
    ///
    /// let chains = ChainMap::new();
    /// ```
    pub fn new() -> Self {
        Self {
            map: FxHashMap::default(),
        }
    }

    /// Get a chain from the map
    ///
    /// # Arguments
    /// * `key` - A chain id
    ///
    /// # Returns
    /// * Option<&Chain>
    ///
    /// # Example
    /// ```
    /// use chaintools::cmap::ChainMap;
    /// use chaintools::io::reader::Reader;
    ///
    /// let chains = Reader::from_file("file.chain").unwrap();
    /// let chain = chains.get(&123);
    /// println!("{:?}", chain);
    ///
    /// >>> Some(Chain { score: 0, refs: ChainHead { size: 0, start: 0, end: 0, strand: 0 },
    /// query: ChainHead { size: 0, start: 0, end: 0, strand: 0 }, alinment: [], id: 123 })
    /// ```
    pub fn get(&self, key: &u32) -> Option<&Chain> {
        self.map.get(key)
    }

    /// Get a mutable chain from the map
    ///
    /// # Arguments
    /// * `key` - A chain id
    ///
    /// # Returns
    /// * Option<&mut Chain>
    ///
    /// # Example
    /// ```
    /// use chaintools::cmap::ChainMap;
    /// use chaintools::io::reader::Reader;
    ///
    /// let chains = Reader::from_file("file.chain").unwrap();
    /// let chain = chains.get_mut(&123);
    /// println!("{:?}", chain);
    ///
    /// >>> Some(Chain { score: 0, refs: ChainHead { size: 0, start: 0, end: 0, strand: 0 },
    /// query: ChainHead { size: 0, start: 0, end: 0, strand: 0 }, alinment: [], id: 123 })
    /// ```
    pub fn get_mut(&mut self, key: &u32) -> Option<&mut Chain> {
        self.map.get_mut(key)
    }

    /// Generalized filtering function
    ///
    /// # Arguments
    /// * `fc` - A closure that takes a reference to a Chain and returns a boolean
    ///
    /// # Returns
    /// * ChainMap
    ///
    /// # Example
    /// ```
    /// use chaintools::cmap::ChainMap;
    /// use chaintools::cmap::chain::Chain;
    /// use chaintools::io::reader::Reader;
    ///
    /// let chains = Reader::from_file("file.chain").unwrap();
    /// let filt_chains = chains.filter(|chain| chain.score > 100);
    /// ```
    pub fn filter<F>(&self, fc: F) -> Self
    where
        F: Fn(&Chain) -> bool,
    {
        let fmap: FxHashMap<u32, Chain> = self
            .map
            .iter()
            .filter(|(_, chain)| fc(chain))
            .map(|(id, chain)| (*id, chain.clone()))
            .collect();

        Self { map: fmap }
    }

    /// Insert a chain into a ChainMap
    ///
    /// # Arguments
    /// * `key` - A chain id
    /// * `value` - A Chain
    ///
    /// # Returns
    /// * &mut ChainMap
    ///
    /// # Example
    /// ```
    /// use chaintools::cmap::ChainMap;
    /// use chaintools::io::reader::Reader;
    ///
    /// let chains = Reader::from_file("file.chain").unwrap();
    /// chains.insert(123, Chain{ score: 0, refs: ChainHead { size: 0, start: 0, end: 0, strand: 0 }, ..., id: 123});
    /// ```
    pub fn insert(&mut self, key: u32, value: Chain) -> &mut ChainMap {
        self.map.insert(key, value);
        self
    }

    /// Remove a chain from a ChainMap
    ///
    /// # Arguments
    /// * `key` - A chain id
    ///
    /// # Returns
    /// * &mut ChainMap
    ///
    /// # Example
    /// ```
    /// use chaintools::cmap::ChainMap;
    /// use chaintools::io::reader::Reader;
    ///
    /// let chains = Reader::from_file("file.chain").unwrap();
    /// chains.remove(&123);
    /// ```
    pub fn remove(&mut self, key: &u32) -> &mut ChainMap {
        self.map.remove(key);
        self
    }

    /// Sort the ChainMap by chain ids
    ///
    /// # Returns
    /// * usize
    ///
    /// # Example
    /// ```
    /// use chaintools::cmap::ChainMap;
    /// use chaintools::io::reader::Reader;
    ///
    /// let chains = Reader::from_file("file.chain").unwrap();
    /// println!("Number of chains: {:?}", chains.len());
    ///
    /// >>> 100
    /// ```
    pub fn len(&self) -> usize {
        self.map.len()
    }

    /// Iterate over the ChainMap
    ///
    /// # Returns
    /// * Iterator<Item = (&u32, &Chain)>
    ///
    /// # Example
    /// ```
    /// use chaintools::cmap::ChainMap;
    /// use chaintools::io::reader::Reader;
    ///
    /// let chains = Reader::from_file("file.chain").unwrap();
    /// for (key, chain) in chains.iter() {
    ///   println!("{:?} {:?}", key, chain);
    ///   break;
    /// }
    /// ```
    pub fn iter(&self) -> impl Iterator<Item = (&u32, &Chain)> {
        self.map.iter()
    }

    /// Iterate over mut ChainMap
    ///
    /// # Returns
    /// * Iterator<Item = (&u32, &mut Chain)>
    ///
    /// # Example
    /// ```
    /// use chaintools::cmap::ChainMap;
    /// use chaintools::io::reader::Reader;
    ///
    /// let chains = Reader::from_file("file.chain").unwrap();
    /// for (key, chain) in chains.iter_mut() {
    ///  println!("{:?} {:?}", key, chain);
    ///  break;
    /// }
    /// ```
    pub fn iter_mut(&mut self) -> impl Iterator<Item = (&u32, &mut Chain)> {
        self.map.iter_mut()
    }

    /// Iterate over the keys of the ChainMap
    ///
    /// # Returns
    /// * Iterator<Item = &u32>
    ///
    /// # Example
    /// ```
    /// use chaintools::cmap::ChainMap;
    /// use chaintools::io::reader::Reader;
    ///
    /// let chains = Reader::from_file("file.chain").unwrap();
    /// for key in chains.keys() {
    ///    println!("{:?}", key);
    /// }
    /// ```
    pub fn keys(&self) -> impl Iterator<Item = &u32> {
        self.map.keys()
    }

    /// Iterate over the values of the ChainMap
    ///
    /// # Returns
    /// * Iterator<Item = &Chain>
    ///
    /// # Example
    /// ```
    /// use chaintools::cmap::ChainMap;
    /// use chaintools::io::reader::Reader;
    ///
    /// let chains = Reader::from_file("file.chain").unwrap();
    /// for chain in chains.values() {
    ///   println!("{:?}", chain);
    /// }
    /// ```
    pub fn values(&self) -> impl Iterator<Item = &Chain> {
        self.map.values()
    }

    /// Iterate over mut values of the ChainMap
    ///
    /// # Returns
    /// * Iterator<Item = &mut Chain>
    ///
    /// # Example
    /// ```
    /// use chaintools::cmap::ChainMap;
    /// use chaintools::io::reader::Reader;
    ///
    /// let chains = Reader::from_file("file.chain").unwrap();
    /// for chain in chains.values_mut() {
    ///     println!("{:?}", chain);
    ///     chain.score = 100;
    ///     println!("{:?}", chain);
    ///     break;
    ///  }
    /// ```
    pub fn values_mut(&mut self) -> impl Iterator<Item = &mut Chain> {
        self.map.values_mut()
    }

    /// Filter the ChainMap by score
    ///
    /// # Arguments
    /// * `score` - An i64
    ///
    /// # Returns
    /// * ChainMap
    ///
    /// # Example
    /// ```
    /// use chaintools::cmap::ChainMap;
    /// use chaintools::io::reader::Reader;
    ///
    /// let chains = Reader::from_file("file.chain").unwrap();
    /// let score = 1000;
    /// let filt_chains = chains.filter_by_score(score);
    /// ```
    pub fn filter_by_score(&self, score: u64) -> Self {
        let map = self
            .map
            .par_iter()
            .filter(|(_, v)| v.score >= score)
            .map(|(k, v)| (*k, v.clone()))
            .collect();
        Self { map }
    }

    /// Filter the ChainMap by reference size
    ///
    /// # Arguments
    /// * `size` - A u64
    ///
    /// # Returns
    /// * ChainMap
    ///
    /// # Example
    /// ```
    /// use chaintools::cmap::ChainMap;
    /// use chaintools::io::reader::Reader;
    ///
    /// let chains = Reader::from_file("file.chain").unwrap();
    /// let size = 1000;
    /// let filt_chains = chains.filter_query_by_size(size);
    /// ```
    pub fn filter_ref_by_size(&self, size: u64) -> Self {
        let map = self
            .map
            .par_iter()
            .filter(|(_, v)| v.refs.size >= size)
            .map(|(k, v)| (*k, v.clone()))
            .collect();
        Self { map }
    }

    /// Filter the ChainMap by query size
    ///
    /// # Arguments
    /// * `size` - A u64
    ///
    /// # Returns
    /// * ChainMap
    ///
    /// # Example
    /// ```
    /// use chaintools::cmap::ChainMap;
    /// use chaintools::io::reader::Reader;
    ///
    /// let chains = Reader::from_file("file.chain").unwrap();
    /// let size = 1000;
    /// let filt_chains = chains.filter_query_by_size(size);
    /// ```
    pub fn filter_query_by_size(&self, size: u64) -> Self {
        let map = self
            .map
            .par_iter()
            .filter(|(_, v)| v.query.size >= size)
            .map(|(k, v)| (*k, v.clone()))
            .collect();
        Self { map }
    }

    /// Filter the ChainMap by chain ids
    ///
    /// # Arguments
    /// * `ids` - A vector of chain ids
    ///
    /// # Returns
    /// * ChainMap
    ///
    /// # Example
    /// ```
    /// use chaintools::cmap::ChainMap;
    /// use chaintools::io::reader::Reader;
    ///
    /// let chains = Reader::from_file("file.chain").unwrap();
    /// let ids = vec![123, 456];
    /// let filt_chains = chains.filter(ids).len();
    /// println!("{:?}", filt_chains);
    ///
    /// >>> 2
    pub fn filter_id(&self, ids: Vec<u32>) -> Self {
        let map = self
            .map
            .par_iter()
            .filter(|(k, _)| ids.contains(k))
            .map(|(k, v)| (*k, v.clone()))
            .collect();
        Self { map }
    }
}
