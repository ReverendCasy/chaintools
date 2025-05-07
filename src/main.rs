use anyhow::Context;
use chaintools as chain;
use cubiculum::structs::structs::Interval;
use fxhash::FxHashMap;
use std::process;
use std::time::Instant;

fn main() {
    let file:  &str = "/beegfs/projects/project-ymalovichko/toga_extension/duplication_tracing/TOGA2.0_tests/TREE_TESTS/mm10_iqtree_all_50/tmp/input_data/genome_alignment.chain";
    // let file: &str = "nanochain.chain";
    let acc_chains: Vec<&str> = vec!["1043254"];
    let chainmap: chain::cmap::map::ChainMap = chain::io::reader::Reader::extract(file, acc_chains).unwrap_or_else(|e| 
        {
            println!("File opening failed {}", e);
            process::exit(1);
        }
    );
    println!("Keys: {:?}", chainmap.map.keys());
    // let mode: chain::cmap::chain::BlockSide = chain::cmap::chain::BlockSide::Both;
    // let blocks: Vec<Box<dyn chain::cmap::chain::ChainBlock>> = chainmap
    //     .map[&1043254]
    //     .to_blocks(mode, false);
    // println!("Blocks: {:#?}", blocks);
    // println!("--------------------------------------------------\n");

    println!("Indexing test:");
    match chain::io::indexer::BinaryIndex::index(file) {
        Ok(()) => {println!("File was successfully indexed")},
        Err(err) => {println!("Indexing error: {:?}", err)}
    }
    println!("--------------------------------------------------\n");

    println!("Index-based extraction test");
    let chains_to_extract: Vec<u64> = vec![12, 38, 687002, 1043255, 1043256];
    let extracted_by_index: chain::map::ChainMap = chain::reader::Reader::extract_ix(file, Some(chains_to_extract)).unwrap_or_else(|err|
        {
            panic!("Index-based extraction failed!: {:?}", err);
        }
    );
    // let mode: chain::cmap::chain::BlockSide = chain::cmap::chain::BlockSide::Both;
    // println!("{:#?}", extracted_by_index.map[&1043255].header());
    // let blocks: Vec<Box <dyn chain::cmap::chain::ChainBlock>> = extracted_by_index
    //     .map[&1043255]
    //     .to_blocks(mode, true);
    // println!("Blocks: {:#?}", blocks);

    println!("Projecting coordinates");
    let mut vec_to_map: Vec<Interval> = vec![
        Interval::new()
    ];
    vec_to_map[0].update_name(String::from("ENST00000358005.7#CCDC32"));
    vec_to_map[0].update_chrom(String::from("chr15"));
    vec_to_map[0].update_start(40553970);
    vec_to_map[0].update_end(40554127);
    let mapped_coords: FxHashMap<&str, Interval> = extracted_by_index
        .map[&38]
        .map_through(&mut vec_to_map, 3000, 2.5)
        .expect("Failed mapping coordinates");
    println!("Mapped coordinates: {:#?}", mapped_coords);

    println!("Projecting a block protruding beyond the chain");
    let mut vec_to_map: Vec<Interval> = vec![
        Interval::new()
    ];
    vec_to_map[0].update_name(String::from("ENST00000216338.9#GZMH_utr2"));
    vec_to_map[0].update_chrom(String::from("chr14"));
    vec_to_map[0].update_start(24606479);
    vec_to_map[0].update_end(24606602);
    let mapped_coords: FxHashMap<&str, Interval> = extracted_by_index
        .map[&687002]
        .map_through(&mut vec_to_map, 3000, 2.5)
        .expect("Failed mapping coordinates");
    println!("Mapped coordinates for ENST00000216338.9#GZMH_utr2: {:#?}", mapped_coords);

    println!("Projecting a block partially covered by double-sided gap");
    let mut vec_to_map: Vec<Interval> = vec![
        Interval::new()
    ];
    vec_to_map[0].update_name(String::from("ENST00000693548.1#MYH15_partial_dside_gap"));
    vec_to_map[0].update_chrom(String::from("chr3"));
    vec_to_map[0].update_start(108454005);
    vec_to_map[0].update_end(108454142);
    let mapped_coords: FxHashMap<&str, Interval> = extracted_by_index
        .map[&12]
        .map_through(&mut vec_to_map, 3000, 2.5)
        .expect("Failed mapping coordinates");
    println!("ENST00000693548.1#MYH15_partial_dside_gap: {:#?}", mapped_coords);
    println!("--------------------------------------------------\n");

    println!("Timing legacy mapper");
    let mut vec_to_map: Vec<Interval> = vec![
        Interval::new()
    ];
    vec_to_map[0].update_name(String::from("ENST00000358005.7#CCDC32"));
    vec_to_map[0].update_chrom(String::from("chr15"));
    vec_to_map[0].update_start(40553970);
    vec_to_map[0].update_end(40554127);
    let now = Instant::now();
    let mapped_coords: FxHashMap<&str, Interval> = extracted_by_index
        .map[&38]
        .map_through_(&mut vec_to_map, 3000, 2.5)
        .expect("Failed mapping coordinates");
    let elapsed = now.elapsed();
    println!("Mapped coordinates: {:#?}", mapped_coords);
    println!("Elapsed time: {:?}", elapsed);
    println!();

    println!("Timing yield-based mapper");
    let now = Instant::now();
    let mapped_coords: FxHashMap<&str, Interval> = extracted_by_index
        .map[&38]
        .map_through(&mut vec_to_map, 3000, 2.5)
        .expect("Failed mapping coordinates");
    let elapsed = now.elapsed();
    println!("Mapped coordinates: {:#?}", mapped_coords);
    println!("Elapsed time: {:?}", elapsed);
    println!("--------------------------------------------------\n");


    println!("Checking intersection functionality");
    let mut vec_for_cov: Vec<Interval> = vec![
        Interval::new(),
        Interval::new(),
        Interval::new(),
        Interval::new(),
        Interval::new()
    ];
    vec_for_cov[0].update_chrom(String::from("chr14"));
    vec_for_cov[0].update_start(24606602);
    vec_for_cov[0].update_end(24606746);
    vec_for_cov[0].update_name(String::from("ENST00000216338.9#GZMH_1"));
    vec_for_cov[1].update_chrom(String::from("chr14"));
    vec_for_cov[1].update_start(24607148);
    vec_for_cov[1].update_end(24607406);
    vec_for_cov[1].update_name(String::from("ENST00000216338.9#GZMH_2"));
    vec_for_cov[2].update_chrom(String::from("chr14"));
    vec_for_cov[2].update_start(24607611);
    vec_for_cov[2].update_end(24607747);
    vec_for_cov[2].update_name(String::from("ENST00000216338.9#GZMH_3"));
    vec_for_cov[3].update_chrom(String::from("chr14"));
    vec_for_cov[3].update_start(24608264);
    vec_for_cov[3].update_end(24608412);
    vec_for_cov[3].update_name(String::from("ENST00000216338.9#GZMH_4"));
    vec_for_cov[4].update_chrom(String::from("chr14"));
    vec_for_cov[4].update_start(24609558);
    vec_for_cov[4].update_end(24609613);
    vec_for_cov[4].update_name(String::from("ENST00000216338.9#GZMH_5"));
    println!("vec_for_cov={:?}", vec_for_cov);

    println!("Timing legacy intersection");
    let now = Instant::now();
    let coverage: FxHashMap<&str, u64> = extracted_by_index
        .map[&687002]
        .alignment_cov_(&mut vec_for_cov)
        .expect("Failed calculating alignment coverage");
    let elapsed = now.elapsed();
    println!("Exon coverage for ENST00000216338.9#GZMH_utr2: {:#?}", coverage);
    println!("Elapsed: {:#?}", elapsed);

    println!("Timing yield-based intersection");
    let now = Instant::now();
    let coverage: FxHashMap<&str, u64> = extracted_by_index
        .map[&687002]
        .alignment_cov(&mut vec_for_cov)
        .expect("Failed calculating alignment coverage");
    let elapsed = now.elapsed();
    println!("Exon coverage for ENST00000216338.9#GZMH_utr2: {:#?}", coverage);
    println!("Elapsed: {:#?}", elapsed);


}