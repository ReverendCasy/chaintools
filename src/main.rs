use chaintools as chain;
use std::process;

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
    let mode: chain::cmap::chain::BlockSide = chain::cmap::chain::BlockSide::Both;
    let blocks: Vec<chain::cmap::chain::ChainBlock> = chainmap
        .map[&1043254]
        .to_blocks(mode, false);
    println!("Blocks: {:#?}", blocks);
    println!("--------------------------------------------------\n");

    println!("Indexing test:");
    match chain::io::indexer::BinaryIndex::index(file) {
        Ok(()) => {println!("File was successfully indexed")},
        Err(err) => {println!("Indexing error: {:?}", err)}
    }
    println!("--------------------------------------------------\n");

    println!("Index-based extraction test");
    let chains_to_extract: Vec<u64> = vec![1043255, 1043256];
    let extracted_by_index: chain::map::ChainMap = chain::reader::Reader::extract_ix(file, chains_to_extract).unwrap_or_else(|err|
        {
            panic!("Index-based extraction failed!: {:?}", err);
        }
    );
    let mode: chain::cmap::chain::BlockSide = chain::cmap::chain::BlockSide::Both;
    println!("{:#?}", extracted_by_index.map[&1043255].header());
    let blocks: Vec<chain::cmap::chain::ChainBlock> = extracted_by_index
        .map[&1043255]
        .to_blocks(mode, true);
    println!("Blocks: {:#?}", blocks);

}