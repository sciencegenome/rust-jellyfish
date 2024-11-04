use clap::Parser;

#[derive(Debug, Parser)]
#[clap(version)]

pub struct KmeroriginArgs {
    /// please provide the kmer to be searched for the origin
    pub kmer_arg: usize,
    /// please provide the path to be searched for the strings containing the kmer
    pub fastqfile_arg: String,
}
