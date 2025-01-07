use clap::Parser;

#[derive(Debug, Parser)]
#[clap(version)]
pub struct BedtoolsArgs {
    /// please provide the path to the first alignment file
    pub alignment: String,
    /// please provide the reference fasta file
    pub fastafile: String,
    /// please provide the path to the prank for the ancestal state
    pub pathprank: String,
}
