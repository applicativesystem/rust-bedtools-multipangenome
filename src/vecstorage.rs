/*
 * holding all the structs in the separate files so that they
 * can be easily called as a reference call in the result.
 *
 *
 * */
#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct BedtoolsRange {
    pub alignedref: String,
    pub alignedstart: usize,
    pub alignedend: usize,
    pub difference: usize,
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Fasta {
    pub header: String,
    pub sequence: String,
}
