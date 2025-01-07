# rust-bedtools-multipangenome
 - rust bedtools multi-pangenome. 
 - analyzes the multiple pangenome alignments of the same query to the multiple reference genomes
 - limits the longest alignment search to the longest alignment.  
 - build a ancestral sequence based on the longest alignment. 
 - general note: Incase of Golang and RUST, please see the last commit message and if it says compiled binary then it is completed or else still in development version.


  ```
  cargo build 
  ```

  ```
  Usage: rust-bedtools-multipangenome <ALIGNMENT> <FASTAFILE> <PATHPRANK>

  Arguments:
  <ALIGNMENT>  please provide the path to the first alignment file
  <FASTAFILE>  please provide the reference fasta file
  <PATHPRANK>  please provide the path to the prank for the ancestal state

  Options:
  -h, --help     Print help
  -V, --version  Print version

  ```
  - to run the binary
  ```
  ./rust-bedtools-multipangenome \ 
      ./sample-files/sample.bed \
         ./sample-files/sample.fasta \ 
                 ./prank \
  ```

  Gaurav Sablok
