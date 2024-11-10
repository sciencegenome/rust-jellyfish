# rust-jellyfish

- a rust implementation of the jellyfish for the counts. 
- outputs both the unique counts, all counts. 
- directly plots the histograms from the counts using the rust plotters.
- it will produce allkmers, uniquekmers, countkmers
- available as crate [rust-jellyfish](https://crates.io/crates/kmerjellyfish)
```
Usage: kmerjellyfish <KMER_ARG> <FASTQFILE_ARG>

Arguments:
  <KMER_ARG>       please provide the kmer to be searched for the origin
  <FASTQFILE_ARG>  please provide the path to be searched for the strings containing the kmer

Options:
  -h, --help     Print help
  -V, --version  Print version

kmerjellyfish 4 ./sample-files/test.fastq
```

Gaurav Sablok
