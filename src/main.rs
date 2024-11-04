mod args;
use args::KmeroriginArgs;
use clap::Parser;
use std::fs::File;
use std::io::{Write, BufReader, BufRead};
use std::collections::HashSet;
#[allow(dead_code)]

/*
 *Author Gaurav Sablok
 *Universitat Potsdam
 *Date 2024-11-4

 * jellyfish counter for the illumina reads for the miseq, nextseq.
 * it outputs the unique kmers as the all profiles kmers and the
 * respective count for making the histograms for the profiled jellyfish counts.
 *
 * */

fn main() {

    let args = KmeroriginArgs::parse();
    kmer_fastq(args.fastqfile_arg, args.kmer_arg);
}

fn kmer_fastq(path: String, kmer:usize) {
    let f = File::open(&path).expect("file not present");
    let read = BufReader::new(f);
    for i in read.lines() {
        let line = i
                   .expect("line not present");
        let mut header: Vec<&str> = vec![];
        let mut sequence:Vec<&str> = vec![];
        if line.starts_with("@") {
            header.push(&line)
        }
        if line.starts_with("A") || line.starts_with("T") || line.starts_with("G") || line.starts_with("C")  {
            sequence.push(&line)
        }
        // illumina struct for faster analysis
        struct Illuminaseq {
        id : String,
        seq: String,
        }

        let illumina = Vec::new();
        for i in 0..header.len(){
            illumina.push(Illuminaseq{
            id:header[i].to_string(),
            seq:sequence[i].to_string(),
            })
        }
        let mut sequence_iter:Vec<&str> = vec![];
        for i in 0..sequence.len() {
            let i = sequence[i];
            for j in 0..i.len() - &kmer {
                sequence_iter.push(&i[j..j+&kmer])
        }
       let mut append = File::create("kmerallfastq.txt").expect("file not present");
       for i in sequence_iter.iter(){
        write!(append, "{}\n", i.to_string());
       }
       println!("file for the all kmers profiled has been written");
       let hash_kmer: HashSet<_> = sequence_iter.iter().collect();
       let mut finalvec: Vec<&str> = Vec::new();
       for i in hash_kmer.into_iter() {
            finalvec.push(i)
        }
        let mut path = File::create("kmeruniquefasta.txt").expect("file not present");
        for i in finalvec.iter() {
            write!(path, "{}\n", i.to_string());
       }
       println!("file for the unique kmers have been written");
       let mut mutunique = Vec::new();
       let uniquehold = File::open("kmeruniquefasta.txt").expect("file not present");
       let uniqueread = BufReader::new(uniquehold);
       for i in uniqueread.lines() {
         let appendline = i.expect("line not present");
         mutunique.push(appendline)
       }

       // opening a new struct for the count and the count is implemented with in the struct.
       #[derive(Debug)]
       struct CountIllumina {
        kmer: String,
        count:usize,
       }
        let mut unique_count = Vec::new();
        for i in mutunique.iter() {
        let count_add = sequence_iter.iter().filter(|&x| *x == i).count();
        unique_count.push(CountIllumina{
        kmer: i.to_string(),
        count: count_add,
        })
      }
      let mut unique_file = File::create("histogram-count.txt");
      for i in unique_count.iter(){
        let id = i.kmer.to_string();
        let kmer = i.count.to_string();
      write!(unique_file, "{}\t{}",id,kmer)}

       // a vec struct holding the other count implementation.
       #[derive(Debug)]
       struct VecStore {
        id: String,
        kmer: String,
        numberstart: usize,
        numberend : usize,
       }
       let mut indexstorestart = Vec::new();
        for i in sequence.iter() {
            for j in mutunique.iter() {
                let indexout = i.find(j).unwrap();
                let indexoutend = i.find(j).unwrap()+j.len();
                indexstorestart.push(VecStore {id: i.to_string(),
                      kmer: j.to_string(),
                      numberstart:indexout,
                      numberend:indexoutend});
            }
        }
          for line in indexstorestart.iter() {
           println!("{}\t{}\t{}\t{}", line.id, line.kmer, line.numberstart, line.numberend)
        }
    }
    for j in header.iter() {
        println!("{}", j)
    }
        }

}
