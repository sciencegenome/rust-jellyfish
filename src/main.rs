mod args;
use args::JellyfishArgs;
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

    let args = JellyfishArgs::parse();
    kmer_jellyfish(&args.fastqfile_arg, args.kmer_arg);
}

fn kmer_jellyfish(path: &str, kmer:usize) {

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
            if line.starts_with("A") || line.starts_with("T") || line.starts_with("G") || line.starts_with("C") {
                sequence.push(&line)
            }  
            let mut sequence_iter:Vec<&str> = vec![];
            for i in 0..sequence.len() {
                let i = sequence[i]; 
                for j in 0..i.len() - &kmer {
                    sequence_iter.push(&i[j..j+&kmer])
            }
            let mut fileall = File::create("allkmerunique.txt").expect("file not present");
            for i in sequence_iter.iter() {
                write!(fileall, "{}\n", i.to_string());
            }
           let hash_kmer: HashSet<_> = sequence_iter.iter().collect();
          let mut finalvec: Vec<&str> = Vec::new();
          for i in hash_kmer.into_iter() {
                finalvec.push(i)
            }
            let mut path = File::create("kmerunique.txt").expect("file not present");
            for i in finalvec.iter() {
                write!(path, "{}\n", i.to_string());
          }
           let mut mutunique = Vec::new(); 
           let uniquehold = File::open("kmerunique.txt").expect("file not present");
           let uniqueread = BufReader::new(uniquehold);
           for i in uniqueread.lines() {
             let appendline = i.expect("line not present");
             mutunique.push(appendline)
           }
        
           struct CountIllumina {
            kmer: String,
            count:usize,
           }
           let mut allkmer = Vec::new(); 
           let allopen = File::open("allkmerunique.txt").expect("file not present");
           let allread = BufReader::new(allopen);
           for i in allread.lines() {
            let appendline = i.expect("line not present");
            allkmer.push(appendline)
           }
           let mut mutunique = Vec::new(); 
           let uniquehold = File::open("kmerunique.txt").expect("file not present");
           let uniqueread = BufReader::new(uniquehold);
           for i in uniqueread.lines() {
             let appendline = i.expect("line not present");
             mutunique.push(appendline)
           }
        
           let mut unique_count = Vec::new();
           for i in mutunique.iter() {
           let count_add = allkmer.iter().filter(|&x| *x == *i).count();
           unique_count.push(CountIllumina{
           kmer: i.to_string(),
           count: count_add,
           })
         }
         let mut unique_file = File::create("histogram-count.txt").expect("file not present");
         for i in unique_count.iter(){
              write!(unique_file, "{}\t{}\n",i.kmer,i.count);
         }
        
           #[derive(Debug)]
           struct VecStore {
            id: String,
            numberstart: usize,
            numberend : usize,
           }
           let mut indexstorestart = Vec::new();   
            for i in sequence.iter() {
                for j in mutunique.iter() {
                    let indexout = i.find(j).unwrap();
                    let indexoutend = i.find(j).unwrap()+j.len();
                    indexstorestart.push(VecStore {id: i.to_string(), 
                          numberstart:indexout, 
                          numberend:indexoutend});
                }
            }
        }
        }
        }
