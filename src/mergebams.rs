/**

cd ~/develop/mergebams
cargo build --release
ml SAMtools
samtools view -h -s 0.0001 /home/sfurlan/scratch/MRB2/ITS_D383_3/outs/per_sample_outs/ITS_D383_3/count/sample_alignments.bam > ~/develop/mergebams/test/bam1.bam
samtools view -h -s 0.0001 /home/sfurlan/scratch/MRB2/ITS_D544_3/outs/per_sample_outs/ITS_D544_3/count/sample_alignments.bam > ~/develop/mergebams/test/bam2.bam
head /home/sfurlan/scratch/MRB2/ITS_D544_3/outs/per_sample_outs/ITS_D544_3/count/sample_barcodes.csv


cd ~/develop/mergebams/test
cargo build --release
../target/release/mergebams -t 4 -i bam1.bam,bam2.bam -l test1_,test2_ -b dummy -o .
samtools view out_bam.bam | head -n 200
**/
extern crate simple_log;
extern crate clap;
extern crate bam;
// extern crate parasailors;

// use simple_log::LogConfigBuilder;

use bam::RecordWriter;
use bam::record::tags::TagValue;
use clap::{App, load_yaml};

use std::str;
use std::io;
use std::io::BufWriter;

// #[derive(Clone)]


struct Params {
    inputs: String,
    labels: String,
    bcs: String,
    out: String,
    threads: usize,
}


fn load_params() -> Params {
    let yaml = load_yaml!("params_mergebams.yml");
    let params = App::from_yaml(yaml).get_matches();
    let inputs = params.value_of("inputs").unwrap();
    let labels = params.value_of("labels").unwrap();
    let threads: usize = params.value_of("threads").unwrap_or("1").parse().unwrap();
    let bcs = params.value_of("bcs").unwrap_or("1");
    let out = params.value_of("out").unwrap();
    Params{
        inputs: inputs.to_string(),
        labels: labels.to_string(),
        bcs: bcs.to_string(),
        out: out.to_string(),
        threads: threads,
    }
}

fn main() {
    let params = load_params();
    addtags(params);
    // combinebcs(&params);
    return;
    
}



fn addtags(params: Params) {
    let out_bam = params.out.to_string()+"/out_bam.bam";
    let s1 = params.inputs.to_string();
    let bamvec = s1.split(",").collect::<Vec<&str>>();
    let labvec1 = params.labels.to_string();
    let labvec = labvec1.split(",").collect::<Vec<&str>>();
    let (read_threads, write_threads) = if (*&params.threads as i8) > 2{
        (((*&params.threads/2) -1) as u16, ((*&params.threads/2) -1) as u16)
    } else {
        (0 as u16, 0 as u16)
    };
    let hreader = bam::BamReader::from_path(bamvec[0], read_threads).unwrap();
    let mut writer = bam::BamWriter::build()
        .write_header(true)
        .additional_threads(write_threads)
        .from_path(out_bam, hreader.header().clone()).unwrap();
    for (pos, inbam) in bamvec.iter().enumerate() {
        eprintln!("Opening: {:?}", inbam);
        let reader = bam::BamReader::from_path(inbam.to_string(), read_threads).unwrap();
        for record in reader {
            let mut newrecord = record.as_ref().unwrap().clone();
            match record.unwrap().tags().get(b"CB") {
                Some(TagValue::String(array_view, _)) => {
                    let labstr = labvec[pos];
                    let preftag = labstr.as_bytes().to_vec();
                    let oldtag = array_view.to_vec();
                    let newCB = &[preftag, oldtag].concat();
                    newrecord.tags_mut().remove(b"CB");
                    newrecord.tags_mut().push_string(b"CB", &newCB);
                    writer.write(&newrecord).unwrap();
                },
                Some(TagValue::Char(value)) => println!("Char = {}", value),
                _ => panic!("'CB' does not appear to have string type"),
            }
        }
    }
}
    

