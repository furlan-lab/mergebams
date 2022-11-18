/**

cd ~/develop/mergebams
cargo build --release
ml SAMtools
samtools view -h -s 0.0001 /home/sfurlan/scratch/MRB2/ITS_D383_3/outs/per_sample_outs/ITS_D383_3/count/sample_alignments.bam > ~/develop/mergebams/test/bam1.bam
samtools view -h -s 0.0001 /home/sfurlan/scratch/MRB2/ITS_D544_3/outs/per_sample_outs/ITS_D544_3/count/sample_alignments.bam > ~/develop/mergebams/test/bam2.bam
head /home/sfurlan/scratch/MRB2/ITS_D544_3/outs/per_sample_outs/ITS_D544_3/count/sample_barcodes.csv


cd ~/develop/mergebams/test
cargo build --release
../target/release/mergebams -i bam1.bam,bam2.bam -l test1,test2 -b dummy -o .
**/
extern crate simple_log;
extern crate clap;
extern crate bam;
// extern crate parasailors;

// use simple_log::LogConfigBuilder;

use bam::RecordWriter;
use clap::{App, load_yaml};

use std::str;

// #[derive(Clone)]


struct Params {
    inputs: String,
    labels: String,
    bcs: String,
    out: String,
    threads: usize,
}


fn load_params() -> Params {
    let yaml = load_yaml!("params_count.yml");
    let params = App::from_yaml(yaml).get_matches();
    let inputs = params.value_of("inputs").unwrap();
    let labels = params.value_of("labels").unwrap();
    let threads: usize = params.value_of("threads").unwrap_or("1").parse().unwrap();
    let bcs = params.value_of("bcs").unwrap();
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
    let labvec = params.labels.to_string().split(",").collect::<Vec<&str>>();
    // println!("include header: {:?}", *&params.include_header);
    // println!("threads: {:?}", *&params.threads);
    let (read_threads, write_threads) = if (*&params.threads as i8) > 2{
        (((*&params.threads/2) -1) as u16, ((*&params.threads/2) -1) as u16)
    } else {
        (0 as u16, 0 as u16)
    };

    // make a header 
    // let mut header_reader = bam::BamReader::from_path(bamvec[0].to_string(), 0);
    // let header = bam::header::Header::from_bam(& mut header_reader).unwrap();
    let header = bam::header::Header::new();
    let mut writer = bam::BamWriter::build()
        .write_header(false)
        .additional_threads(write_threads)
        .from_path(out_bam, header).unwrap();
    for (pos, inbam) in bamvec.iter().enumerate() {
        println!("Opening: {:?}", inbam);
        let reader = bam::BamReader::from_path(inbam.to_string(), read_threads).unwrap();
        for record in reader {
            let mut newrecord = record.unwrap();
            let newtag = b"test";
            newrecord.tags_mut().push_string(b"CB", newtag);
            writer.write(&newrecord).unwrap();
        }
    }
}
    


// fn count_variants_mm(params: &Params, variant: Variant){
//     eprintln!("Processing using cb and umi in BAM tags");
//     let mut total: usize = 0;
//     let seqname = variant.seq;
//     let start = variant.start.parse::<u32>().unwrap();
//     let vname = variant.name;
//     let mut reader = bam::IndexedReader::build()
//         .additional_threads(*&params.threads as u16)
//         .from_path(&params.ibam).unwrap();
//     let mut seqnames = Vec::new();
//     let mut cb;
//     let mut umi;
//     let mut result = "null";
//     let query_nt = variant.query_nt as char;
//     let header = reader.header().clone();
//     let data = header.reference_names();
//     for seq in data {
//         seqnames.push(seq)
//     }
//     let ref_id = seqnames.iter().position(|&r| r == &seqname).unwrap();
//     let region = process_variant(ref_id as u32, start);
//     for record in reader.fetch_by(&&region, |record| record.mapq() > 4 && (record.flag().all_bits(0 as u16) || record.flag().all_bits(16 as u16))).unwrap(){
//         total+=1;
//         match record.as_ref().unwrap().tags().get(b"CB") {
//             Some( bam::record::tags::TagValue::String(cba, _)) => {
//                 cb = str::from_utf8(&cba).unwrap().to_string();
//             },
//             _ => panic!("Unexpected type"),
//         }
//         match record.as_ref().unwrap().tags().get(b"UB") {
//             Some( bam::record::tags::TagValue::String(uba, _)) => {
//                 // assert!(string == string);
//                 umi = str::from_utf8(&uba).unwrap().to_string();
//                 // write!(writer, cb);
//             },
//             _ => panic!("Unexpected type"),
//         }
//         for entry in record.as_ref().unwrap().alignment_entries().unwrap() {
//             if let Some((ref_pos, ref_nt)) = entry.ref_pos_nt() {
//                 if region.start() == ref_pos {
//                     // if entry.is_insertion() || entry.is_insertion(){
//                     //     println!("{}", "Indel");
//                     // }
//                     if let Some((_record_pos, record_nt)) = entry.record_pos_nt() {
//                         // let result = match record_nt as char {
//                         //     ref_nt as char => "ref",
//                         //     query_nt as char => "query",
//                         //     _ => "other",
//                         // };
//                         if ref_nt as char == record_nt as char {
//                             result = "ref";
//                         } else if record_nt as char == query_nt{
//                             result = "query";
//                         } else {
//                             result = "other";
//                         }
//                             println!("{} {} {} {} {} {}",&cb, &umi, seqname, ref_pos, vname, result);
//                         }
//                     } else {
//                         continue
//                     }
//             } else {
//                 continue
//             }        }
//     }
//     eprintln!("Found {} reads spanning this variant!", total);

// }
