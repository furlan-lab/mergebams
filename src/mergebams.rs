
/**

*********create test files*********

cd ~/develop/mergebams
cargo build --release
ml SAMtools
samtools view -b -h -s 0.0001 /home/sfurlan/scratch/MRB2/ITS_D383_3/outs/per_sample_outs/ITS_D383_3/count/sample_alignments.bam > ~/develop/mergebams/test/bam1.bam
samtools view -b -h -s 0.0001 /home/sfurlan/scratch/MRB2/ITS_D544_3/outs/per_sample_outs/ITS_D544_3/count/sample_alignments.bam > ~/develop/mergebams/test/bam2.bam
samtools view -b -h -s 0.0001 /home/sfurlan/scratch/MRB2/ITS_D544_3/outs/per_sample_outs/ITS_D544_3/count/sample_alignments.bam > ~/develop/mergebams/test/bam3.bam
samtools view -b -h -s 0.0001 /home/sfurlan/scratch/NB1/AML_601_34/outs/per_sample_outs/AML_601_34/count/sample_alignments.bam > ~/develop/mergebams/test/bam4.bam


head /home/sfurlan/scratch/MRB2/ITS_D544_3/outs/per_sample_outs/ITS_D544_3/count/sample_barcodes.csv

************************************


*********installation of mergebams*************
* 
* 
* *********************************************

*********usage of mergebams*********

cd ~/develop/mergebams/test
cargo build --release
ls../target/release/mergebams -i bam1.bam,bam2.bam -l test1_,test2_ -b barcodes1.tsv.gz,barcodes2.tsv.gz -o .
../target/release/mergebams -i bam2.bam,bam3.bam -l test1_,test2_ -b barcodes1.tsv.gz,barcodes2.tsv.gz -o .
../target/release/mergebams -i bam3.bam,bam4.bam -l test1_,test2_ -b barcodes1.tsv.gz,barcodes2.tsv.gz -o .
../target/release/mergebams -i bam1.bam,bam2.bam,bam3.bam,bam4.bam -l test1_,test2_,test3_,test4_ -b barcodes1.tsv.gz,barcodes2.tsv.gz,barcodes3.tsv.gz,barcodes4.tsv.gz -o .
../target/release/mergebams -i bam1.bam,bam2.bam,bam3.bam -l test1_,test2_,test3_ -b barcodes1.tsv.gz,barcodes2.tsv.gz,barcodes3.tsv.gz -o .


zcat < out_barcodes.tsv.gz
samtools view out_bam.bam | head -n 200
samtools sort out_bam.bam > out_bam.sorted.bam
samtools view out_bam.sorted.bam | head -n 200



************************************


**/


extern crate simple_log;
extern crate clap;
extern crate bam;
extern crate csv;
extern crate flate2;
extern crate differ;
extern crate itertools;

use itertools::Itertools;
use differ::{Differ, Tag};
use flate2::GzBuilder;
use flate2::Compression;
use bam::RecordWriter;
use bam::record::tags::TagValue;
use clap::{App, load_yaml};
use std::str;
use flate2::{read};
use std::{
    // error::Error,
    ffi::OsStr,
    fs::File,
    io::{self, BufReader, BufRead, Write},
    path::Path,
};


struct Params {
    inputs: String,
    labels: String,
    bcs: String,
    out: String,
    threads: usize,
}


fn main() {
    let params = load_params();
    // let header_result = checkheaders(params);
    if let Ok((header, params)) = checkheaders(params){
        let params = addtags(params, header);
        if params.bcs == "0"{
            return;
        } else {
            let _result = lines_from_file(params);
            return;
        }
    }else{
        eprintln!("ERROR: BAM header sequences do not match - you will need to fix this before merging bams");
    }

}


fn load_params() -> Params {
    let yaml = load_yaml!("params_mergebams.yml");
    let params = App::from_yaml(yaml).get_matches();
    let inputs = params.value_of("inputs").unwrap();
    let labels = params.value_of("labels").unwrap();
    let threads: usize = params.value_of("threads").unwrap_or("1").parse().unwrap();
    let bcs = params.value_of("bcs").unwrap_or("0");
    let out = params.value_of("out").unwrap();
    Params{
        inputs: inputs.to_string(),
        labels: labels.to_string(),
        bcs: bcs.to_string(),
        out: out.to_string(),
        threads: threads,
    }
}


fn lines_from_file(params: Params) -> io::Result<()>{
    let out_csv = params.out.to_string()+"/out_barcodes.tsv.gz";
    let s1 = params.bcs.to_string();
    let tsvvec = s1.split(",").collect::<Vec<&str>>();
    let labvec1 = params.labels.to_string();
    let labvec = labvec1.split(",").collect::<Vec<&str>>();
    let f = File::create(out_csv)?;
    let mut gz = GzBuilder::new()
                    .filename("/out_barcodes.tsv.gz")
                    .write(f, Compression::default());

    for (pos, filename) in tsvvec.iter().enumerate() {
        let path = Path::new(filename);
        let file = match File::open(&path) {
            Err(_why) => panic!("couldn't open {}", path.display()),
            Ok(file) => file,
        };
        if path.extension() == Some(OsStr::new("gz")){
            let buf = BufReader::new(read::GzDecoder::new(file));
            for line in buf.lines(){
                writeln!(&mut gz, "{}", labvec[pos].to_string() + &line.unwrap())?;
            }
        }else{
            let buf = BufReader::new(file);
            for line in buf.lines(){
                writeln!(&mut gz, "{}", labvec[pos].to_string() + &line.unwrap())?;
            }
        }
    }
    return Ok(());
}


fn checkheaders(params: Params) -> Result<(bam::Header, Params), &'static str>{
    let inputs = params.inputs.to_string();
    let bam_vec = inputs.split(",").collect::<Vec<&str>>();
    let bam_vec_to_header = inputs.split(",").collect::<Vec<&str>>();
    let bam_it = bam_vec.into_iter().combinations(2);
    let mut grand_count = 0;
    for bam_pairs in bam_it {
        let mut header_vec = Vec::new();
        for inbam in bam_pairs.iter() {
            let hreader = bam::BamReader::from_path(inbam, 0).unwrap();
            let header = hreader.header().clone();
            header_vec.push(header);
        }
        let mut count = 0;
        let a = header_vec[0].reference_names();
        let b = header_vec[1].reference_names();
        let differ = Differ::new(&a, &b);
        for span in differ.spans() {
            match span.tag {
                Tag::Equal => (),
                _ => count+=1,
            }
        }
        eprintln!("Found {} discrepencies between sequence names in the header of:\n\n{}\nand:\n\n{}\n", count, bam_pairs[0], bam_pairs[1]);
        grand_count+=count;
    }
    if grand_count == 0{
        let new_header = make_new_header(bam_vec_to_header);
        return Ok((new_header, params));
    } else {
        return Err("Discrepent headers")
    }
    
}

fn make_new_header(bam_vec: Vec<&str>) -> bam::Header {
    // assumes no discrepencies in the headers across bams in bam_vec
    let hreader = bam::BamReader::from_path(bam_vec[0], 0).unwrap();
    let mut out_header = hreader.header().clone();
    let mergebam_line = bam_vec.join(", ");
    let _msg = out_header.push_line(&("@CO\tmergebams has included the BAM records from the following files (using the header from the first): ".to_owned()+&mergebam_line));

    // for (pos, inbam) in bam_vec.iter().enumerate()  {
    //     // let hreader = bam::BamReader::from_path(inbam, 0).unwrap();
    //     // let header = hreader.header();
    //     // eprintln!("{}, {}", pos, inbam);
    //     let hreader = bam::BamReader::from_path(inbam, 0).unwrap();
    //     let header = hreader.header();
    //     eprintln!("{:?}", header.reference_names());
    //     // if pos == 0 {
    //     //     let mut hreader = std::io::Read::read(in_bam).unwrap();
    //     //     out_header = bam::Header::from_bam(& mut hreader).unwrap();
    //     //     // for line in header.lines() {
    //     //     //     eprintln!("{:?}", line);
    //     //     // }
    //     //     // let mut header_line = HeaderEntry::header_line("1.6".to_string());
    //     //     // header_line.push(b"SO", "Coordinate".to_string());
    //     //     // header.push_entry(header_line).unwrap();
    //     //     // // Reference line       "@SQ  SN:chr1  LN:10000".
    //     //     // header.push_entry(HeaderEntry::ref_sequence("chr1".to_string(), 10000)).unwrap();
 
    //     // }
    // }
    return out_header;
}

fn addtags(params: Params, header: bam::Header) -> Params{
    let out_bam = params.out.to_string()+"/out_bam.bam";
    let fail_bam = params.out.to_string()+"/fail_bam.bam";
    let out_bam_msg = out_bam.clone();
    let inputs = params.inputs.to_string();
    let bam_vec = inputs.split(",").collect::<Vec<&str>>();
    let bam_vec_msg = inputs.split(",").join(" and ");
    let labels = params.labels.to_string();
    let lab_vec = labels.split(",").collect::<Vec<&str>>();
    let (read_threads, write_threads) = if (*&params.threads as i8) > 2{
        (((*&params.threads/2) -1) as u16, ((*&params.threads/2) -1) as u16)
    } else {
        (0 as u16, 0 as u16)
    };
    let mut fail_count = 0;
    let mut pass_count = 0;
    let mut other_count = 0;
    // let hreader = bam::BamReader::from_path(bam_vec[0], read_threads).unwrap();
    let mut pass_writer = bam::BamWriter::build()
        .write_header(true)
        .additional_threads(write_threads)
        .from_path(out_bam, header.clone()).unwrap();
    let mut fail_writer = bam::BamWriter::build()
        .write_header(true)
        .additional_threads(0)
        .from_path(fail_bam, header.clone()).unwrap();
    eprintln!("Headers ok\nWriting:\n\n{}\nfrom:\n\n{}\n", out_bam_msg, bam_vec_msg);
    for (pos, inbam) in bam_vec.iter().enumerate() {
        let reader = bam::BamReader::from_path(inbam.to_string(), read_threads).unwrap();
        for record in reader {
            let mut newrecord = record.as_ref().unwrap().clone();
            match record.unwrap().tags().get(b"CB") {
                Some(TagValue::String(array_view, _)) => {
                    let labstr = lab_vec[pos];
                    let preftag = labstr.as_bytes().to_vec();
                    let oldtag = array_view.to_vec();
                    let new_cb = &[preftag, oldtag].concat();
                    newrecord.tags_mut().remove(b"CB");
                    newrecord.tags_mut().push_string(b"CB", &new_cb);
                    pass_writer.write(&newrecord).unwrap();
                    pass_count+=1;
                },
                Some(TagValue::Char(value)) => {
                    let labstr = lab_vec[pos];
                    let preftag = labstr.as_bytes().to_vec();
                    let oldtag = value.to_string().as_bytes().to_vec();
                    let new_cb = &[preftag, oldtag].concat();
                    newrecord.tags_mut().remove(b"CB");
                    newrecord.tags_mut().push_string(b"CB", &new_cb);
                    pass_writer.write(&newrecord).unwrap();
                    other_count+=1;
                },
                _ => {
                    // eprintln!("ERROR: 'CB' not found");
                    fail_writer.write(&newrecord).unwrap();
                    fail_count+=1;
                }
            }
        }
    }
    eprintln!("Processed all reads!!\nFound:\n{} - reads PASSING\n{} - reads PASSING but with issues\n{} - reads FAILING", pass_count, other_count, fail_count);
    return params;
}
    

