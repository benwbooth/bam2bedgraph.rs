use linked_hash_map::LinkedHashMap;
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::BTreeSet;
use std::hash::{Hash, Hasher};
use regex::Regex;
use regex::Captures;
use url::percent_encoding::{percent_decode, utf8_percent_encode, SIMPLE_ENCODE_SET};
use bio::data_structures::interval_tree::IntervalTree;
use bio::utils::Interval;
use bio::alphabets::dna;
use std;
use std::fs::File;
use std::io::{BufReader, BufWriter, BufRead, Write};
use std::io::stdout;
use std::ops::Range;
use std::path::Path;
use ::errors::*;
use itertools::Itertools;
use unindent::unindent;

#[derive(Default, Clone)]
pub struct Record {
    pub row: usize,
    pub seqname: String,
    pub source: String,
    pub feature_type: String,
    pub start: u64,
    pub end: u64,
    pub score: String,
    pub strand: String,
    pub frame: String,
    pub attributes: LinkedHashMap<String, String>,
}

impl PartialEq for Record {
    fn eq(&self, other: &Self) -> bool {
        self.seqname == other.seqname && self.source == other.source &&
        self.feature_type == other.feature_type && self.start == other.start &&
        self.end == other.end && self.score == other.score &&
        self.strand == other.strand && self.frame == other.frame
    }
}
impl Eq for Record {}

impl Hash for Record {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.seqname.hash(state);
        self.source.hash(state);
        self.feature_type.hash(state);
        self.start.hash(state);
        self.end.hash(state);
        self.score.hash(state);
        self.strand.hash(state);
        self.frame.hash(state);
    }
}

impl Record {
    pub fn new() -> Record {
        Record { ..Default::default() }
    }
    pub fn from_row(row: usize, line: &str, filetype: &str, chrmap: &HashMap<String,String>) -> Result<Record> {
        lazy_static! {
            static ref COMMENT: Regex = Regex::new(r"^#").unwrap();
            static ref GTF_ATTR: Regex = Regex::new(r#"^(?P<key>\S+)\s+(?:"(?P<qval>[^"]*)"|(?P<val>\S+));\s*"#).unwrap();
        }
        if COMMENT.is_match(line) {
            return Err("Comment".into());
        }
        let fields: Vec<_> = line.split('\t').collect();
        let seqname = String::from(*fields.get(0).unwrap_or(&""));
        let seqname = chrmap.get(&seqname).unwrap_or(&seqname);
        Ok(Record {
            row: row,
            seqname: seqname.clone(),
            source: String::from(*fields.get(1).unwrap_or(&"")),
            feature_type: String::from(*fields.get(2).unwrap_or(&"")),
            start: fields.get(3).unwrap_or(&"0").parse::<u64>()?,
            end: fields.get(4).unwrap_or(&"0").parse::<u64>()?,
            score: String::from(*fields.get(5).unwrap_or(&"")),
            strand: String::from(*fields.get(6).unwrap_or(&"")),
            frame: String::from(*fields.get(7).unwrap_or(&"")),
            attributes: if filetype == "gff" {
                fields.get(8).unwrap_or(&"").split(';').map(|a| {
                    let kv: Vec<_> = a.splitn(2, '=').collect();
                    let key = String::from(percent_decode(kv[0].as_bytes()).decode_utf8().unwrap());
                    let value = String::from(percent_decode(kv.get(1).unwrap_or(&"").as_bytes()).decode_utf8().unwrap());
                    (key, value)
                }).collect()
            } else if filetype == "gtf" {
                GTF_ATTR.captures_iter(fields.get(8).unwrap_or(&"")).map(|caps| {
                    (String::from(caps["key"].to_string()), 
                     String::from(caps.name("qval").unwrap_or_else(|| caps.name("val").unwrap()).as_str().to_string()))
                }).collect()
            } else {
                return Err(format!("Don't know how to read filetype {}", filetype).into());
            },
        })
    }
    
    pub fn to_gtf(&self) -> Result<String> {
        Ok(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
            self.seqname, 
            self.source, 
            self.feature_type,
            self.start,
            self.end,
            self.score,
            self.strand,
            self.frame,
            self.attributes.iter().map(|(k,v)| format!("{} \"{}\";", k, v)).join(" ")))
        
    }
    
    pub fn to_gff(&self) -> Result<String> {
        define_encode_set! {
            pub GFF_ENCODE_SET = [SIMPLE_ENCODE_SET] | {'\t', '\r', '\n', ';', '%', '='}
        }
        Ok(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
            self.seqname, 
            self.source, 
            self.feature_type,
            self.start,
            self.end,
            self.score,
            self.strand,
            self.frame,
            self.attributes.iter().map(|(k,v)| 
                format!("{}={}", 
                    utf8_percent_encode(k, GFF_ENCODE_SET), 
                    utf8_percent_encode(v, GFF_ENCODE_SET))).join(";")))
    }
}

pub struct IndexedAnnotation {
    pub rows: Vec<Record>,
    pub id2row: HashMap<String, usize>,
    pub row2parents: HashMap<usize, Vec<usize>>,
    pub row2children: HashMap<usize, Vec<usize>>,
    pub tree: HashMap<String, IntervalTree<u64, usize>>,
    pub chrmap: HashMap<String, String>,
    pub vizchrmap: HashMap<String, String>,
    pub refs: LinkedHashMap<String, u64>,
}

impl IndexedAnnotation {
    pub fn from_gtf(
        annotfile: &str, 
        gene_type: &str, 
        transcript_type: &str,
        chrmap_file: &Option<String>,
        vizchrmap_file: &Option<String>) 
        -> Result<IndexedAnnotation> 
    {
        IndexedAnnotation::from_file(annotfile, "gtf", gene_type, transcript_type, chrmap_file, vizchrmap_file)
    }
    pub fn from_gff(
        annotfile: &str, 
        chrmap_file: &Option<String>,
        vizchrmap_file: &Option<String>) 
        -> Result<IndexedAnnotation> 
    {
        IndexedAnnotation::from_file(annotfile, "gff", "", "", chrmap_file, vizchrmap_file)
    }
    pub fn from_file(
        annotfile: &str, 
        filetype: &str, 
        gene_type: &str, 
        transcript_type: &str, 
        chrmap_file: &Option<String>,
        vizchrmap_file: &Option<String>) 
        -> Result<IndexedAnnotation> 
    {
        // read in the optional chr map
        let mut chrmap = HashMap::<String,String>::new();
        if let Some(charmap_file) = chrmap_file.clone() {
            let f = File::open(&charmap_file)?;
            let mut file = BufReader::new(&f);
            let mut line = String::new();
            while file.read_line(&mut line)? > 0 {
                let line = line.trim_right_matches('\n').trim_right_matches('\r');
                let cols: Vec<&str> = line.split('\t').collect();
                if let Some(key) = cols.get(0) {
                    if let Some(value) = cols.get(1) {
                        chrmap.insert(String::from(*key), String::from(*value));
                    }
                }
            }
        }
        // read in the optional visualization chr map
        let mut vizchrmap = HashMap::<String,String>::new();
        if let Some(vizcharmap_file) = vizchrmap_file.clone() {
            let f = File::open(&vizcharmap_file)?;
            let mut file = BufReader::new(&f);
            let mut line = String::new();
            while file.read_line(&mut line)? > 0 {
                let line = line.trim_right_matches('\n').trim_right_matches('\r');
                let cols: Vec<&str> = line.split('\t').collect();
                if let Some(key) = cols.get(0) {
                    if let Some(value) = cols.get(1) {
                        vizchrmap.insert(String::from(*key), String::from(*value));
                    }
                }
            }
        }
        
        let mut id2row = HashMap::<String, usize>::new();
        let mut rows = Vec::<Record>::new();
        let f = File::open(&annotfile)?;
        let mut file = BufReader::new(&f);
        let mut refs = HashMap::<String,u64>::new();
        let mut line = String::new();
        while file.read_line(&mut line)? > 0 {
            let line = line.trim_right_matches('\n').trim_right_matches('\r');
            let row = rows.len();
            if let Ok(record) = Record::from_row(row, &line, filetype, &chrmap) {
                if let Some(id) = record.attributes.get("ID") {
                    id2row.insert(id.clone(), row);
                }
                // get the max ref lengths
                let mut reflength = refs.entry(record.seqname.clone()).or_insert(record.end);
                if *reflength < record.end { *reflength = record.end }
                rows.push(record);
            }
        }
        // sort the refs
        let mut sorted_refs = LinkedHashMap::<String,u64>::new();
        let chrs = refs.keys().collect::<Vec<_>>();
        for chr in chrs {
            let length = refs[chr];
            sorted_refs.insert(chr.clone(), length);
        }
        
        // populate row2parents and row2children indices
        let mut row2parents = HashMap::<usize, BTreeSet<usize>>::new();
        let mut row2children = HashMap::<usize, BTreeSet<usize>>::new();
        // make fake rows for GTF gene_id and transcript_id attributes
        let mut fake_rows = Vec::<Record>::new();
        
        match filetype {
            "gtf" => {
                let mut firstrow = HashMap::<&Record, usize>::new();
                for (row, record) in rows.iter().enumerate() {
                    // in GTF, if two rows both have identical fields (except attributes), they
                    // are the same feature.
                    let row = *firstrow.entry(record).or_insert(row);
                    if let Some(transcript_id) = record.attributes.get("transcript_id") {
                        let mut transcript_record = Record::new();
                        transcript_record.feature_type = String::from(transcript_type);
                        transcript_record.attributes
                            .insert(String::from("ID"), transcript_id.clone());
                        if let Some(transcript_name) = record.attributes.get("transcript_name") {
                            transcript_record.attributes
                                .insert(String::from("Name"), transcript_name.clone());
                        }
                        let transcriptrow = rows.len() + fake_rows.len();
                        id2row.insert(transcript_id.clone(), transcriptrow);
                        row2parents.entry(row).or_insert_with(BTreeSet::new).insert(transcriptrow);
                        row2children.entry(transcriptrow)
                            .or_insert_with(BTreeSet::new)
                            .insert(row);
                        fake_rows.push(transcript_record);

                        if let Some(gene_id) = record.attributes.get("gene_id") {
                            let mut gene_record = Record::new();
                            gene_record.feature_type = String::from(gene_type);
                            gene_record.attributes.insert(String::from("ID"), gene_id.clone());
                            if let Some(gene_name) = record.attributes.get("gene_name") {
                                gene_record.attributes
                                    .insert(String::from("Name"), gene_name.clone());
                            }
                            let generow = rows.len() + fake_rows.len();
                            id2row.insert(gene_id.clone(), generow);
                            row2parents.entry(transcriptrow)
                                .or_insert_with(BTreeSet::new)
                                .insert(generow);
                            row2children.entry(generow)
                                .or_insert_with(BTreeSet::new)
                                .insert(transcriptrow);
                            fake_rows.push(gene_record);
                        }
                    }
                }
            }
            "gff" => {
                for (row, record) in rows.iter().enumerate() {
                    if let Some(parent) = record.attributes.get("Parent") {
                        for p in parent.split(',') {
                            if let Some(parentrow) = id2row.get(p) {
                                row2parents.entry(row)
                                    .or_insert_with(BTreeSet::new)
                                    .insert(*parentrow);
                                row2children.entry(*parentrow)
                                    .or_insert_with(BTreeSet::new)
                                    .insert(row);
                            }
                        }
                    }
                }
            }
            _ => (),
        };
        rows.append(&mut fake_rows);

        // fill in missing information, build interval tree
        let mut updates = Vec::<(usize, Record)>::new();
        for (row, record) in rows.iter().enumerate() {
            // if the seqname is missing then the record is not complete
            if record.seqname.is_empty() {
                // recursively get all child nodes
                let mut children = HashSet::<usize>::new();
                let mut addchildren = BTreeSet::<usize>::new();
                if let Some(cs) = row2children.get_mut(&row) {
                    addchildren.append(cs);
                }
                while !addchildren.is_empty() {
                    let mut morechildren = BTreeSet::<usize>::new();
                    for ac in &mut addchildren.iter() {
                        if !children.contains(ac) {
                            children.insert(*ac);
                            if let Some(cs) = row2children.get_mut(ac) {
                                morechildren.append(cs);
                            }
                        }
                    }
                    addchildren.clear();
                    addchildren.append(&mut morechildren);
                }
                // aggregate the child nodes to get the missing information
                let mut seqname = HashMap::<String, u64>::new();
                let mut source = HashMap::<String, u64>::new();
                let mut start: Option<u64> = Option::None;
                let mut end: Option<u64> = Option::None;
                let mut strand = HashMap::<String, u64>::new();
                for c in &children {
                    let child = &rows[*c];
                    if start.is_none() || child.start < start.r()? {
                        start = Some(child.start);
                    }
                    if end.is_none() || child.end > end.r()? {
                        end = Some(child.end);
                    }
                    *seqname.entry(child.seqname.clone()).or_insert(0u64) += 1;
                    *source.entry(child.source.clone()).or_insert(0u64) += 1;
                    *strand.entry(child.strand.clone()).or_insert(0u64) += 1;
                }
                // enqueue the missing information
                let mut update = Record::new();
                if let Some((seqname, _)) = seqname.iter().max_by_key(|a| a.1) {
                    update.seqname = seqname.clone();
                }
                if let Some((source, _)) = source.iter().max_by_key(|a| a.1) {
                    update.source = source.clone();
                }
                if let Some(start) = start {
                    update.start = start;
                }
                if let Some(end) = end {
                    update.end = end;
                }
                if let Some((strand, _)) = strand.iter().max_by_key(|a| a.1) {
                    update.strand = strand.clone();
                }
                updates.push((row, update));
            }
        }
        // update the missing information
        for u in &mut updates {
            let update = &mut u.1;
            let mut record = &mut rows[u.0];
            record.seqname = update.seqname.clone();
            record.source = update.source.clone();
            record.start = update.start;
            record.end = update.end;
            record.strand = update.strand.clone();
        }
        // sort row2parents by coordinates
        let mut row2parents_update = HashMap::<usize, Vec<usize>>::new();
        for (row, parents) in row2parents {
            let mut up = row2parents_update.entry(row).or_insert_with(Vec::new);
            let mut vec: Vec<_> = parents.into_iter().collect();
            vec.sort_by(|a, b| {
                rows[*a]
                    .seqname.as_bytes()
                    .cmp(&rows[*b].seqname.as_bytes())
                    .then_with(|| rows[*a].start.cmp(&rows[*b].start))
            });
            *up = vec;
        }
        // sort row2children by coordinates
        let mut row2children_update = HashMap::<usize, Vec<usize>>::new();
        for (row, children) in row2children {
            let mut up = row2children_update.entry(row).or_insert_with(Vec::new);
            let mut vec: Vec<_> = children.into_iter().collect();
            vec.sort_by(|a, b| {
                rows[*a]
                    .seqname.as_bytes()
                    .cmp(&rows[*b].seqname.as_bytes())
                    .then_with(|| rows[*a].start.cmp(&rows[*b].start))
            });
            *up = vec;
        }
        // build the interval tree
        let mut tree = HashMap::<String, IntervalTree<u64, usize>>::new();
        for (row, record) in rows.iter().enumerate() {
            let mut t = tree.entry(record.seqname.clone()).or_insert_with(IntervalTree::new);
            t.insert(Interval::new((record.start - 1)..(record.end))?, row);
        }
        
        Ok(IndexedAnnotation {
            rows: rows,
            id2row: id2row,
            row2children: row2children_update,
            row2parents: row2parents_update,
            tree: tree,
            chrmap: chrmap,
            vizchrmap: vizchrmap,
            refs: sorted_refs,
        })
    }
    
    pub fn to_gtf(&self, filename: &str) -> Result<()> {
        let mut output: BufWriter<Box<Write>> = BufWriter::new(
            if filename == "-" { Box::new(stdout()) } 
            else { Box::new(File::create(filename)?) });
            
        for (row, record) in self.rows.iter().enumerate() {
            if let Some(transcript_rows) = self.row2parents.get(&row) {
                for transcript_row in transcript_rows {
                    let transcript = &self.rows[*transcript_row];
                    if let Some(transcript_id) = transcript.attributes.get("ID") {
                        if let Some(gene_rows) = self.row2parents.get(transcript_row) {
                            for gene_row in gene_rows {
                                let gene = &self.rows[*gene_row];
                                if let Some(gene_id) = gene.attributes.get("ID") {
                                    let mut rec = Record {
                                        row: row,
                                        seqname: record.seqname.clone(),
                                        source: record.source.clone(),
                                        feature_type: record.feature_type.clone(),
                                        start: record.start,
                                        end: record.end,
                                        score: record.score.clone(),
                                        strand: record.strand.clone(),
                                        frame: record.frame.clone(),
                                        attributes: LinkedHashMap::new(),
                                    };
                                    rec.attributes.insert(String::from("gene_id"), gene_id.clone());
                                    rec.attributes.insert(String::from("transcript_id"), transcript_id.clone());
                                    for (k,v) in &record.attributes {
                                        if k != "gene_id" && k != "transcript_id" {
                                            rec.attributes.insert(k.clone(), v.clone());
                                        }
                                    }
                                    output.write_fmt(format_args!("{}\n", rec.to_gtf()?))?;
                                }
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }
    
    pub fn to_gff(&self, filename: &str) -> Result<()> {
        let mut output: BufWriter<Box<Write>> = BufWriter::new(
            if filename == "-" { Box::new(stdout()) } 
            else { Box::new(File::create(filename)?) });
            
        for (row, record) in self.rows.iter().enumerate() {
            let mut rec = Record {
                row: row,
                seqname: record.seqname.clone(),
                source: record.source.clone(),
                feature_type: record.feature_type.clone(),
                start: record.start,
                end: record.end,
                score: record.score.clone(),
                strand: record.strand.clone(),
                frame: record.frame.clone(),
                attributes: LinkedHashMap::new(),
            };
            let mut parents = Vec::<String>::new();
            if let Some(parent_rows) = self.row2parents.get(&row) {
                for parent_row in parent_rows {
                    let parent = &self.rows[*parent_row];
                    if let Some(parent_id) = parent.attributes.get("ID") {
                        parents.push(parent_id.clone());
                    }
                }
            }
            for (k,v) in &record.attributes {
                if k == "Parent" {
                    rec.attributes.insert(k.clone(), String::from(parents.join(",")));
                }
                else {
                    rec.attributes.insert(k.clone(), v.clone());
                }
            }
            output.write_fmt(format_args!("{}\n", rec.to_gff()?))?;
        }
        Ok(())
    }
    
    pub fn to_bed(&self, 
        bed_file: &str, 
        exon_types: &[String],
        cds_types: &[String],
        transcript_types: &[String],
        gene_types: &[String])
        -> Result<()>
    {
        let exon_types = exon_types.iter().map(|t| String::from(t.as_ref())).collect::<HashSet<String>>();
        let mut cds_types = cds_types.iter().map(|t| String::from(t.as_ref())).collect::<HashSet<String>>();
        if cds_types.is_empty() { cds_types.insert(String::from("CDS")); }
        let transcript_types = transcript_types.iter().map(|t| String::from(t.as_ref())).collect::<HashSet<String>>();
        let gene_types = gene_types.iter().map(|t| String::from(t.as_ref())).collect::<HashSet<String>>();
        let mut seen_transcript = HashSet::<usize>::new();
        let mut transcript_names = HashSet::<String>::new();
        
        let unsorted_bed_file = format!("{}.unsorted.bed", bed_file);
        {   let mut bw = BufWriter::new(File::create(&unsorted_bed_file)?);
            let mut chrs = self.tree.keys().collect::<Vec<_>>();
            chrs.sort_by_key(|a| a.as_bytes());
            for chr in chrs {
                for node in self.tree[chr].find(0..std::u64::MAX) {
                    let gene_row = node.data();
                    let gene = &self.rows[*gene_row];
                    if gene_types.is_empty() || gene_types.contains(&gene.feature_type) {
                        if let Some(ref gene_children) = self.row2children.get(gene_row) {
                            for transcript_row in gene_children.iter() {
                                if !seen_transcript.contains(transcript_row) {
                                    let transcript = &self.rows[*transcript_row];
                                    if transcript_types.is_empty() || transcript_types.contains(&transcript.feature_type) {
                                        seen_transcript.insert(*transcript_row);
                                        
                                        // get the exon and CDS features
                                        let mut exons = Vec::<usize>::new();
                                        let mut cds_features = Vec::<usize>::new();
                                        let mut exon_starts = HashSet::<u64>::new();
                                        let mut exon_ends = HashSet::<u64>::new();
                                        if let Some(child_rows) = self.row2children.get(transcript_row) {
                                            for child_row in child_rows {
                                                let child = &self.rows[*child_row];
                                                if child.seqname == transcript.seqname {
                                                    if cds_types.contains(&child.feature_type) {
                                                        cds_features.push(*child_row);
                                                    }
                                                    else if exon_types.is_empty() || exon_types.contains(&child.feature_type) {
                                                        exon_starts.insert(self.rows[*child_row].start-1);
                                                        exon_ends.insert(self.rows[*child_row].end);
                                                        exons.push(*child_row);
                                                    }
                                                }
                                            }
                                        }
                                        if exons.is_empty() { exons.push(*transcript_row); }
                                        exons.sort_by_key(|a| self.rows[*a].start);
                                        cds_features.sort_by_key(|a| self.rows[*a].start);
                                        
                                        // get the cds ranges by unsplicing the CDS features
                                        let mut cdss = Vec::<Range<u64>>::new();
                                        for cds_feature_row in cds_features {
                                            let cds_feature = &self.rows[cds_feature_row];
                                            if !cdss.is_empty() && 
                                                exon_ends.contains(&cdss[cdss.len()-1].end) && 
                                                exon_starts.contains(&(cds_feature.start-1)) 
                                            {
                                                // extend current CDS if this is a splice
                                                if let Some(last) = cdss.last_mut() {
                                                    (*last).end = cds_feature.end;
                                                }
                                            }
                                            else {
                                                // otherwise add new CDS
                                                cdss.push((cds_feature.start-1)..cds_feature.end);
                                            }
                                        }
                                        cdss.sort_by_key(|a| a.start);
                                        
                                        // write the bed record
                                        let score = 0;
                                        let item_rgb = &[0,0,0].iter().map(|v| v.to_string()).join(",");
                                        let start = self.rows[exons[0]].start-1;
                                        let end = self.rows[exons[exons.len()-1]].end;
                                        // make a separate bed record for each cds
                                        if cdss.is_empty() { cdss.push(start..start) }
                                        for cds in cdss {
                                            // choose a unique transcript name
                                            let transcript_name = 
                                                transcript.attributes.get("transcript_name").or_else(||
                                                transcript.attributes.get("Name").or_else(||
                                                transcript.attributes.get("ID")));
                                            let mut transcript_name = match transcript_name {
                                                Some(t) => t.clone(),
                                                None => String::from(format!("{}:{}..{}:{}", 
                                                    transcript.seqname, transcript.start-1, transcript.end, transcript.strand)),
                                            };
                                            while transcript_names.contains(&transcript_name) {
                                                lazy_static! {
                                                    static ref TRANSCRIPT_RENAME: Regex = Regex::new(r"(?:\.([0-9]+))?$").unwrap();
                                                }
                                                transcript_name = String::from(TRANSCRIPT_RENAME.replace(&transcript_name.as_ref(), |caps: &Captures| {
                                                    format!(".{}", caps.get(1).map(|m| m.as_str().parse::<u64>().unwrap()).unwrap_or(0)+1)}));
                                            }
                                            transcript_names.insert(transcript_name.clone());
                                            
                                            let chr = self.vizchrmap.get(&transcript.seqname).unwrap_or(&transcript.seqname);
                                            // write the bed line
                                            let line = &[
                                                &chr,
                                                &start.to_string(),
                                                &end.to_string(),
                                                &transcript_name,
                                                &score.to_string(),
                                                &transcript.strand,
                                                &cds.start.to_string(),
                                                &cds.end.to_string(),
                                                item_rgb,
                                                &exons.len().to_string(),
                                                &exons.iter().map(|e| (self.rows[*e].end-self.rows[*e].start+1).to_string()).join(","),
                                                &exons.iter().map(|e| ((self.rows[*e].start-1)-start).to_string()).join(","),
                                            ].iter().join("\t");
                                            writeln!(bw, "{}", line)?;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // sort the bed file
        cmd!("sort","-k1,1","-k2,2n","-o",&bed_file,&unsorted_bed_file).env("LC_COLLATE","C").run()?;
        // remove the unsorted bed file
        std::fs::remove_file(&unsorted_bed_file)?;
        
        Ok(())
    }
    
    pub fn to_bigbed(&self, 
        file: &str, 
        exon_types: &[String],
        cds_types: &[String],
        transcript_types: &[String],
        gene_types: &[String],
        trackdb: &mut BufWriter<Box<Write>>)
        -> Result<()> 
    {
        // write the bed file
        let bed_file = format!("{}.bed", file);
        self.to_bed(&bed_file, exon_types, cds_types, transcript_types, gene_types)?;
        
        // write the genome file
        let genome_filename = format!("{}.genome", file);
        {   let mut genome_fh = BufWriter::new(File::create(&genome_filename)?);
            for (chr, length) in &self.refs {
                let chr = self.vizchrmap.get(chr).unwrap_or(chr);
                writeln!(genome_fh, "{}\t{}", chr, length)?;
            }
        }
        
        // call bedToBigBed
        cmd!("bedToBigBed","-type=bed12","-tab","-extraIndex=name", &bed_file, &genome_filename, &file).run()?;
        // remove the bed file
        //std::fs::remove_file(&bed_file)?;
        // remove the genome file
        std::fs::remove_file(&genome_filename)?;
        
        define_encode_set! {
            pub PATH_ENCODE_SET = [SIMPLE_ENCODE_SET] | {'+', '?', '&'}
        }
        let url = utf8_percent_encode(file, PATH_ENCODE_SET);
        let track_name = Path::new(file).file_stem().r()?.to_str().r()?;
        // write to the trackDb.txt file
        trackdb.write_fmt(format_args!("{}", unindent(&format!(r##"
        
            track {}
            type bigBed 12
            bigDataUrl {}
            shortLabel {}
            longLabel {}
            itemRgb On
            visibility hide
        "##, track_name, url, track_name, track_name))))?;
        trackdb.flush()?;
        
        Ok(())
    }
    
    pub fn read_fasta(fasta_file: &str) -> Result<LinkedHashMap<String,(String,String)>> {
        let mut fasta = LinkedHashMap::<String,(String,String)>::new();
        lazy_static! {
            static ref HEADER: Regex = Regex::new(r"^>(\S*)([^\n]*)").unwrap();
        }
        let f = File::open(&fasta_file)?;
        let mut file = BufReader::new(&f);
        let mut header: Option<String> = None;
        let mut attrs: Option<String> = None;
        let mut sequence: Option<String> = None;
        let mut line = String::new();
        while file.read_line(&mut line)? > 0 {
            let line = line.trim_right_matches('\n').trim_right_matches('\r');
            if let Some(cap) = HEADER.captures(&line) {
                if let Some(header) = header {
                    if let Some(attrs) = attrs {
                        if let Some(sequence) = sequence {
                            fasta.insert(header.clone(), (attrs.clone(), sequence.clone()));
                        }
                    }
                }
                header = Some(cap[1].to_string());
                attrs = Some(cap[2].to_string());
                sequence = Some(String::new());
            }
            else if let Some(ref mut sequence) = sequence {
                sequence.push_str(line.trim());
            }
        }
        if let Some(header) = header {
            if let Some(attrs) = attrs {
                if let Some(sequence) = sequence {
                    fasta.insert(header, (attrs, sequence));
                }
            }
        }
        Ok(fasta)
    }
    
    pub fn to_fasta(&self, 
        fasta_file: &str, 
        genome_file: &str,
        exon_types: &[String],
        transcript_types: &[String],
        gene_types: &[String])
        -> Result<()> 
    {
        let genome = IndexedAnnotation::read_fasta(genome_file)?;
        let exon_types = exon_types.iter().map(|t| String::from(t.as_ref())).collect::<HashSet<String>>();
        let transcript_types = transcript_types.iter().map(|t| String::from(t.as_ref())).collect::<HashSet<String>>();
        let gene_types = gene_types.iter().map(|t| String::from(t.as_ref())).collect::<HashSet<String>>();
        
        let mut bw = BufWriter::new(File::create(&fasta_file)?);
        let mut chrs = self.tree.keys().collect::<Vec<_>>();
        chrs.sort_by_key(|a| a.as_bytes());
        for chr in chrs {
            for node in self.tree[chr].find(0..std::u64::MAX) {
                let gene_row = node.data();
                let gene = &self.rows[*gene_row];
                if gene_types.is_empty() || gene_types.contains(&gene.feature_type) {
                    if let Some(ref gene_children) = self.row2children.get(gene_row) {
                        for transcript_row in gene_children.iter() {
                            let transcript = &self.rows[*transcript_row];
                            if transcript_types.is_empty() || transcript_types.contains(&transcript.feature_type) {
                                let mut exons = Vec::<usize>::new();
                                if let Some(child_rows) = self.row2children.get(transcript_row) {
                                    for child_row in child_rows {
                                        let child = &self.rows[*child_row];
                                        if child.seqname == transcript.seqname {
                                            if exon_types.is_empty() || exon_types.contains(&child.feature_type) {
                                                exons.push(*child_row);
                                            }
                                        }
                                    }
                                }
                                if exons.is_empty() { exons.push(*transcript_row); }
                                exons.sort_by_key(|a| self.rows[*a].start);
                                let mut transcript_seq = String::new();
                                for exon_row in &exons {
                                    let exon = &self.rows[*exon_row];
                                    if let Some(seq) = genome.get(&exon.seqname) {
                                        if exon.end as usize > seq.1.len() {
                                            return Err(format!("to_fasta: Range {}..{} of exon at row {} exceeds sequence for {} in file {}!", 
                                                exon.start-1, exon.end, exon_row, exon.seqname, genome_file).into());
                                        }
                                        let exon_seq = &seq.1[(exon.start-1) as usize..exon.end as usize];
                                        transcript_seq.push_str(&exon_seq);
                                    }
                                    else {
                                        return Err(format!("to_fasta: Could not find fasta sequence for {} in file {}!",
                                            exon.seqname, genome_file).into());
                                    }
                                }
                                if transcript.strand == "-" {
                                    let seq = transcript_seq.into_bytes();
                                    transcript_seq = String::from_utf8(dna::revcomp(&seq))?;
                                }
                                
                                let transcript_name: Option<String> = 
                                    transcript.attributes.get("transcript_id").or_else(||
                                    transcript.attributes.get("ID")).
                                    map(|a| a.to_owned());
                                if transcript_name.is_none() {
                                    return Err(format!("Could not get ID for transcript at row {}", transcript_row).into());
                                }
                                let transcript_name = transcript_name.r()?;
                                lazy_static! {
                                    static ref FASTA_FORMAT: Regex = Regex::new(r".{1,72}").unwrap();
                                }
                                transcript_seq = FASTA_FORMAT.replace_all(&transcript_seq,"$0\n").into_owned();
                                
                                define_encode_set! {
                                    pub GFF_ENCODE_SET = [SIMPLE_ENCODE_SET] | {'\t', '\r', '\n', ';', '%', '='}
                                }
                                let attrs = transcript.attributes.iter().map(|(k,v)| 
                                    format!("{}={}", 
                                        utf8_percent_encode(k, GFF_ENCODE_SET), 
                                        utf8_percent_encode(v, GFF_ENCODE_SET))).join("; ");
                                write!(bw, ">{} {}\n{}", transcript_name, attrs, transcript_seq)?;
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }
}
