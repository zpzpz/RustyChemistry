use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::helpers::*;

pub struct inputFile {

    pub basis: String,
    pub nuclei: Vec<Nucleus>,
    pub dim: usize,
    pub Nelec: usize,


}

pub fn testingInputFile()->inputFile{
    let basis = "STO3G".to_string();
    let nuclei = vec![
        Nucleus{
            atom_no: 1,
            x: 0_f64,
            y: 0_f64,
            z: 0_f64
        },
        Nucleus{
            atom_no: 0,
            x: 0_f64,
            y: 0_f64,
            z: 1.511_f64
        }
    ];
    let dim = 2;
    let Nelec = 2;

    return inputFile{basis, nuclei, dim, Nelec}

}

pub fn printInputFile(file: &inputFile){
    println!("Basis Set: {}", file.basis);
    println!("Number of basis functions: {}", file.dim);
    println!("Number of electrons: {}", file.Nelec);
    println!("Number of atoms: {}", file.nuclei.len());
    for nuc in &file.nuclei{
        println!("{} atom at position {}, {}, {}", atomicSymbol(nuc.atom_no), nuc.x, nuc.y, nuc.z);
    }
}

pub fn ReadInputFile(FileName: &String)->inputFile{
    let file = File::open(FileName).expect("input file not found");;
    let reader = BufReader::new(file);
    
    let mut lines = vec![];
    for line in reader.lines(){
        let line = line.expect("error reading lines");
        lines.push(line.clone());
    }
    let basis = &lines[0];
    let dim: usize = lines[1].parse().expect("error reading input file");
    let Nelec: usize = lines[2].parse().expect("error reading input file");
    let mut nuclei = vec![];
    for i in 3..lines.len(){
        let values: Vec<&str> = lines[i].split_whitespace().collect();
        let atom_no: usize = values[0].parse().expect("error reading input file");
        let x: f64 = values[1].parse().expect("error reading input file");
        let y: f64 = values[2].parse().expect("error reading input file");
        let z: f64 = values[3].parse().expect("error reading input file");
        let nuc = Nucleus{atom_no, x, y, z};
        nuclei.push(nuc);
    } 
    return inputFile{basis: basis.to_string(), nuclei, dim, Nelec};
}
