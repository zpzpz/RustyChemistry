#![allow(warnings)]
/* RustyChemistry v0.0.1 
 * Licensed under GNU Public License V3
 *
 * Current Basis Sets Implemented: STO-3G (first row only)
 *
 * Current Methods Implemented: Hartree-Fock (first row only)
 *
 * Current bugs: energy is not calculated properly, and SCF convergence does not work properly.
 *
 * TODO: 
 *  -implement 2p orbitals
 *  -finish STO-3G down to neon
 *  -implement json for basis sets
 *  -implement Pople basis sets
 *  -implement geometry optimization
 *  -implement output as XYZ files
 *
 *
 */
const version_string: &str = "0.0.1";
extern crate openblas_src;
use integrate::newton_cotes::trapezoidal_rule;
use num_traits::pow::Pow;
use libm::erf;
use libm::sqrt;
use num::integer::div_floor;
use ndarray::arr2;
use ndarray::{Array2, Array1};
use ndarray_linalg::Eig;
use fast_inv_sqrt::InvSqrt64;
use std::env;

pub mod helpers;
pub mod sto3g;
pub mod hartreefock;
pub mod inputfile;

use crate::helpers::*;
use crate::sto3g::*;
use crate::hartreefock::*;
use crate::inputfile::*;

pub const E: f64 = 2.71828182845904523536028747135266250_f64;
pub const pi: f64 = 3.14159265358979323846264338327950288_f64;

fn main() {
    println!("RustyChemistry v{} ", version_string);
    //read args
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        println!("usage: {} inputfile.txt", args[0]);
        std::process::exit(1);
    }
    let mut verbose = false;
    let mut inputfile = String::new();
    let mut i = 1; // Skip program name at index 0
    while i < args.len() {
        match args[i].as_str() {
            "-v" | "--verbose" => {
                verbose = true;
                i += 1;
            }
            arg => {
                inputfile = arg.to_string();
                i += 1;
            }
        }
    }

    if verbose{
        println!("{}",inputfile);
    }

    //read input file
    let input = ReadInputFile(&inputfile);
    println!("\ninput file successfully read as: ");
    printInputFile(&input);

    let Nelec = input.Nelec.clone(); 
    let dim = input.dim.clone();
    
    //initialize orbital and nuclei vectors
    let mut orbitals = vec![];
    let mut nuclei = vec![];
    let mut basis = vec![];
    for i in 0..input.nuclei.len(){
        nuclei.push(&input.nuclei[i]);
        basis.push(build_atomic_basis(&input.nuclei[i]));
    }
    //borrowing fuckery
    for i in 0..basis.len(){
        for j in 0..basis[i].len(){
            orbitals.push(&basis[i][j]);
        }
    }

    if verbose{
        for i in 0..orbitals.len(){
            print_basis_function_sto3g(&orbitals[i]);
        }
    }
    //overlap matrix
    let mut S = vec![];

    //kinetic matrix
    let mut T = vec![];
    
    //potential matrix
    let mut V = vec![];
    
    //initialize 1 electron matrices
    for i in 0..dim{
        S.push(vec![]);
        T.push(vec![]);
        V.push(vec![]);
        for j in 0..dim{
            S[i].push(0_f64);
            T[i].push(0_f64);
            V[i].push(0_f64);
        }
    }

    //two electron vector (yoshimine indexing)
    let mut twoe: Vec<f64> = vec![];
    //initialize
    let num_twoe = (dim*(dim+1)*(dim*dim+dim+2))/8;
    for i in 0..num_twoe{
        twoe.insert(i,0_f64);
    }
    
    //list orbitals for convenience of indexing
    let NUC = nuclear_energy(&nuclei);
    if verbose{
        println!("\nnuclear repulsion energy: {} \n", NUC);
    }
    //calculate 1 electron integrals
    if verbose{
        println!("calculating 1 electron integrals");
    }
    for i in 0..dim{
        for j in 0..dim{
            S[i][j] = calc_overlap_integral_sto3g(orbitals[i], orbitals[j]);
            T[i][j] = calc_kinetic_energy_integral_sto3g(orbitals[i], orbitals[j]);
            for k in 0..nuclei.len(){//num of nuclei
                V[i][j] += calc_attraction_integral_sto3g(orbitals[i], orbitals[j], nuclei[k]);
            }
        }
    }
    
    let H = matrix_add(&T, &V);
    if verbose{
        println!("calculating 2 electron integrals");
    }
    let mut used_indexs: Vec<usize> = vec![];
    let mut cur_key = 0;
    for i in 0..dim{
        for j in 0..dim{
            for k in 0..dim{
                for l in 0..dim{
                    cur_key = eint(i,j,k,l);
                    if !(used_indexs.contains(&cur_key)){
                        twoe[cur_key] = calc_electron_repulsion_integral_sto3g(orbitals[i], orbitals[j], orbitals[k], orbitals[l]);
                        if verbose{
                            print!("i: {}, j: {}, k: {}, l: {} ", i,j,k,l);
                            print!("int: {} \n", twoe[cur_key]);
                        }
                        used_indexs.push(cur_key);
                    }
                }
            }
        }
    }
    if verbose{
        println!("\ntesting self overlap");
        for i in 0..dim{
            let s_rr = calc_overlap_integral_sto3g(&orbitals[i], &orbitals[i]);
            println!("S_rr = {}", s_rr);  // Should be very close to 1.0
        }
    
        println!("\none electron matrices");
        println!("overlap: ");
        matrix_print(&S);
        println!("kinetic: ");
        matrix_print(&T);
        println!("potential: ");
        matrix_print(&V);
        println!("hamiltonian: ");
        matrix_print(&H);
        println!();
        println!("\ntwo electron integrals: ");
        for i in 0..twoe.len(){
            print!("yoshi: {} ", i);
            print!("val: {} ", twoe[i]);
            println!();
        }
    }
    let X = diagonalize_matrix(&S);//diagnalzie overlap
    if verbose{
        println!("\ncheck overlap transform");
        matrix_print(&X);
        println!();
        matrix_print(&S);
        println!();
        let identity = matrix_mult(&matrix_mult(&transpose(&X), &S), &X);
        matrix_print(&identity);
        println!();
    }
    //initialize guess for density matrix (empty)
    let mut P: Vec<Vec<f64>> = [].to_vec();
    let mut Pold: Vec<Vec<f64>> = [].to_vec();
    for i in 0..dim{
        P.push([].to_vec());
        Pold.push([].to_vec());
        for j in 0..dim{
            P[i].push(0_f64);
            Pold[i].push(0_f64);
        }
    }
    let mut P_vec = vec![P.clone(),Pold.clone()];
    let mut F = make_fock_matrix(&H, &P_vec[0], dim, &twoe);
    if verbose{
        println!("\ninitial fock guess");
        matrix_print(&F);
    }
    let mut Fprime = matrix_mult(&transpose(&X), &matrix_mult(&F, &X));

    let (mut eigenvalues, mut diag_Fprime) = calc_eigenvec(&Fprime);
    let mut C = matrix_mult(&X, &diag_Fprime);
    if verbose{
        println!("\ncoefficient matrix");
        matrix_print(&C);
    }
    P_vec = make_density(&C, &P, dim, Nelec);
    P = P_vec[0].clone();
    Pold = P_vec[1].clone();
    if verbose{
        println!("\n initial density matrix");
        matrix_print(&P);
    }
    println!("\nstarting self-consistent field iteration: \n");
    let mut count = 0;
    let mut DELTA = 1_f64;
    while DELTA > 0.000001 && count < 25{
        count += 1;

        F = make_fock_matrix(&H,&P,dim,&twoe);
        Fprime = matrix_mult(&transpose(&X), &matrix_mult(&F, &X));
        (eigenvalues, diag_Fprime) = calc_eigenvec(&Fprime);
        C = matrix_mult(&X, &diag_Fprime);
        P_vec = make_density(&C, &P, dim, Nelec);
        P = P_vec[0].clone();
        Pold = P_vec[1].clone();

        DELTA = delta_p(&P, &Pold);
        println!("{}: energy: {1:.10}, RMSD: {2:.10}", count, energy(&P, &H, &F, dim)+NUC, DELTA);
    }
    println!("\nSCF successfully converged!");
    println!("final energy: {} hartrees in {} cycles", energy(&P, &H, &F, dim)+NUC, count);
    
}
