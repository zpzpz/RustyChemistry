//sto3g.rs
//implements 1 and 2 electron integrals in sto3g basis set
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::HashMap;
use crate::helpers::*;
use crate::pi;
use libm::sqrt;

pub struct BasisFunction_STO3G {
    pub a: Gauss_Prim,
    pub b: Gauss_Prim,
    pub c: Gauss_Prim,
    pub ang_x: usize,
    pub ang_y: usize,
    pub ang_z: usize
}

fn get_prim_from_basis_sto3g(principal: usize, momentum: usize, primitive_index: usize, nuc: &Nucleus) -> Gauss_Prim {
    let basis_text = include_str!("sto-3g.txt");
    let element = &atomicSymbol(nuc.atom_no);

    let mut prims = Vec::new();
    let mut in_element = false;
    let mut shell_count = 0;
    let mut lines = basis_text.lines().peekable();

    while let Some(line) = lines.next() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('!') {
            continue;
        }

        // Check if this is the element header
        if line.split_whitespace().next() == Some(element) {
            in_element = true;
            shell_count = 0;
            continue;
        }

        // End of this element's basis
        if in_element && line == "****" {
            break;
        }

        // Parse orbital headers (S or SP)
        if in_element {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 3 && (parts[0] == "S" || parts[0] == "SP") {
                shell_count += 1;
                let num_prims: usize = parts[1].parse().unwrap_or(0);
                let is_sp = parts[0] == "SP";

                // If this is the shell we want, collect its primitives
                if shell_count == principal {
                    for _ in 0..num_prims {
                        while let Some(prim_line) = lines.next() {
                            let prim_line = prim_line.trim();
                            if prim_line.is_empty() || prim_line.starts_with('!') {
                                continue;
                            }

                            let prim_parts: Vec<&str> = prim_line.split_whitespace().collect();
                            if prim_parts.len() >= 2 {
                                let exp = prim_parts[0].replace('D', "E").parse::<f64>().unwrap_or(0.0);
                                let coeff = if !is_sp || momentum == 0 {
                                    prim_parts[1].replace('D', "E").parse::<f64>().unwrap_or(0.0)
                                } else {
                                    prim_parts[2].replace('D', "E").parse::<f64>().unwrap_or(0.0)
                                };
                                prims.push((exp, coeff));
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    if primitive_index < prims.len() {
        let (exponent, coefficient) = prims[primitive_index];
        return Gauss_Prim { c: coefficient, a: exponent, x: nuc.x, y: nuc.y, z: nuc.z };
    }

    Gauss_Prim { c: 0.0, a: 0.0, x: nuc.x, y: nuc.y, z: nuc.z }
}

/*fn get_prim_from_basis_sto3g(principal: usize, momentum: usize, magnetic: usize, primitive_index: usize, nuc: &Nucleus)->Gauss_Prim{
    let basis = load_sto_3g();
    let coefficient = basis[nuc.atom_no][principal][momentum][0][primitive_index];
    let exponent = basis[nuc.atom_no][principal][momentum][1][primitive_index];
    return Gauss_Prim{c: coefficient, a: exponent, x: nuc.x, y: nuc.y, z: nuc.z};
}*/

fn Build_BasisFunction_STO3G (nuc: &Nucleus, principal: usize, momentum: usize, Ang_x: usize, Ang_y: usize, Ang_z: usize)->BasisFunction_STO3G {
    let A = get_prim_from_basis_sto3g(principal, momentum, 0, &nuc);
    let B = get_prim_from_basis_sto3g(principal, momentum, 1, &nuc);
    let C = get_prim_from_basis_sto3g(principal, momentum, 2, &nuc);

    return BasisFunction_STO3G{a: A, b: B, c: C, ang_x: Ang_z, ang_y: Ang_y, ang_z: Ang_z};
}

pub fn build_atomic_basis (nuc: &Nucleus)->Vec<BasisFunction_STO3G>{
    let mut output = vec![];
    if nuc.atom_no < 2 { //hydrogen and helium get one basis function
        output.push(Build_BasisFunction_STO3G(nuc, 1, 0, 0, 0, 0));
    } else if nuc.atom_no >= 2 && nuc.atom_no < 4 { //lithium and berylium get two
        output.push(Build_BasisFunction_STO3G(nuc, 1, 0, 0, 0, 0));
        output.push(Build_BasisFunction_STO3G(nuc, 2, 0, 0, 0, 0));
    } else if nuc.atom_no >=4 && nuc.atom_no < 9 { //boron->neon get five
        output.push(Build_BasisFunction_STO3G(nuc, 1, 0, 0, 0, 0));
        output.push(Build_BasisFunction_STO3G(nuc, 2, 0, 0, 0, 0));
        output.push(Build_BasisFunction_STO3G(nuc, 2, 1, 1, 0, 0)); //P_x
        output.push(Build_BasisFunction_STO3G(nuc, 2, 1, 0, 1, 0)); //P_y
        output.push(Build_BasisFunction_STO3G(nuc, 2, 1, 0, 0, 1)); //P_z
    }
    output
}



//overlap integral
fn calc_overlap_integral_prim(r: &Gauss_Prim, s: &Gauss_Prim)->f64{
    let overlap_prim = Gauss_Prod(r, s);
    return normalize_gaussian(r.a) * normalize_gaussian(s.a) * overlap_prim.c * ((pi)/(overlap_prim.a)).powf(3_f64/2_f64);
}

//sum of the integrals is the integral of the sum
pub fn calc_overlap_integral_sto3g(r: &BasisFunction_STO3G, s: &BasisFunction_STO3G)->f64{

    //list prims so we can iterate over them nicely: todo: rewrite
    let primsR = [
        &r.a,
        &r.b,
        &r.c
    ];
    let primsS = [
        &s.a,
        &s.b,
        &s.c
    ];

    let mut sum = 0.0_f64;
    for i in 0..3{
        for j in 0..3{
            sum += primsR[i].c * primsS[j].c * calc_overlap_integral_prim(primsR[i], primsS[j]);
        }
    }

    return sum;
}

//s orbitals only implemented right now
fn calc_kinetic_energy_integral_prim(r: &Gauss_Prim, s: &Gauss_Prim)->f64{
    let overlap_prim = Gauss_Prod(r, s);
    let AB2 = (r.x-s.x).powf(2_f64) + (r.y-s.y).powf(2_f64) + (r.z-s.z).powf(2_f64);

    let pre = normalize_gaussian(r.a) * normalize_gaussian(s.a) * (pi / overlap_prim.a).powf(1.5_f64) * overlap_prim.c;
    return r.a * s.a / overlap_prim.a * (3_f64 - 2_f64 * r.a * s.a / overlap_prim.a * AB2) * pre;
}
pub fn calc_kinetic_energy_integral_sto3g(r: &BasisFunction_STO3G, s: &BasisFunction_STO3G)->f64{
    let primsR = [
        &r.a,
        &r.b,
        &r.c
    ];
    let primsS = [
        &s.a,
        &s.b,
        &s.c
    ];

    let mut sum = 0.0_f64;
    for i in 0..3{
        for j in 0..3{
            sum += primsR[i].c * primsS[j].c * calc_kinetic_energy_integral_prim(primsR[i], primsS[j]); 
        }
    }

    return sum;
}

//for now only 1s orbitals implemented until i do the calculus xD
fn calc_attraction_integral_prim(r: &Gauss_Prim, s: &Gauss_Prim, nuc: &Nucleus)->f64{
    let overlap_prim = Gauss_Prod(r, s);
    let PC2 = (overlap_prim.x-nuc.x).powf(2_f64) + (overlap_prim.y-nuc.y).powf(2_f64) + (overlap_prim.z-nuc.z).powf(2_f64);
    let Z: f64 = (nuc.atom_no+1) as f64;
    return (normalize_gaussian(r.a) * normalize_gaussian(s.a) * (-2_f64 * pi * Z / overlap_prim.a) * boys_zero(overlap_prim.a * PC2) * overlap_prim.c);

}

//only 1s integrals implemented right now
pub fn calc_attraction_integral_sto3g(r: &BasisFunction_STO3G, s: &BasisFunction_STO3G, nuc: &Nucleus)->f64{
    let primsR = [
        &r.a,
        &r.b,
        &r.c
    ];
    let primsS = [
        &s.a,
        &s.b,
        &s.c
    ];

    let mut sum = 0.0_f64;
    for i in 0..3{
        for j in 0..3{
            sum += primsR[i].c * primsS[j].c * calc_attraction_integral_prim(primsR[i], primsS[j], &nuc);
        }
    }
    return sum;
}

//electron repulsion integral - two electron
//only s integrals implemented right now
fn calc_electron_repulsion_integral_prim(A: &Gauss_Prim, B: &Gauss_Prim, C: &Gauss_Prim, D: &Gauss_Prim)->f64{
    let P = Gauss_Prod(A,B);
    let Q = Gauss_Prod(C,D);
    let PQ2 = (P.x-Q.x).powf(2_f64) + (P.y-Q.y).powf(2_f64) + (P.z-Q.z).powf(2_f64);
    let zeta = (P.a*Q.a)/(P.a+Q.a);
    let pre = (normalize_gaussian(A.a) * normalize_gaussian(B.a) * normalize_gaussian(C.a) * normalize_gaussian(D.a)
        * 2_f64 * pi.powf(2.5_f64) / (P.a * Q.a * sqrt(P.a + Q.a))
        * P.c * Q.c);

    return pre * boys_zero(zeta * PQ2);
}

pub fn calc_electron_repulsion_integral_sto3g(r: &BasisFunction_STO3G, s: &BasisFunction_STO3G, u: &BasisFunction_STO3G, v: &BasisFunction_STO3G)->f64{
    let mut sum = 0.0_f64;
    
    let primsR = [&r.a, &r.b, &r.c];
    let primsS = [&s.a, &s.b, &s.c];
    let primsU = [&u.a, &u.b, &u.c];
    let primsV = [&v.a, &v.b, &v.c];

    for i in 0..3{
        for j in 0..3{
            for k in 0..3{
                for l in 0..3{
                        sum += primsR[i].c * primsS[j].c * primsU[k].c * primsV[l].c * calc_electron_repulsion_integral_prim(primsR[i], primsS[j], primsU[k], primsV[l]);
                }
            }
        }
    }
    
    return sum;
}

//pretty print macro
pub fn print_basis_function_sto3g(A: &BasisFunction_STO3G){
    println!("STO3G basis function made of 3 gaussian primitives");
    print_gauss_prim(&A.a);
    print_gauss_prim(&A.b);
    print_gauss_prim(&A.c);
}

