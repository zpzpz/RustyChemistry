//hartreefock.rs
//contains hartree fock helper and matrix builder functions

use ndarray::arr2;
use ndarray::{Array2, Array1};

use crate::helpers::*;


//makes fock matrix
//H: hamiltonian matrix
//P: density matrixx
//dim: number of basis functions
//tei: two electron integral vector using yoshimine indexing
pub fn make_fock_matrix(H: &Vec<Vec<f64>>, P: &Vec<Vec<f64>>, dim: usize, tei: &Vec<f64>)->Vec<Vec<f64>>{
    let mut output = vec![];
    //initialize fock matrix as hamiltonian
    for i in 0..dim{
        output.push(vec![]);
        for j in 0..dim{
            output[i].push(H[i][j]);
        }
    }

    for i in 0..dim{
        for j in 0..dim{
            for k in 0..dim{
                for l in 0..dim{
                    output[k][l] = output[k][l] + P[i][j]*(tei[eint(i,j,k,l)]-0.5*tei[eint(i,k,j,l)]);
                }
            }
        }
    }

    return output;
}

//returns the transformed fock matrix, Fprime
//X: transform matrix
//F: fock matrix
pub fn transform_fock(X: &Vec<Vec<f64>>, F: &Vec<Vec<f64>>)->Vec<Vec<f64>>{
    //cast to arrays
    let rows_X = X.len();
    let cols_X = X[0].len();

    let flat_x: Vec<f64> = X.iter().flat_map(|row| row.iter().cloned()).collect();
    let x_matrix = Array2::from_shape_vec((rows_X, cols_X), flat_x).unwrap();

    let rows_F = F.len();
    let cols_F = F[0].len();

    let flat_f: Vec<f64> = F.iter().flat_map(|row| row.iter().cloned()).collect();
    let f_matrix = Array2::from_shape_vec((rows_F, cols_F), flat_f).unwrap();

    let dot = f_matrix.dot(&x_matrix);
    let prod = x_matrix.t().dot(&dot);

    prod
        .outer_iter()
        .map(|row| row.to_vec())
        .collect()
}

//returns the updated hartree-fock density matrix along with the old density matrix
//C: coefficient matrix
//D: old density matrix
//dim: number of basis functions
//Nelec: number of electrons
pub fn make_density(C: &Vec<Vec<f64>>, D: &Vec<Vec<f64>>, dim: usize, Nelec: usize)->Vec<Vec<Vec<f64>>>{
    let mut Dold = vec![];
    for mu in 0..dim{
        Dold.push([].to_vec());
        for nu in 0..dim{
            Dold[mu].push(0_f64);
        }
    }
    let mut output = vec![];
    for mu in 0..dim{
        output.push([].to_vec());
        for nu in 0..dim{
            output[mu].push(D[mu][nu]);
        }
    }

    for mu in 0..dim{
        for nu in 0..dim{
            Dold[mu][nu] = D[mu][nu];
            output[mu][nu] = 0_f64;
            for m in 0..Nelec/2{
                output[mu][nu] = output[mu][nu] + 2_f64*C[mu][m]*C[nu][m];
            }
        }
    }
    return vec![output,Dold];
}

//calculates RMSD change in density matrix to determine convergence
//D: new density matrix
//Dold: old density matrix
pub fn delta_p(D: &Vec<Vec<f64>>, Dold: &Vec<Vec<f64>>)->f64{
    let mut output = 0_f64;
    let dim = D.len();
    for i in 0..dim{
        for j in 0..dim{
            output = output + ((D[i][j]-Dold[i][j]).powf(2_f64));
        }
    }
    output
}

//calculates hartree fock energy
//D: density matrix
//H: hamiltonian
//F: fock matrix
//dim: number of basis functions
pub fn energy(D: &Vec<Vec<f64>>, H: &Vec<Vec<f64>>, F: &Vec<Vec<f64>>, dim: usize)->f64{
    let mut output = 0_f64;
    for mu in 0..dim{
        for nu in 0..dim{
            output += 0.5 * D[mu][nu] * (H[mu][nu] + F[mu][nu]);
        }
    }
    output
}
