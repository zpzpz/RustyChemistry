//helpers.rs
//contains various helper functions
use crate::pi;
use ndarray::{Array2, Array1};
use ndarray_linalg::Eig;
use libm::erf;
use libm::sqrt;

pub struct Gauss_Prim {
   pub c: f64,
   pub a: f64,
   pub x: f64,
   pub y: f64,
   pub z: f64
}

pub struct Nucleus {
    pub atom_no: usize,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

//am i stupid?
pub fn atomicSymbol(n: usize)->String{
    if n == 0 { return "H".to_string();};
    if n == 1 { return "He".to_string();};
    if n == 2 { return "Li".to_string();};
    if n == 3 { return "Be".to_string();};
    if n == 4 { return "B".to_string();};
    if n == 5 { return "C".to_string();};
    if n == 6 { return "N".to_string();};
    if n == 7 { return "O".to_string();};
    if n == 8 { return "F".to_string();};
    if n == 9 { return "Ne".to_string();};
    return "null".to_string();
}

//returns the nuclear repulsion energy
pub fn nuclear_energy(nucs: &Vec<&Nucleus>)->f64{
    let mut output = 0_f64;

    for i in 0..nucs.len(){
        for j in (i+1)..nucs.len(){
            let dist = sqrt((nucs[i].x-nucs[j].x).powf(2_f64)+(nucs[i].y-nucs[j].y).powf(2_f64)+(nucs[i].z-nucs[j].z).powf(2_f64));
            output += ((nucs[i].atom_no+1) as f64 * (nucs[j].atom_no+1) as f64)/dist;
        }
    }

    return output;
}

//gaussian product theorem: the product of two gaussians is itself a gaussian
//returns the gaussian product of two gaussians A and B
pub fn Gauss_Prod(A: &Gauss_Prim, B: &Gauss_Prim)->Gauss_Prim{
    let AB2 = (A.x-B.x).powf(2_f64) + (A.y-B.y).powf(2_f64) + (A.z-B.z).powf(2_f64);
    let prod_c = (-((A.a * B.a) / (A.a + B.a)) * AB2).exp();
    let prod_a = A.a + B.a;
    let prod_x = (A.a * A.x + B.a * B.x)/(A.a + B.a);
    let prod_y = (A.a * A.y + B.a * B.y)/(A.a + B.a);
    let prod_z = (A.a * A.z + B.a * B.z)/(A.a + B.a);

    return Gauss_Prim{c: prod_c, a: prod_a, x: prod_x, y: prod_y, z: prod_z}
}

//pretty print macro
pub fn print_gauss_prim(A: &Gauss_Prim){
    print!("coefficient: {} ", A.c);
    print!("exp: {} ", A.a);
    print!("x: {} ", A.x);
    print!("y: {} ", A.y);
    print!("z: {} ", A.z);
    println!();
}

//zeroth order boys function
pub fn boys_zero(t: f64)->f64{
    if (t < 1e-7_f64){
        return 1.0_f64;
    }
    return 0.5 * sqrt(pi/t) * erf(sqrt(t));
}

//yoshimine compound index, thanks Adam Baskerville
pub fn eint(a: usize, b: usize, c: usize, d: usize)->usize{
    let mut ab = 0_usize;
    let mut cd = 0_usize;
        let ab = if a > b {
        a * (a + 1) / 2 + b
    } else {
        b * (b + 1) / 2 + a
    };

    let cd = if c > d {
        c * (c + 1) / 2 + d
    } else {
        d * (d + 1) / 2 + c
    };

    if ab > cd {
        ab * (ab + 1) / 2 + cd
    } else {
        cd * (cd + 1) / 2 + ab
    }
}


//crummy vector matrix implementation :: todo :: get away from this crap
pub fn transpose(matrix: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let rows = matrix.len();
    let cols = matrix[0].len();

    let mut transposed = vec![vec![0.0; rows]; cols];

    for i in 0..rows {
        for j in 0..cols {
            transposed[j][i] = matrix[i][j];
        }
    }

    transposed
}

//crummy vector matrix implementation
pub fn matrix_mult(A: &Vec<Vec<f64>>, B: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let rows_A = A.len();
    let cols_A = A[0].len();
    let cols_B = B[0].len();
    if cols_A != cols_B{
        println!("failed matrix multiplication!");
        return vec![vec![0_f64,0_f64]];
    }
    let mut result = vec![vec![0.0; cols_B]; rows_A];

    for i in 0..rows_A {
        for j in 0..cols_B {
            let mut sum = 0.0;
            for k in 0..cols_A {
                sum += A[i][k] * B[k][j];
            }
            result[i][j] = sum;
        }
    }

    result
}

//crummy vector matrix implementation
pub fn matrix_add(A: &Vec<Vec<f64>>, B: &Vec<Vec<f64>>)->Vec<Vec<f64>>{
    let mut output: Vec<Vec<f64>> = vec![];
    if A.len() != B.len(){
        return vec![vec![0_f64],vec![0_f64]];
    }
    for i in 0..A.len(){
        if A[i].len() != B[i].len(){
            return vec![vec![0_f64],vec![0_f64]];
        }
        output.insert(i,vec![]); //add a row
        for j in 0..A[i].len(){
            output[i].insert(j,A[i][j]+B[i][j]);
        }
    }
    return output;
}

//crummy vector matrix implementation
pub fn matrix_print(A: &Vec<Vec<f64>>){
    for row in A {
        for value in row {
            print!("{0:.8} ", value);
        }
        println!();
    }
}

//returns the 1s normalization constant for a gaussian primitive
pub fn normalize_gaussian(exponent: f64) -> f64 {
    // Normalization constant for 1s Gaussian: (2α/π)^(3/4)
    let norm = (2.0 * exponent / std::f64::consts::PI).powf(0.75);
    norm
}

//returns a sorted list of eigenvalues and a sorted matrix of eigenvectors
pub fn calc_eigenvec(A: &Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let rows = A.len();
    let cols = A[0].len();
    
    let flat: Vec<f64> = A.iter().flat_map(|row| row.iter().cloned()).collect();
    let matrix = Array2::from_shape_vec((rows, cols), flat).unwrap();
    
    let (eigenvalues, eigenvectors) = matrix.eig().unwrap();

    // Extract real parts
    let mut eigenvalues_real: Vec<f64> = eigenvalues.iter().map(|x| x.re).collect();
    let eigenvectors_real_raw = eigenvectors.mapv(|x| x.re);
    
    // Create pairs of (eigenvalue, index) and sort by eigenvalue
    let mut indexed_evals: Vec<(f64, usize)> = eigenvalues_real
        .iter()
        .enumerate()
        .map(|(i, &val)| (val, i))
        .collect();

    indexed_evals.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    // Reorder eigenvalues and eigenvectors by sorted indices
    let sorted_evals: Vec<f64> = indexed_evals.iter().map(|(val, _)| *val).collect();
    
    let sorted_evecs: Vec<Vec<f64>> = indexed_evals
        .iter()
        .map(|(_, idx)| {
            eigenvectors_real_raw
                .column(*idx)
                .to_vec()
        })
        .collect();
        
    // Transpose to get (N_basis × N_basis) matrix where each column is an eigenvector
    let n = sorted_evecs.len();
    let mut transposed = vec![vec![0.0; n]; n];
    for (i, evec) in sorted_evecs.iter().enumerate() {
        for (j, &val) in evec.iter().enumerate() {
            transposed[j][i] = val;
        }
    }

    (sorted_evals, transposed)
}

//returns overlap transformation matrix i.e.
//the matrix X s.t. X^t S X = 1
pub fn diagonalize_matrix(A: &Vec<Vec<f64>>)->Vec<Vec<f64>>{

    let rows = A.len();
    let cols = A[0].len();

    let flat: Vec<f64> = A.iter().flat_map(|row| row.iter().cloned()).collect();
    let matrix = Array2::from_shape_vec((rows, cols), flat).unwrap(); //cast to an ndarray

    let (eigenvalues, eigenvectors) = matrix.eig().unwrap();

    let invsqrt_eigenvals: Array1<f64> = eigenvalues.mapv(|x| 1.0 / x.re.sqrt());

    let mut diag_matrix = Array2::zeros((rows, cols));
    for i in 0..rows.min(cols) {
        diag_matrix[[i, i]] = invsqrt_eigenvals[i];
    }

    let eigenvectors_real = eigenvectors.mapv(|x| x.re);
    //let output = eigenvectors_real.dot(&diag_matrix).dot(&eigenvectors_real.t());
    let output = eigenvectors_real.dot(&diag_matrix.dot(&eigenvectors_real.t()));

    output
        .outer_iter()
        .map(|row| row.to_vec())
        .collect()
}

fn hermite_polynomial(n: usize, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => 2.0 * x,
        _ => {
            let mut h_prev = 1.0;      // H_0
            let mut h_curr = 2.0 * x;  // H_1
            
            for k in 1..n {
                let h_next = 2.0 * x * h_curr - 2.0 * (k as f64) * h_prev;
                h_prev = h_curr;
                h_curr = h_next;
            }
            h_curr
        }
    }
}

