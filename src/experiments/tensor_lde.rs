use ff::*;

fn query_tensor_decomposition_element_matrix_over_identyty_matrix<F: PrimeField>(
    submatrix: &(Vec<F>, (usize, usize)),
    idx: (usize, usize)
) -> F {
    // this is a tensor decomposition query that takes matrix
    // in a form of
    // |a|0|0|0|                |1|0|0|0|
    // |0|a|0|0|   =  a \cross  |0|1|0|0|
    // |0|0|a|0|                |0|0|1|0|
    // |0|0|0|a|                |0|0|0|1|


    if idx.0 / (submatrix.1).0 != idx.1 / (submatrix.1).1 {
        return F::zero();
    }

    let sub_row = idx.0 % (submatrix.1).0;
    let sub_column = idx.1 % (submatrix.1).1;

    let row_len = (submatrix.1).1;

    (submatrix.0)[row_len * sub_row + sub_column]
}

fn query_tensor_decomposition_element_matrix_over_diagonal_matrix<F: PrimeField>(
    submatrix: &(Vec<F>, (usize, usize)),
    diagonal_matrix: &(Vec<F>, usize),
    idx: (usize, usize)
) -> F {
    // this is a tensor decomposition query that takes matrix
    // in a form of
    // |a|0|0|0|                |b|0|0|0|
    // |0|a|0|0|   =  a \cross  |0|c|0|0|
    // |0|0|a|0|                |0|0|d|0|
    // |0|0|0|a|                |0|0|0|e|


    if idx.0 / (submatrix.1).0 != idx.1 / (submatrix.1).1 {
        return F::zero();
    }

    let diagonal_matrix_element_index = idx.0 / (submatrix.1).0;

    let sub_row = idx.0 % (submatrix.1).0;
    let sub_column = idx.1 % (submatrix.1).1;

    let row_len = (submatrix.1).1;

    let mut subelement = (submatrix.0)[row_len * sub_row + sub_column];

    subelement.mul_assign(&(diagonal_matrix.0)[diagonal_matrix_element_index]);

    subelement
}

fn query_tensor_decomposition_element_vector_over_vector<F: PrimeField>(
    subvector_1: &(Vec<F>, usize),
    subvector_2: &(Vec<F>, usize),
    idx: usize
) -> F {
    // this is a tensor decomposition query that takes matrix
    // in a form of
    // |a|                    
    // |b|   =  [a, b] \cross [c', d']
    // |c|                    
    // |d|                   

    let element_index_0 = idx % (subvector_1.1);
    let element_index_1 = idx / (subvector_1.1);

    assert!(element_index_1 < subvector_2.1);

    let mut subelement = (subvector_1.0)[element_index_0];

    subelement.mul_assign(&(subvector_2.0)[element_index_1]);

    subelement
}

fn decompose_lde_generator_for_vector_over_vector<F: PrimeField>(
    lde_factor: usize,
    domain_size: usize,
    // natural_index: usize,
    decomposition_info: (usize, usize),
    omega: &F,
    coset_generator: &F,
) -> ((F, F, usize), (F, F, usize)) {

    // one can think about the LDE of some power series of values
    // 1, alpha, alpha^2
    // as of evaluations of the polynomial that interpolates this series
    // at the points {1, v, v^2, ..., v^{l-1} } * {1, omega, omega^2, ...},
    // where for interpolation f(1) = 1, f(omega) = alpha, ... (natural enumeration)
    // with additional requirement that v^l = omega
    // so when one wants to query the LDE using LDEs of tensor decomposition it's necessary
    // to send the some part of the LDE basis element v^k * omega^n into
    // (v^k', omega^n', v^k'', omega^n'') (NOTE: omegas in every case should still form the basis!)
    // different parts of tensor decomposition, so when those parts are combined together (by multiplication)
    // they would form the element that would be LDE(v^k * omega^n)


    // let coset_index = natural_index % lde_factor;
    // let domain_index = natural_index / domain_size; 

    // for multiplicative domain we still have that if 
    // omega is a generator for size 2k,
    // then omega^2 is a generator for size k

    let subvector_1_domain_size = decomposition_info.0;
    let subvector_2_domain_size = decomposition_info.1;

    assert!(subvector_1_domain_size * subvector_2_domain_size == domain_size);

    // v^{lde factor} = omega

    // easiest options is that 
    // values of one of LDEs of the subvector should capture 
    // v^k only,
    // and values on another should capture excessive powers of omega

    let coset_generator_1 = *coset_generator;
    let omega_1 = omega.pow([(domain_size / subvector_1_domain_size) as u64]);
    let lde_factor_1 = lde_factor;

    let coset_generator_2 = *omega;
    let omega_2 = omega.pow([(domain_size / subvector_2_domain_size) as u64]);
    let lde_factor_2 = domain_size / subvector_2_domain_size; // == subvector_1_domain_size


    ((coset_generator_1, omega_1, lde_factor_1), (coset_generator_2, omega_2, lde_factor_2))
}

#[test]
fn test_query_matrix_by_identity() {
    use super::Fr;
    let zero = Fr::zero();
    let one = Fr::one();
    let two = Fr::from_str("2").unwrap();
    let three = Fr::from_str("3").unwrap();

    let original_matrix = vec![
        zero, one, zero, zero,
        two, three, zero, zero,
        zero, zero, zero, one,
        zero, zero, two, three
    ];

    let submatrix = vec![
        zero, one, 
        two, three
    ];

    let submatrix_dim = (2, 2);
    let matrix_dim = (4, 4);

    let submatrix_info = (submatrix, submatrix_dim);

    for query_index in 0..original_matrix.len() {
        let query_row = query_index / matrix_dim.1;
        let query_column = query_index % matrix_dim.1;

        let expected = original_matrix[query_index];
        let value = query_tensor_decomposition_element_matrix_over_identyty_matrix(&submatrix_info, (query_row, query_column));
        assert!(value == expected);
    }
}

#[test]
fn test_query_matrix_by_diagonal_identity() {
    use super::Fr;
    let zero = Fr::zero();
    let one = Fr::one();
    let two = Fr::from_str("2").unwrap();
    let three = Fr::from_str("3").unwrap();

    let original_matrix = vec![
        zero, one, zero, zero,
        two, three, zero, zero,
        zero, zero, zero, one,
        zero, zero, two, three
    ];

    let submatrix = vec![
        zero, one, 
        two, three
    ];

    let submatrix_dim = (2, 2);
    let matrix_dim = (4, 4);

    let submatrix_info = (submatrix, submatrix_dim);

    let diagonal_identity_matrix = vec![one, one];
    let diagonal_identity_matrix_info = (diagonal_identity_matrix, 2);

    for query_index in 0..original_matrix.len() {
        let query_row = query_index / matrix_dim.1;
        let query_column = query_index % matrix_dim.1;

        let expected = original_matrix[query_index];
        let value = query_tensor_decomposition_element_matrix_over_diagonal_matrix(
            &submatrix_info, 
            &diagonal_identity_matrix_info,
            (query_row, query_column));

        assert!(value == expected);
    }
}

#[test]
fn test_query_matrix_by_diagonal_matrix() {
    use super::Fr;
    let zero = Fr::zero();
    let one = Fr::one();
    let two = Fr::from_str("2").unwrap();
    let three = Fr::from_str("3").unwrap();

    let four = Fr::from_str("4").unwrap();
    let six = Fr::from_str("6").unwrap();


    let original_matrix = vec![
        zero, one, zero, zero,
        two, three, zero, zero,
        zero, zero, zero, two,
        zero, zero, four, six
    ];

    let submatrix = vec![
        zero, one, 
        two, three
    ];

    let submatrix_dim = (2, 2);
    let matrix_dim = (4, 4);

    let submatrix_info = (submatrix, submatrix_dim);

    let diagonal_matrix = vec![one, two];
    let diagonal_matrix_info = (diagonal_matrix, 2);

    for query_index in 0..original_matrix.len() {
        let query_row = query_index / matrix_dim.1;
        let query_column = query_index % matrix_dim.1;

        let expected = original_matrix[query_index];
        let value = query_tensor_decomposition_element_matrix_over_diagonal_matrix(
            &submatrix_info, 
            &diagonal_matrix_info,
            (query_row, query_column));

        assert!(value == expected);
    }
}

#[test]
fn test_query_vector_of_powers() {
    use super::Fr;

    let alpha = Fr::from_str("123").unwrap();
    let num_powers = 16;

    let mut original_powers = vec![];
    for i in 0..num_powers {
        let tmp = alpha.pow([i as u64]);
        original_powers.push(tmp);
    }

    let subvector_1_dim = 2;
    let mut subvector_1 = vec![];
    for i in 0..subvector_1_dim {
        let tmp = alpha.pow([i as u64]);
        subvector_1.push(tmp);
    }

    let subvector_2_dim = num_powers / subvector_1_dim;
    let mut subvector_2 = vec![];
    for i in 0..subvector_2_dim {
        let tmp = alpha.pow([(i*subvector_1_dim) as u64]);
        subvector_2.push(tmp);
    }

    let subvector_1_info = (subvector_1, subvector_1_dim);
    let subvector_2_info = (subvector_2, subvector_2_dim);

    for query_index in 0..original_powers.len() {
        let expected = original_powers[query_index];
        let value = query_tensor_decomposition_element_vector_over_vector(
            &subvector_1_info, 
            &subvector_2_info,
            query_index
        );

        assert!(value == expected);
    }
}

#[test]
fn test_tensor_lde() {
    use super::Fr;

    let alpha = Fr::from_str("123").unwrap();
    let num_powers = 16;

    let lde_factor = 16;

    let mut original_powers = vec![];
    for i in 0..num_powers {
        let tmp = alpha.pow([i as u64]);
        original_powers.push(tmp);
    }

    let subvector_1_dim = 2;
    let mut subvector_1 = vec![];
    for i in 0..subvector_1_dim {
        let tmp = alpha.pow([i as u64]);
        subvector_1.push(tmp);
    }

    let subvector_2_dim = num_powers / subvector_1_dim;
    let mut subvector_2 = vec![];
    for i in 0..subvector_2_dim {
        let tmp = alpha.pow([(i*subvector_1_dim) as u64]);
        subvector_2.push(tmp);
    }

    let subvector_1_info = (subvector_1, subvector_1_dim);
    let subvector_2_info = (subvector_2, subvector_2_dim);

    use crate::domains::*;

    let lde_domain = Domain::<Fr>::new_for_size((lde_factor * num_powers) as u64).unwrap();

    let coset_generator = lde_domain.generator;

    let main_domain = Domain::<Fr>::new_for_size(num_powers as u64).unwrap();

    let omega = main_domain.generator;

    let info = decompose_lde_generator_for_vector_over_vector(
        lde_factor,
        num_powers,
        (subvector_1_info.1, subvector_2_info.1),
        &omega,
        &coset_generator
    );

    use crate::fft::multicore::*;
    use crate::polynomials::*;

    let worker = Worker::new();

    let poly_1 = Polynomial::<Fr, Values>::from_values(subvector_1_info.0.clone()).unwrap();
    let poly_2 = Polynomial::<Fr, Values>::from_values(subvector_2_info.0.clone()).unwrap();

    let poly_1 = poly_1.ifft(&worker);
    let poly_2 = poly_2.ifft(&worker);

    let mut poly_1_lde = Polynomial::<Fr, Values>::new_for_size((info.0).2).unwrap();

    let mut poly_2_lde = Polynomial::<Fr, Values>::new_for_size((info.1).2).unwrap();

    assert!((info.0).1 == poly_1.omega);
    assert!((info.1).1 == poly_2.omega);

    let coset_size = (info.0).2;

    for coset_idx in 0..coset_size {
        let poly_1_coset_values = poly_1.clone().coset_fft_for_generator(&worker, (info.0).0);
        assert!(poly_1_coset_values.size() == poly_1.size());
        for j in 0..poly_1.size() {
            poly_1_lde.as_mut()[j * coset_size + coset_idx] = poly_1_coset_values.as_ref()[j];
        }
    }

    let coset_size = (info.1).2;

    for coset_idx in 0..coset_size {
        let poly_2_coset_values = poly_2.clone().coset_fft_for_generator(&worker, (info.1).0);
        for j in 0..poly_2.size() {
            poly_2_lde.as_mut()[j * coset_size + coset_idx] = poly_2_coset_values.as_ref()[j];
        }
    }


}










