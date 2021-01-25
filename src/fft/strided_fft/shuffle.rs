use std::{fmt::Debug, ptr::swap};

use crate::fft::fft::{
    copy_mut_alias
};

use super::super::mem_utils::prefetch_slice_element;
use crate::fft::multicore::Worker;
// assuming matrix is row-major ordered
pub fn transpose<T: Sized>(matrix: &mut [T], size: usize, stretch: usize) {
    assert_eq!(matrix.len(), size * size * stretch);
    match stretch {
        1 => transpose_square(matrix, size),
        2 => transpose_n_by_2n(matrix, size),
        _ => unimplemented!(
            "Only stretch sizes 1 and 2 are supported, received stretch {}",
            stretch
        ),
    }
}

pub fn transpose_square<T: Sized>(matrix: &mut [T], size: usize) {
    const PREFETCH_STRIDE: usize = 4;
    // let step_by: usize = 64 / std::mem::size_of::<T>();
    let step_by = 2;
    debug_assert_eq!(matrix.len(), size * size);
    assert_eq!(size & 1, 0);

    // Iterate over upper-left triangle, working in 2x2 blocks
    // Stretches of two are useful because they span a 64B cache line when T is 32
    // bytes.

    for row in (0..size).step_by(step_by) {
        let i = row * size + row;
        matrix.swap(i + 1, i + size);
        for col in (row..size).step_by(step_by).skip(1) {
            let i = row * size + col;
            let j = col * size + row;
            if PREFETCH_STRIDE > 0 && col + PREFETCH_STRIDE * 2 < size {
                prefetch_slice_element(matrix, i + PREFETCH_STRIDE * 2);
                prefetch_slice_element(matrix, i + PREFETCH_STRIDE * 2 + size);
                prefetch_slice_element(matrix, j + PREFETCH_STRIDE * 2 * size);
                prefetch_slice_element(matrix, j + PREFETCH_STRIDE * 2 * size + size);
            }
            matrix.swap(i, j);
            matrix.swap(i + 1, j + size);
            matrix.swap(i + size, j + 1);
            matrix.swap(i + size + 1, j + size + 1);
        }
    }
}

pub fn transpose_n_by_2n<T: Sized>(matrix: &mut [T], size: usize) {
    const PREFETCH_STRIDE: usize = 8;
    debug_assert_eq!(matrix.len(), 2 * size * size);

    // Iterate over upper-left triangle, working in 1x2 blocks
    for row in 0..size {
        for col in (row..size).skip(1) {
            let i = (row * size + col) * 2;
            let j = (col * size + row) * 2;
            if PREFETCH_STRIDE > 0 && col + PREFETCH_STRIDE < size {
                prefetch_slice_element(matrix, i + PREFETCH_STRIDE * 2 * size);
            }
            matrix.swap(i, j);
            matrix.swap(i + 1, j + 1);
        }
    }
}

pub fn transpose_square_with_square_tiles<T: Sized, const N: usize>(matrix: &mut [T], size: usize) {
    debug_assert_eq!(size % N, 0);
    let num_tiles = size / N;

    for column_square_idx in 0..num_tiles {
        let base = column_square_idx * N * size;
        // swap diagonal
        for row in 1..N {
            for column in row..N {
                matrix.swap(base + row * size + column, base + column * size + row);
            }
        }
        // prefetch
        if column_square_idx != 0 {
            let row_square_idx = column_square_idx - 1;
            let src_start = (row_square_idx + 1) * N * size + column_square_idx * N;
            let dst_start = column_square_idx * N * size + (row_square_idx + 1) * N;
            for row in 0..N {
                prefetch_slice_element(&*matrix, src_start + row * size);
                prefetch_slice_element(&*matrix, dst_start + row * size);
            }
        }

        // swap the rest
        for row_square_idx in (1..column_square_idx).rev() {
            // prefetch for next iteration
            {
                let src_start = (row_square_idx + 1) * N * size + column_square_idx * N;
                let dst_start = column_square_idx * N * size + (row_square_idx + 1) * N;
                for row in 0..N {
                    prefetch_slice_element(&*matrix, src_start + row * size);
                    prefetch_slice_element(&*matrix, dst_start + row * size);
                }
            }
            let src_start = row_square_idx * N * size + column_square_idx * N;
            let dst_start = column_square_idx * N * size + row_square_idx * N;
            for row in 0..N {
                for column in 0..N {
                    // walk over row of src, fill the columns of dst
                    matrix.swap(
                        src_start + row * size + column,
                        dst_start + column * size + row,
                    );
                }
            }
        }

        // final swap
        let row_square_idx = 0;
        let src_start = row_square_idx * N * size + column_square_idx * N;
        let dst_start = column_square_idx * N * size + row_square_idx * N;
        for row in 0..N {
            for column in 0..N {
                // walk over row of src, fill the columns of dst
                matrix.swap(
                    src_start + row * size + column,
                    dst_start + column * size + row,
                );
            }
        }
    }
}

#[inline]
pub fn transpose_square_segment_ttt<T: Sized + Copy, const N: usize>(
    matrix: &mut [T],
    src_start: usize,
    dst_start: usize,
    row_width: usize,
) {
    use std::mem::MaybeUninit;
    let mut src_array: [MaybeUninit<&mut [T; N]>; N] = MaybeUninit::uninit_array();
    let mut dst_array: [MaybeUninit<&mut [T; N]>; N] = MaybeUninit::uninit_array();
    let mut view = matrix;
    let mut src_start = src_start;
    for i in 0..N {
        let (new_view, row_start) =
            super::super::fft::get_mut_assuming_not_overlapping(view, src_start);
        super::super::mem_utils::prefetch_element(row_start);
        unsafe {
            src_array[i] = MaybeUninit::new(&mut *(row_start as *mut T as *mut [T; N]));
        }
        view = new_view;
        src_start += row_width;
    }

    let mut dst_start = dst_start;
    for i in 0..N {
        let (new_view, row_start) =
            super::super::fft::get_mut_assuming_not_overlapping(view, dst_start);
        super::super::mem_utils::prefetch_element(row_start);
        unsafe {
            dst_array[i] = MaybeUninit::new(&mut *(row_start as *mut T as *mut [T; N]));
        }
        view = new_view;
        dst_start += row_width;
    }

    let src_array = unsafe { MaybeUninit::array_assume_init(src_array) };
    let dst_array = unsafe { MaybeUninit::array_assume_init(dst_array) };

    for row in 0..N {
        for column in 0..N {
            // walk over row of src, fill the columns of dst
            std::mem::swap(&mut src_array[row][column], &mut dst_array[column][row]);
        }
    }
}

#[inline]
pub fn transpose_square_segment<T: Sized + Copy, const N: usize>(
    matrix: &mut [T],
    src_start: usize,
    dst_start: usize,
    row_width: usize,
) {
    use std::mem::MaybeUninit;
    let mut src_array: [MaybeUninit<&mut [T; N]>; N] = MaybeUninit::uninit_array();
    let mut dst_array: [MaybeUninit<&mut [T; N]>; N] = MaybeUninit::uninit_array();
    let mut view = matrix;
    let mut src_start = src_start;
    for i in 0..N {
        let (new_view, row_start) =
            super::super::fft::get_mut_assuming_not_overlapping(view, src_start);
        super::super::mem_utils::prefetch_element(row_start);
        unsafe {
            src_array[i] = MaybeUninit::new(&mut *(row_start as *mut T as *mut [T; N]));
        }
        view = new_view;
        src_start += row_width;
    }

    let mut dst_start = dst_start;
    for i in 0..N {
        let (new_view, row_start) =
            super::super::fft::get_mut_assuming_not_overlapping(view, dst_start);
        super::super::mem_utils::prefetch_element(row_start);
        unsafe {
            dst_array[i] = MaybeUninit::new(&mut *(row_start as *mut T as *mut [T; N]));
        }
        view = new_view;
        dst_start += row_width;
    }

    let src_array = unsafe { MaybeUninit::array_assume_init(src_array) };
    let dst_array = unsafe { MaybeUninit::array_assume_init(dst_array) };

    for row in 0..N {
        for column in 0..N {
            // walk over row of src, fill the columns of dst
            std::mem::swap(&mut src_array[row][column], &mut dst_array[column][row]);
        }
    }
}

pub fn transpose_square_with_chunks<T: Sized + Copy, const N: usize>(
    matrix: &mut [T],
    size: usize,
) {
    debug_assert_eq!(matrix.len(), size * size);
    use super::super::fft::*;
    let view_as_chunked_rows = subview_mut::<_, N>(matrix);
    assert_eq!(size % N, 0);

    // we work on "squares" of size N by N, where N is chosen in a way such that
    // N field elements fit into the cache line of the processor

    let num_subsquares = size / N;

    for column_square_idx in 0..num_subsquares {
        // transpose diagonal element
        let (src, dst) = copy_mut_alias(view_as_chunked_rows);
        let src: [&mut [T; N]; N] =
            get_as_square(src, column_square_idx, column_square_idx, num_subsquares);
        let dst: [&mut [T; N]; N] =
            get_as_square(dst, column_square_idx, column_square_idx, num_subsquares);
        for row in 0..N {
            for column in row..N {
                std::mem::swap(&mut src[row][column], &mut dst[column][row]);
            }
        }
        for row_square_idx in (0..column_square_idx).rev() {
            let (src, dst) = copy_mut_alias(view_as_chunked_rows);
            let src: [&mut [T; N]; N] =
                get_as_square(src, row_square_idx, column_square_idx, num_subsquares);
            let dst: [&mut [T; N]; N] =
                get_as_square(dst, column_square_idx, row_square_idx, num_subsquares);
            for row in 0..N {
                for column in 0..N {
                    // walk over row of src, fill the columns of dst
                    std::mem::swap(&mut src[row][column], &mut dst[column][row]);
                }
            }
        }
    }
}

fn print_matrix<T: Sized + Copy + Debug>(matrix: &mut [T], dim: usize) {

    for chunk in matrix.chunks_exact(dim) {
        for val in chunk {
            print!("{:?}\t", val);
        }
        println!("")
    }
}

fn dialate_encode(val: u64) -> u64 {
    let mut val = val;
    val &= 0x0000ffff;
    val = (val ^ (val << 8)) & 0x00ff00ff;
    val = (val ^ (val << 4)) & 0x0f0f0f0f;
    val = (val ^ (val << 2)) & 0x33333333;
    val = (val ^ (val << 1)) & 0x55555555;

    val
}

fn dialate_decode(val: u64) -> u64 {
    let mut val = val;
    val &= 0x55555555;
    val = (val ^ (val >> 1)) & 0x33333333;
    val = (val ^ (val >> 2)) & 0x0f0f0f0f;
    val = (val ^ (val >> 4)) & 0x00ff00ff;
    val = (val ^ (val >> 8)) & 0x0000ffff;

    val
}

fn zorder_encode(x: u64, y: u64) -> u64 {
    (dialate_encode(y) << 1) | dialate_encode(x)
}

fn zorder_decode(k: u64) -> (u64, u64) {
    let x = dialate_decode(k);
    let y = dialate_decode(k >> 1);

    (x, y)
}

fn construct_matrix_with_morton_order(dim: usize) -> Vec<usize> {
    let mut matrix = vec![0; dim * dim];
    for i in 0..dim {
        for j in 0..dim {
            let z_idx = zorder_encode(j as u64, i as u64) as usize;
            matrix[z_idx] = j + i * dim;
        }
    }

    matrix
}

pub fn natural_to_morton<T: Sized + Copy>(src: &[T], dim: usize) -> Vec<T> {
    let mut result = src.to_vec();
    for row in 0..dim {
        for col in 0..dim {
            let k = zorder_encode(col as u64, row as u64) as usize;
            result[k] = src[row * dim + col];
        }
    }

    result
}

pub fn morton_to_natural<T: Sized + Copy>(src: &[T], dim: usize) -> Vec<T> {
    let mut result = src.to_vec();
    for i in 0..src.len() {
        let (row, col) = zorder_decode(i as u64);
        let id = col as usize * dim + row as usize;
        result[id] = src[i];
    }

    result
}

fn is_power_of_4(n: u64) -> bool {
    let mut n = n;

    if n == 0 {
        return false;
    }

    while n != 1 {
        if n % 4 != 0 {
            return false;
        }
        n = n / 4;
    }
    return true;
}

pub fn recursive_morton<T: Sized + Copy + Debug>(matrix: &mut [T], dim: usize) {
    assert!(is_power_of_4(matrix.len() as u64));
    assert_eq!(matrix.len(), dim * dim);
    let size = dim * dim;

    if dim == 2 {
        matrix.swap(1, 2);
        return;
    }

    // assuma matrix is in morton order

    // we split matirx into four sub squares
    // by doing so, dimension will be 1/2 of original dimension
    let new_dim = dim / 2;

    // number of elements in matrix will be new_dim*new_dim
    let new_size = new_dim * new_dim;
    
    // a00
    recursive_morton(&mut matrix[0 * new_size..1 * new_size], new_dim);

    // a01
    recursive_morton(&mut matrix[1 * new_size..2 * new_size], new_dim);

    // a10
    recursive_morton(&mut matrix[2 * new_size..3 * new_size], new_dim);

    // a11
    recursive_morton(&mut matrix[3 * new_size..4 * new_size], new_dim);

    // swap 4x4
    let number_of_elems_in_block = 4;

    let number_of_rows_for_single_sub_matrix = dim / 4;

    let number_of_element_in_single_row = number_of_rows_for_single_sub_matrix * dim;
    

    // find start of first row and iterate over then swap first blocks;
    let number_of_blocks_in_row = number_of_element_in_single_row / number_of_elems_in_block;
    // swap 2  and 3 row by 4x4 blocks
    let start_of_first_row = 1 * number_of_element_in_single_row;
    let start_of_second_row = 2 * number_of_element_in_single_row;

    let (matrix, matrix_alias) = copy_mut_alias(matrix);

    for block_idx in 0..number_of_blocks_in_row {
        let first_idx = start_of_first_row + block_idx * number_of_elems_in_block;
        let second_idx = start_of_second_row + block_idx * number_of_elems_in_block;

        matrix[first_idx..first_idx + number_of_elems_in_block]
            .swap_with_slice(&mut matrix_alias[second_idx..second_idx + number_of_elems_in_block]);

        // iterate through each block and swap them
        // for idx in 0..number_of_elems_in_block{
        //     prefetch_slice_element(&*matrix, first_idx + idx);
        //     prefetch_slice_element(&*matrix, second_idx + idx);
        // }
        // for idx in 0..number_of_elems_in_block{
        //     matrix.swap(
        //         first_idx + idx,
        //         second_idx + idx
        //     );
        // }
    }
}


pub fn recursive_morton_with_multicore<T: Sized + Copy + Debug + Send>(
    matrix: &mut [T],
    dim: usize,
    worker: &Worker,
) {
    assert!(is_power_of_4(matrix.len() as u64));
    assert_eq!(matrix.len(), dim * dim);
    let size = dim * dim;

    if dim == 2 {
        matrix.swap(1, 2);
        return;
    }

    // assuma matrix is in morton order

    // we split matirx into four sub squares
    // by doing so, dimension will be 1/2 of original dimension
    let new_dim = dim / 2;

    // number of elements in matrix will be new_dim*new_dim
    let new_size = new_dim * new_dim;

    recursive_morton(&mut matrix[0 * new_size..1 * new_size], new_dim);

    recursive_morton(&mut matrix[1 * new_size..2 * new_size], new_dim);

    recursive_morton(&mut matrix[2 * new_size..3 * new_size], new_dim);

    recursive_morton(&mut matrix[3 * new_size..4 * new_size], new_dim);

    // swap 4x4
    let number_of_elems_in_block = 1;

    let number_of_rows_for_single_sub_matrix = dim / 4;

    let number_of_element_in_single_row = number_of_rows_for_single_sub_matrix * dim;

    // find start of first row and iterate over then swap first blocks;
    let number_of_blocks_in_row = number_of_element_in_single_row / number_of_elems_in_block;

    // swap 2  and 3 row by 4x4 blocks
    let start_of_first_row = 1 * number_of_element_in_single_row;
    let start_of_second_row = 2 * number_of_element_in_single_row;

    let (matrix, matrix_alias) = copy_mut_alias(matrix);

    worker.scope(number_of_blocks_in_row, |scope, chunk_size| {        
        for (chunk_idx, (first, second)) in matrix
            [start_of_first_row..start_of_first_row + number_of_blocks_in_row * number_of_elems_in_block]
            .chunks_mut(chunk_size)
            .zip(
                matrix_alias
                    [start_of_second_row..start_of_second_row + number_of_blocks_in_row * number_of_elems_in_block]
                    .chunks_mut(chunk_size),
            )
            .enumerate()
        {

            scope.spawn(move |_| {
                // first[..].swap_with_slice(
                //     &mut second[..],
                // );
                // let offset = chunk_idx*chunk_size;
                let num_elements_per_slice = first.len();
                for idx in 0..num_elements_per_slice{
                    prefetch_slice_element(first, idx);
                    prefetch_slice_element(second, idx);
                }
                for idx in 0..num_elements_per_slice{
                    let tmp = first[idx];
                    first[idx] = second[idx];
                    second[idx] = tmp;
                }
            });
        }
    });
}

fn transpose_naive<T: Sized + Copy + Debug, const DIM: usize>(matrix: &mut [T]) {
    for row in 0..DIM {
        for col in row..DIM {
            let tmp = matrix[row + col * DIM];
            matrix[row + col * DIM] = matrix[col + row * DIM];
            matrix[col + row * DIM] = tmp;
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        fft::fft::{get_as_square, subview},
        optimized_fields::naive_f125::Fr,
    };
    use crate::fft::multicore::Worker;
    use ff::{Field, PrimeField};

    #[test]
    fn test_matrix_tranpose_with_outplace() {
        let worker = Worker::new();
        const DIM: usize = 8;

        let mut matrix = vec![];
        for i in 0..(DIM * DIM) {
            // let el = Fr::from_str(&i.to_string()).unwrap();
            // matrix.push(el);
            matrix.push(i);
        }

        let mut dst_expected = matrix.clone();
        transpose_naive::<_, DIM>(&mut dst_expected);

        let mut matrix_in_morton_order = natural_to_morton(&matrix, DIM as usize);

        recursive_morton_with_multicore(&mut matrix_in_morton_order, DIM as usize, &worker);        

        let result_in_natural = morton_to_natural(&mut matrix_in_morton_order, DIM as usize);
        assert_eq!(dst_expected, result_in_natural);
    }

    #[test]
    fn test_zorder() {
        let x = 2;
        let y = 3;
        let k = zorder_encode(x, y);
        let (_x, _y) = zorder_decode(k);
        assert_eq!(k, 14);
        assert_eq!(x, _x);
        assert_eq!(y, _y);
    }
    
    #[test]
    fn test_transpose_single() {
        use crate::optimized_fields::naive_f125::Fr;
        const N: usize = 4;

        let size = N * 1;

        let mut matrix = vec![];
        for i in 0..(size * size) {
            let el = Fr::from_str(&i.to_string()).unwrap();
            matrix.push(el);
        }

        println!("Original: {:?}", &matrix);

        let mut naive_transposed = matrix.clone();
        for row in 0..size {
            for column in row..size {
                let tmp = naive_transposed[column * size + row];
                naive_transposed[column * size + row] = naive_transposed[row * size + column];
                naive_transposed[row * size + column] = tmp;
            }
        }

        transpose_square_with_chunks::<_, N>(&mut matrix, size);

        println!("Transposed: {:?}", &matrix);

        assert_eq!(naive_transposed, matrix);
    }

    #[test]
    fn test_transpose_multiple() {
        use crate::optimized_fields::naive_f125::Fr;
        const N: usize = 4;

        let size = N * 2;

        let mut matrix = vec![];
        for i in 0..(size * size) {
            let el = Fr::from_str(&i.to_string()).unwrap();
            matrix.push(el);
        }

        println!("Original: {:?}", &matrix);

        let mut naive_transposed = matrix.clone();
        for row in 0..size {
            for column in row..size {
                let tmp = naive_transposed[column * size + row];
                naive_transposed[column * size + row] = naive_transposed[row * size + column];
                naive_transposed[row * size + column] = tmp;
            }
        }

        transpose_square_with_chunks::<_, N>(&mut matrix, size);

        println!("Transposed: {:?}", &matrix);

        for row in 0..size {
            for column in 0..size {
                let idx = row * size + column;
                if naive_transposed[idx] != matrix[idx] {
                    panic!(
                        "Invalid element for row {}, column {}: expected {}, got {}",
                        row, column, naive_transposed[idx], matrix[idx]
                    );
                }
            }
        }
    }

    #[test]
    fn test_diagonal() {
        let n = 4;
        let mut m = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

        for i in 0..n {
            for j in 0..n {
                m[i * n + j] = m[j * n + i];
            }
        }

        m.iter().for_each(|e| println!("{}", e));
    }
}
