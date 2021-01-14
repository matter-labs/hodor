use super::super::mem_utils::prefetch_slice_element;

// assuming matrix is row-major ordered
pub fn transpose<T: Sized>(matrix: &mut [T], size: usize, stretch: usize) {
    assert_eq!(matrix.len(), size * size * stretch);
    match stretch {
        1 => transpose_square(matrix, size),
        2 => transpose_n_by_2n(matrix, size),
        _ => unimplemented!("Only stretch sizes 1 and 2 are supported, received stretch {}", stretch),
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
                        dst_start + column * size + row
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
                    dst_start + column * size + row
                );
            }
        }
    }
}

#[inline]
pub fn transpose_square_segment_ttt<T: Sized + Copy, const N: usize>(matrix: &mut [T], src_start: usize, dst_start: usize, row_width: usize) {
    use std::mem::MaybeUninit;
    let mut src_array: [MaybeUninit<&mut [T; N]>; N] = MaybeUninit::uninit_array();
    let mut dst_array: [MaybeUninit<&mut [T; N]>; N] = MaybeUninit::uninit_array();
    let mut view = matrix;
    let mut src_start = src_start;
    for i in 0..N {
        let (new_view, row_start) = super::super::fft::get_mut_assuming_not_overlapping(view, src_start);
        super::super::mem_utils::prefetch_element(row_start);
        unsafe {
            src_array[i] = MaybeUninit::new(&mut *(row_start as *mut T as *mut [T; N]));
        }
        view = new_view;
        src_start += row_width;
    }

    let mut dst_start = dst_start;
    for i in 0..N {
        let (new_view, row_start) = super::super::fft::get_mut_assuming_not_overlapping(view, dst_start);
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
pub fn transpose_square_segment<T: Sized + Copy, const N: usize>(matrix: &mut [T], src_start: usize, dst_start: usize, row_width: usize) {
    use std::mem::MaybeUninit;
    let mut src_array: [MaybeUninit<&mut [T; N]>; N] = MaybeUninit::uninit_array();
    let mut dst_array: [MaybeUninit<&mut [T; N]>; N] = MaybeUninit::uninit_array();
    let mut view = matrix;
    let mut src_start = src_start;
    for i in 0..N {
        let (new_view, row_start) = super::super::fft::get_mut_assuming_not_overlapping(view, src_start);
        super::super::mem_utils::prefetch_element(row_start);
        unsafe {
            src_array[i] = MaybeUninit::new(&mut *(row_start as *mut T as *mut [T; N]));
        }
        view = new_view;
        src_start += row_width;
    }

    let mut dst_start = dst_start;
    for i in 0..N {
        let (new_view, row_start) = super::super::fft::get_mut_assuming_not_overlapping(view, dst_start);
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

pub fn transpose_square_with_chunks<T: Sized + Copy, const N: usize>(matrix: &mut [T], size: usize) {
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
        let src: [&mut [T; N]; N] = get_as_square(src, column_square_idx, column_square_idx, num_subsquares);
        let dst: [&mut [T; N]; N] = get_as_square(dst, column_square_idx, column_square_idx, num_subsquares);
        for row in 0..N {
            for column in row..N {
                std::mem::swap(&mut src[row][column], &mut dst[column][row]);
            }
        }
        for row_square_idx in (0..column_square_idx).rev() {
            let (src, dst) = copy_mut_alias(view_as_chunked_rows);
            let src: [&mut [T; N]; N] = get_as_square(src, row_square_idx, column_square_idx, num_subsquares);
            let dst: [&mut [T; N]; N] = get_as_square(dst, column_square_idx, row_square_idx, num_subsquares);
            for row in 0..N {
                for column in 0..N {
                    // walk over row of src, fill the columns of dst 
                    std::mem::swap(&mut src[row][column], &mut dst[column][row]);
                }
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::ff::PrimeField;

    #[test]
    fn test_transpose_single() {
        use crate::optimized_fields::naive_f125::Fr;
        const N: usize = 4;

        let size = N * 1;

        let mut matrix = vec![];
        for i in 0..(size*size) {
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
        for i in 0..(size*size) {
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
                    panic!("Invalid element for row {}, column {}: expected {}, got {}", row, column, naive_transposed[idx], matrix[idx]);
                }
            }
        }
    }
}