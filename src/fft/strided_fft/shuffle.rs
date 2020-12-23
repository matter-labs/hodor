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