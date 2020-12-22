use ff::PrimeField;
use super::multicore::*;
use crate::utils::*;

pub(crate) fn best_fft<F: PrimeField>(a: &mut [F], worker: &Worker, omega: &F, log_n: u32, use_cpus_hint: Option<usize>)
{
    let log_cpus = if let Some(hint) = use_cpus_hint {
        assert!(hint <= worker.cpus);
        log2_floor(hint)
    } else {
        worker.log_num_cpus()
    };

    if log_cpus == 0 || log_n <= log_cpus {
        serial_fft(a, omega, log_n);
    } else {
        parallel_fft(a, worker, omega, log_n, log_cpus);
    }
}

pub(crate) fn serial_fft<F: PrimeField>(a: &mut [F], omega: &F, log_n: u32)
{
    #[inline(always)]
    fn bitreverse(mut n: u32, l: u32) -> u32 {
        let mut r = 0;
        for _ in 0..l {
            r = (r << 1) | (n & 1);
            n >>= 1;
        }
        r
    }

    let n = a.len() as u32;
    assert_eq!(n, 1 << log_n);

    for k in 0..n {
        let rk = bitreverse(k, log_n);
        if k < rk {
            a.swap(rk as usize, k as usize);
        }
    }

    let mut m = 1;
    for _ in 0..log_n {
        let w_m = omega.pow(&[(n / (2*m)) as u64]);

        let mut k = 0;
        while k < n {
            let mut w = F::one();
            for j in 0..m {
                
                let mut t = a[(k+j+m) as usize];
                t.mul_assign(&w);
                let mut tmp = a[(k+j) as usize];
                tmp.sub_assign(&t);
                a[(k+j+m) as usize] = tmp;
                a[(k+j) as usize].add_assign(&t);
                w.mul_assign(&w_m);
            }

            k += 2*m;
        }
        
        m *= 2;
    }
}

pub(crate) fn parallel_fft<F: PrimeField>(
    a: &mut [F],
    worker: &Worker,
    omega: &F,
    log_n: u32,
    log_cpus: u32
)
{
    assert!(log_n >= log_cpus);

    let num_cpus = 1 << log_cpus;
    let log_new_n = log_n - log_cpus;
    let mut tmp = vec![vec![F::zero(); 1 << log_new_n]; num_cpus];
    let new_omega = omega.pow(&[num_cpus as u64]);

    worker.scope(0, |scope, _| {
        let a = &*a;

        for (j, tmp) in tmp.iter_mut().enumerate() {
            scope.spawn(move |_| {
                // Shuffle into a sub-FFT
                let omega_j = omega.pow(&[j as u64]);
                let omega_step = omega.pow(&[(j as u64) << log_new_n]);

                let mut elt = F::one();
                for i in 0..(1 << log_new_n) {
                    for s in 0..num_cpus {
                        let idx = (i + (s << log_new_n)) % (1 << log_n);
                        let mut t = a[idx];

                        t.mul_assign(&elt);
                        tmp[i].add_assign(&t);
                        elt.mul_assign(&omega_step);
                    }
                    elt.mul_assign(&omega_j);
                }

                // Perform sub-FFT
                serial_fft(tmp, &new_omega, log_new_n);
            });
        }
    });

    worker.scope(a.len(), |scope, chunk| {
        let tmp = &tmp;

        for (idx, a) in a.chunks_mut(chunk).enumerate() {
            scope.spawn(move |_| {
                let mut idx = idx * chunk;
                let mask = (1 << log_cpus) - 1;
                for a in a {
                    *a = tmp[idx & mask][idx >> log_cpus];
                    idx += 1;
                }
            });
        }
    });
}

fn subview<T: Sized, const N: usize>(slice: &[T]) -> &[[T; N]] {
    assert_eq!(slice.len() % N, 0);
    let len = slice.len() / N;
    let array_slice: &[[T; N]] = unsafe { std::slice::from_raw_parts(slice.as_ptr().cast(), len) };

    array_slice
}

fn subview_mut<T: Sized, const N: usize>(slice: &mut [T]) -> &mut [[T; N]] {
    assert_eq!(slice.len() % N, 0);
    let len = slice.len() / N;
    let array_slice: &mut [[T; N]] = unsafe { std::slice::from_raw_parts_mut(slice.as_mut_ptr().cast(), len) };

    array_slice
}

fn allocate_zeroable_array_without_stack<T: Sized, const N: usize>(size: usize) -> Vec<[T; N]> {
    let allocate_bytes = N * std::mem::size_of::<T>() * size;
    unsafe {
        let t = vec![0u8; allocate_bytes];
        assert_eq!(t.capacity() % std::mem::size_of::<T>() / N, 0);
        assert_eq!(t.capacity() / std::mem::size_of::<T>() / N, size);
        let mut t = std::mem::ManuallyDrop::new(t);
        let t: Vec<[T; N]> = Vec::from_raw_parts(t.as_mut_ptr().cast(),size, size);

        t
    }
}

fn allocate_zeroable<T: Sized>(size: usize) -> Vec<T> {
    let allocate_bytes = std::mem::size_of::<T>() * size;
    unsafe {
        let t = vec![0u8; allocate_bytes];
        assert_eq!(t.capacity() % std::mem::size_of::<T>(), 0);
        assert_eq!(t.capacity() / std::mem::size_of::<T>(), size);
        let mut t = std::mem::ManuallyDrop::new(t);
        let t: Vec<T> = Vec::from_raw_parts(t.as_mut_ptr().cast(),size, size);

        t
    }
}


// we expect that a sub-fft fits into the L2 or L3 cache line
// #[cfg(feature = "const_generic_fft")]
pub(crate) fn parallel_cache_friendly_fft<F: PrimeField, const N: usize>(
    a: &mut [F],
    worker: &Worker,
    omega: &F,
    log_n: u32,
    // log_cpus: u32
)
{
    use std::time::Instant;
    let begin = Instant::now();

    let log_new_n = log2_floor(N);
    assert!(log_n >= log_new_n);

    let log_num_subworks = log_n - log_new_n;
    let num_subworks = 1 << log_num_subworks;

    println!("Processing {} subworks on {} cores", num_subworks, worker.num_cpus());

    let mut tmp = allocate_zeroable_array_without_stack::<F, N>(num_subworks);

    // let mut tmp = vec![[F::zero(); N]; num_subworks];
    let new_omega = omega.pow(&[num_subworks as u64]);

    let mask = (1 << log_n) - 1;

    let start = Instant::now();

    let mut aa = allocate_zeroable::<F>(a.len());
    let l = log_new_n;
    let h = log_num_subworks;

    worker.scope(aa.len(), |scope, chunk| {
        let a = &*a;

        for (idx, aa) in aa.chunks_mut(chunk).enumerate() {
            scope.spawn(move |_| {
                let mut idx = idx * chunk;
                for aa in aa {
                    *aa = a[make_lowest_bits_to_be_highest(idx, l, h) & mask];
                    idx += 1;
                }
            });
        }
    });

    println!("Initial shuffle taken {:?}", start.elapsed());

    println!("Before main loop: {:?}", begin.elapsed());

    let start = Instant::now();

    // we have num_subworks outer loops

    worker.scope(num_subworks, |scope, chunk| {
        let a = &*aa;

        let min_stack_size = N * std::mem::size_of::<F>() * 110 / 100;

        for (chunk_idx, tmp) in tmp.chunks_mut(chunk).enumerate() {
            let start = chunk_idx * chunk;
            let builder = scope.builder().stack_size(min_stack_size).name(format!("FFT thread for {} chunk", chunk_idx));
            builder.spawn(move |_| {
                let mut omega_j = omega.pow(&[start as u64]);
                let mut omega_step = omega.pow(&[(start as u64) << log_new_n]);
                let omega_step_multiplier = omega.pow(&[(1 << log_new_n) as u64]);

                let mut subwork = [F::zero(); N];

                // each inner loop runs over N * num_subworks == n iterations
                for (_j, tmp) in tmp.iter_mut().enumerate() {
                    let start = Instant::now();

                    // Shuffle into a sub-FFT
                    let mut elt = F::one();
                    let mut offset = 0;
                    for i in 0..N {
                        let src = &a[offset..(offset + num_subworks)];
                        for t in src.iter() {
                            let mut t = *t;

                            t.mul_assign(&elt);
                            subwork[i].add_assign(&t);
                            elt.mul_assign(&omega_step);
                        }
                        elt.mul_assign(&omega_j);
                        offset += num_subworks;
                    }

                    if chunk_idx == 0 && _j == 0 {
                        println!("Shuffle taken {:?}", start.elapsed());
                    }

                    let start = Instant::now();

                    // Perform sub-FFT
                    serial_cache_friendly_fft(tmp, &new_omega);

                    if chunk_idx == 0 && _j == 0 {
                        println!("Sub-FFT taken {:?}", start.elapsed());
                    }

                    *tmp = subwork;

                    omega_j.mul_assign(&omega);
                    omega_step.mul_assign(&omega_step_multiplier);
                }
            }).expect("inner thread must spawn successfully");
        }
    });

    println!("Main loop taken {:?}", start.elapsed());

    let start = Instant::now();

    worker.scope(a.len(), |scope, chunk| {
        let tmp = &tmp;
        let mask = (1 << log_num_subworks) - 1;

        for (idx, a) in a.chunks_mut(chunk).enumerate() {
            scope.spawn(move |_| {
                let mut idx = idx * chunk;
                for a in a {
                    *a = tmp[idx & mask][idx >> log_num_subworks];
                    idx += 1;
                }
            });
        }
    });

    println!("Final copy taken {:?}", start.elapsed());
}


// we expect that a sub-fft fits into the L2 or L3 cache line
// #[cfg(feature = "const_generic_fft")]
pub(crate) fn parallel_partitioned_fft<F: PrimeField, const N: usize, const K: usize>(
    a: &mut [F],
    worker: &Worker,
    omega: &F,
    log_n: u32,
)
{
    use std::time::Instant;

    let n = a.len();
    assert_eq!(N*K, n);
    let log_new_n = log2_floor(N);
    let log_num_subworks = log2_floor(K);

    assert!(log_n >= log_new_n);

    // unsafe trick to avoid overflow
    let allocate_bytes = N * std::mem::size_of::<F>() * K;
    let mut tmp: Vec<[F; N]> = unsafe {
        let t = vec![0u8; allocate_bytes];
        let mut t: Vec<[F; N]> = std::mem::transmute(t);
        t.set_len(K);

        t
    };

    let new_omega = omega.pow(&[K as u64]);

    let mask = (1 << log_n) - 1;

    // subshuffle for better memory locality of inner loops
    let mut aa: Vec<F> = unsafe {
        let t = vec![0u8; allocate_bytes];
        let mut t: Vec<F> = std::mem::transmute(t);
        t.set_len(a.len());

        t
    };
    let l = log_new_n;
    let h = log_num_subworks;

    worker.scope(aa.len(), |scope, chunk| {
        let a = &*a;

        for (idx, aa) in aa.chunks_mut(chunk).enumerate() {
            scope.spawn(move |_| {
                let mut idx = idx * chunk;
                for aa in aa {
                    *aa = a[make_lowest_bits_to_be_highest(idx, l, h) & mask];
                    idx += 1;
                }
            });
        }
    });

    let start = Instant::now();

    // we have num_subworks outer loops

    worker.scope(K, |scope, chunk| {
        let a = subview::<F, K>(&aa);

        let min_stack_size = N * std::mem::size_of::<F>() * 110 / 100;

        for (chunk_idx, tmp) in tmp.chunks_mut(chunk).enumerate() {
            let start = chunk_idx * chunk;
            let builder = scope.builder().stack_size(min_stack_size).name(format!("FFT thread for {} chunk", chunk_idx));
            builder.spawn(move |_| {
                let mut omega_j = omega.pow(&[start as u64]);
                let mut omega_step = omega.pow(&[(start as u64) << log_new_n]);
                let omega_step_multiplier = omega.pow(&[(1 << log_new_n) as u64]);

                let mut subwork = [F::zero(); N];

                // each inner loop runs over N * num_subworks == n iterations
                for (_j, tmp) in tmp.iter_mut().enumerate() {
                    let start = Instant::now();
                    // Shuffle into a sub-FFT
                    let mut elt = F::one();
                    for i in 0..N {
                        let src: [F; K] = a[i];
                        for t in src.iter() {
                            let mut t = *t;

                            t.mul_assign(&elt);
                            subwork[i].add_assign(&t);
                            elt.mul_assign(&omega_step);
                        }
                        elt.mul_assign(&omega_j);
                    }

                    if chunk_idx == 0 && _j == 0 {
                        println!("Shuffle taken {:?}", start.elapsed());
                    }

                    let start = Instant::now();

                    // Perform sub-FFT
                    serial_cache_friendly_fft(tmp, &new_omega);

                    if chunk_idx == 0 && _j == 0 {
                        println!("Sub-FFT taken {:?}", start.elapsed());
                    }

                    *tmp = subwork;

                    omega_j.mul_assign(&omega);
                    omega_step.mul_assign(&omega_step_multiplier);
                }
            }).expect("inner thread must spawn successfully");
        }
    });

    println!("Main loop taken {:?}", start.elapsed());

    worker.scope(a.len(), |scope, chunk| {
        let tmp = &tmp;
        let mask = K - 1;

        for (idx, a) in a.chunks_mut(chunk).enumerate() {
            scope.spawn(move |_| {
                let mut idx = idx * chunk;
                for a in a {
                    *a = tmp[idx & mask][idx >> log_num_subworks];
                    idx += 1;
                }
            });
        }
    });
}


// we expect that a sub-fft fits into the L2 or L3 cache line
// #[cfg(feature = "const_generic_fft")]
pub(crate) fn parallel_cache_friendly_fft_without_stack_spilling<F: PrimeField, const N: usize>(
    a: &mut [F],
    worker: &Worker,
    omega: &F,
    log_n: u32,
    // log_cpus: u32
)
{    
    use std::time::Instant;
    let begin = Instant::now();

    assert!(N.is_power_of_two());
    let log_new_n = log2_floor(N);
    assert!(log_n >= log_new_n);

    let log_num_subworks = log_n - log_new_n;
    let num_subworks = 1 << log_num_subworks;

    // let mut tmp = vec![F::zero(); a.len()];
    let allocate_bytes = N * std::mem::size_of::<F>() * num_subworks;
    let mut tmp: Vec<F> = unsafe {
        let t = vec![0u8; allocate_bytes];
        let mut t: Vec<F> = std::mem::transmute(t);
        t.set_len(a.len());
        // capacity is garbage here, but we do not care as we do not reallocate it ever

        t
    };

    let new_omega = omega.pow(&[num_subworks as u64]);

    let mask = (1 << log_n) - 1;

    // subshuffle
    let mut aa: Vec<F> = unsafe {
        let t = vec![0u8; allocate_bytes];
        let mut t: Vec<F> = std::mem::transmute(t);
        t.set_len(a.len());
        // capacity is garbage here, but we do not care as we do not reallocate it ever

        t
    };
    let l = log_new_n;
    let h = log_num_subworks;

    let start = Instant::now();

    worker.scope(aa.len(), |scope, chunk| {
        let a = &*a;

        for (idx, aa) in aa.chunks_mut(chunk).enumerate() {
            scope.spawn(move |_| {
                let mut idx = idx * chunk;
                for aa in aa {
                    *aa = a[make_lowest_bits_to_be_highest(idx, l, h) & mask];
                    idx += 1;
                }
            });
        }
    });

    println!("Initial shuffle taken {:?}", start.elapsed());

    println!("Before main loop: {:?}", begin.elapsed());

    let start = Instant::now();

    worker.scope(num_subworks, |scope, chunk| {
        let a = &aa;

        let tmp = subview_mut::<F, N>(&mut tmp);

        for (chunk_idx, tmp) in tmp.chunks_mut(chunk).enumerate() {
            let start = chunk_idx * chunk;
            scope.spawn(move |_| {
                let mut omega_j = omega.pow(&[start as u64]);
                let mut omega_step = omega.pow(&[(start as u64) << log_new_n]);
                let omega_step_multiplier = omega.pow(&[(1 << log_new_n) as u64]);

                for (_j, tmp) in tmp.iter_mut().enumerate() {
                    let start = Instant::now();

                    // Shuffle into a sub-FFT
                    let mut elt = F::one();
                    for i in 0..N {
                        let start = i << log_new_n;
                        for s in 0..num_subworks {
                            // if we would work over a then
                            // we fix i \in 0..N (lowest log_new_n bits) 
                            // and then jump with intervals N (highest bits), that
                            // that is suboptimal from memory locality perspective
                            // let idx = (i + (s << log_new_n)) & mask;
                            // let mut t = a[idx];

                            // so we did our shuffle of aa and got this
                            let idx = (start + s) & mask;
                            let mut t = a[idx];

                            t.mul_assign(&elt);
                            tmp[i].add_assign(&t);
                            elt.mul_assign(&omega_step);
                        }
                        elt.mul_assign(&omega_j);
                    }

                    if chunk_idx == 0 && _j == 0 {
                        println!("Shuffle taken {:?}", start.elapsed());
                    }

                    let start = Instant::now();

                    // Perform sub-FFT
                    serial_fft(tmp, &new_omega, log_new_n);

                    if chunk_idx == 0 && _j == 0 {
                        println!("Sub-FFT taken {:?}", start.elapsed());
                    }

                    omega_j.mul_assign(&omega);
                    omega_step.mul_assign(&omega_step_multiplier);
                }
            });
        }
    });

    println!("Main loop taken {:?}", start.elapsed());

    let start = Instant::now();

    worker.scope(a.len(), |scope, chunk| {
        let tmp = subview::<F, N>(&tmp);
        let mask = (1 << log_num_subworks) - 1;

        for (idx, a) in a.chunks_mut(chunk).enumerate() {
            scope.spawn(move |_| {
                let mut idx = idx * chunk;
                for a in a {
                    *a = tmp[idx & mask][idx >> log_num_subworks];
                    idx += 1;
                }
            });
        }
    });

    println!("Final copy taken {:?}", start.elapsed());
}

#[inline(always)]
fn bitreverse_only_lowest(n: u32, l: u32) -> u32 {
    let shift_back = 32u32 - l;
    let fully_reversed = n.reverse_bits() >> shift_back;
    let mask = (1u32 << l) - 1u32;
    let result = (fully_reversed & mask) | (n & (!mask));

    result
}

#[inline(always)]
fn make_lowest_bits_to_be_highest(n: usize, l: u32, h: u32) -> usize {
    // we make lowest l bits now as highest l bits,
    // and highest h bits (assuming number fits into l+h bits) as lowest bits
    let mask_low = (1usize << l) - 1;

    ((n & mask_low) << h) | (n >> l)
}

pub(crate) fn serial_cache_friendly_fft<F: PrimeField, const N: usize>(
    a: &mut [F; N], 
    omega: &F, 
) {
    #[inline(always)]
    fn bitreverse(mut n: u32, l: u32) -> u32 {
        let mut r = 0;
        for _ in 0..l {
            r = (r << 1) | (n & 1);
            n >>= 1;
        }
        r
    }

    let n = a.len();
    assert_eq!(n, N);

    let log_n = log2_floor(N);

    for k in 0..n {
        let rk = bitreverse(k as u32, log_n) as usize;
        if k < rk {
            a.swap(rk as usize, k as usize);
        }
    }

    let mut m = 1;
    for _ in 0..log_n {
        let w_m = omega.pow(&[(N / (2*m)) as u64]);

        let mut k = 0;
        while k < n {
            let mut w = F::one();
            for j in 0..m {
                
                let mut t = a[(k+j+m) as usize];
                t.mul_assign(&w);
                let mut tmp = a[(k+j) as usize];
                tmp.sub_assign(&t);
                a[(k+j+m) as usize] = tmp;
                a[(k+j) as usize].add_assign(&t);
                w.mul_assign(&w_m);
            }

            k += 2*m;
        }
        
        m *= 2;
    }
}