use super::*;

#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
pub(crate) fn add_nocarry_impl(a: [u64; 2], b: [u64; 2]) -> [u64; 2] {
    // we use ADD/ADC only here as it's the same latency, but shorter encoding
    
    let [mut r0, mut r1] = a;
    let [b0, b1] = b;

    unsafe {
        asm!(
            "add {in_out_0}, {b_0}",
            "adc {in_out_1}, {b_1}",
            in_out_0 = inlateout(reg) r0,
            in_out_1 = inlateout(reg) r1,
            b_0 = in(reg) b0,
            b_1 = in(reg) b1,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1]
}

// #[allow(clippy::too_many_lines)]
// #[inline(always)]
// #[cfg(target_arch = "x86_64")]
// pub(crate) fn sub_noborrow_impl(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
//     let mut r0: u64;
//     let mut r1: u64;
//     let mut r2: u64;
//     let mut r3: u64;

//     unsafe {
//         asm!(
//             "xor r12d, r12d",
//             "mov r12, qword ptr [{a_ptr} + 0]",
//             "sub r12, qword ptr [{b_ptr} + 0]",
//             "mov r13, qword ptr [{a_ptr} + 8]",
//             "sbb r13, qword ptr [{b_ptr} + 8]",
//             "mov r14, qword ptr [{a_ptr} + 16]",
//             "sbb r14, qword ptr [{b_ptr} + 16]",
//             "mov r15, qword ptr [{a_ptr} + 24]",
//             "sbb r15, qword ptr [{b_ptr} + 24]",
//             a_ptr = in(reg) a.as_ptr(),
//             b_ptr = in(reg) b.as_ptr(),
//             out("r12") r0, 
//             out("r13") r1, 
//             out("r14") r2, 
//             out("r15") r3,
//             options(pure, readonly, nostack)
//         );
//     }

//     [r0, r1, r2, r3]
// }


// #[allow(clippy::too_many_lines)]
// #[inline(always)]
// #[cfg(target_arch = "x86_64")]
// pub(crate) fn double_nocarry_impl(a: &[u64; 4]) -> [u64; 4] {
//     // we use ADD/ADC only here as it's the same latency, but shorter encoding
    
//     let mut r0: u64;
//     let mut r1: u64;
//     let mut r2: u64;
//     let mut r3: u64;

//     unsafe {
//         asm!(
//             "xor r12d, r12d",
//             "mov r12, qword ptr [{a_ptr} + 0]",
//             "add r12, r12",
//             "mov r13, qword ptr [{a_ptr} + 8]",
//             "adc r13, r13",
//             "mov r14, qword ptr [{a_ptr} + 16]",
//             "adc r14, r14",
//             "mov r15, qword ptr [{a_ptr} + 24]",
//             "adc r15, r15",
//             a_ptr = in(reg) a.as_ptr(),
//             out("r12") r0, 
//             out("r13") r1, 
//             out("r14") r2, 
//             out("r15") r3,
//             options(pure, readonly, nostack)
//         );
//     }

//     [r0, r1, r2, r3]
// }

// #[allow(clippy::too_many_lines)]
// #[inline(always)]
// #[cfg(target_arch = "x86_64")]
// // #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
// pub(crate) fn add_with_reduction_impl(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
//     let mut r0: u64;
//     let mut r1: u64;
//     let mut r2: u64;
//     let mut r3: u64;

//     unsafe {
//         asm!(
//             // we sum (a+b) using addition chain with OF
//             // and sum (a+b) - p using addition chain with CF
//             // if (a+b) does not overflow the modulus
//             // then sum (a+b) - p will produce CF
//             "xor r12d, r12d",
//             "mov r12, qword ptr [{a_ptr} + 0]",
//             "mov r13, qword ptr [{a_ptr} + 8]",
//             "mov r14, qword ptr [{a_ptr} + 16]",
//             "mov r15, qword ptr [{a_ptr} + 24]",
//             "adox r12, qword ptr [{b_ptr} + 0]",
//             "mov r8, r12",
//             "adcx r8, qword ptr [rip + {q0_ptr}]",
//             "adox r13, qword ptr [{b_ptr} + 8]",
//             "mov r9, r13",
//             "adcx r9, qword ptr [rip + {q1_ptr}]",
//             "adox r14, qword ptr [{b_ptr} + 16]",
//             "mov r10, r14",
//             "adcx r10, qword ptr [rip + {q2_ptr}]",
//             "adox r15, qword ptr [{b_ptr} + 24]",
//             "mov r11, r15",
//             "adcx r11, qword ptr [rip + {q3_ptr}]",

//             // if CF = 0 then take value (a+b) from [r12, .., r15]
//             // otherwise take (a+b) - p

//             "cmovc r12, r8",
//             "cmovc r13, r9",
//             "cmovc r14, r10",
//             "cmovc r15, r11",  

//             q0_ptr = sym MODULUS_NEGATED_STATIC_0,
//             q1_ptr = sym MODULUS_NEGATED_STATIC_1,
//             q2_ptr = sym MODULUS_NEGATED_STATIC_2,
//             q3_ptr = sym MODULUS_NEGATED_STATIC_3,
//             // end of reduction
//             a_ptr = in(reg) a.as_ptr(),
//             b_ptr = in(reg) b.as_ptr(),
//             out("r8") _, 
//             out("r9") _, 
//             out("r10") _, 
//             out("r11") _, 
//             out("r12") r0, 
//             out("r13") r1, 
//             out("r14") r2, 
//             out("r15") r3,
//             options(pure, readonly, nostack)
//         );
//     }

//     // unsafe {
//     //     asm!(
//     //         "xor r12d, r12d",
//     //         "mov r12, qword ptr [{a_ptr} + 0]",
//     //         "mov r13, qword ptr [{a_ptr} + 8]",
//     //         "mov r14, qword ptr [{a_ptr} + 16]",
//     //         "mov r15, qword ptr [{a_ptr} + 24]",
//     //         "add r12, qword ptr [{b_ptr} + 0]",
//     //         "adc r13, qword ptr [{b_ptr} + 8]",
//     //         "adc r14, qword ptr [{b_ptr} + 16]",
//     //         "adc r15, qword ptr [{b_ptr} + 24]",

//     //         "mov r8, r12",
//     //         "mov rdx, qword ptr [rip + {q0_ptr}]",
//     //         "sub r8, rdx",
//     //         "mov r9, r13",
//     //         "mov rdx, qword ptr [rip + {q1_ptr}]",
//     //         "sbb r9, rdx",
//     //         "mov r10, r14",
//     //         "mov rdx, qword ptr [rip + {q2_ptr}]",
//     //         "sbb r10, rdx",
//     //         "mov r11, r15",
//     //         "mov rdx, qword ptr [rip + {q3_ptr}]",
//     //         "sbb r11, rdx",

//     //         "cmovnc r12, r8",
//     //         "cmovnc r13, r9",
//     //         "cmovnc r14, r10",
//     //         "cmovnc r15, r11",  

//     //         q0_ptr = sym #m0,
//     //         q1_ptr = sym #m1,
//     //         q2_ptr = sym #m2,
//     //         q3_ptr = sym #m3,
//     //         // end of reduction
//     //         a_ptr = in(reg) a.as_ptr(),
//     //         b_ptr = in(reg) b.as_ptr(),
//     //         out("rdx") _, 
//     //         out("r8") _, 
//     //         out("r9") _, 
//     //         out("r10") _, 
//     //         out("r11") _, 
//     //         out("r12") r0, 
//     //         out("r13") r1, 
//     //         out("r14") r2, 
//     //         out("r15") r3,
//     //         options(pure, readonly, nostack)
//     //     );
//     // }

//     [r0, r1, r2, r3]
// }

// #[allow(clippy::too_many_lines)]
// #[inline(always)]
// #[cfg(target_arch = "x86_64")]
// // #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
// pub(crate) fn double_with_reduction_impl(a: &[u64; 4]) -> [u64; 4] {
//     let mut r0: u64;
//     let mut r1: u64;
//     let mut r2: u64;
//     let mut r3: u64;

//     unsafe {
//         asm!(
//             // we sum (a+b) using addition chain with OF
//             // and sum (a+b) - p using addition chain with CF
//             // if (a+b) does not overflow the modulus
//             // then sum (a+b) will produce CF
//             "xor r12d, r12d",
//             "mov r12, qword ptr [{a_ptr} + 0]",
//             "mov r13, qword ptr [{a_ptr} + 8]",
//             "mov r14, qword ptr [{a_ptr} + 16]",
//             "mov r15, qword ptr [{a_ptr} + 24]",
//             "adox r12, r12",
//             "mov r8, r12",
//             "adcx r8, qword ptr [rip + {q0_ptr}]",
//             "adox r13, r13",
//             "mov r9, r13",
//             "adcx r9, qword ptr [rip + {q1_ptr}]",
//             "adox r14, r14",
//             "mov r10, r14",
//             "adcx r10, qword ptr [rip + {q2_ptr}]",
//             "adox r15, r15",
//             "mov r11, r15",
//             "adcx r11, qword ptr [rip + {q3_ptr}]",

//             // if CF = 0 then take value (a+b) from [r12, .., r15]
//             // otherwise take (a+b) - p

//             "cmovc r12, r8",
//             "cmovc r13, r9",
//             "cmovc r14, r10",
//             "cmovc r15, r11",  

//             q0_ptr = sym MODULUS_NEGATED_STATIC_0,
//             q1_ptr = sym MODULUS_NEGATED_STATIC_1,
//             q2_ptr = sym MODULUS_NEGATED_STATIC_2,
//             q3_ptr = sym MODULUS_NEGATED_STATIC_3,
//             // end of reduction
//             a_ptr = in(reg) a.as_ptr(),
//             out("r8") _, 
//             out("r9") _, 
//             out("r10") _, 
//             out("r11") _, 
//             out("r12") r0, 
//             out("r13") r1, 
//             out("r14") r2, 
//             out("r15") r3,
//             options(pure, readonly, nostack)
//         );
//     }

//     [r0, r1, r2, r3]
// }

// #[allow(clippy::too_many_lines)]
// #[inline(always)]
// #[cfg(target_arch = "x86_64")]
// // #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
// pub(crate) fn sub_with_reduction_impl(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
//     let mut r0: u64;
//     let mut r1: u64;
//     let mut r2: u64;
//     let mut r3: u64;

//     unsafe {
//         asm!(
//             "xor r12d, r12d",
//             "mov r12, qword ptr [{a_ptr} + 0]",
//             "sub r12, qword ptr [{b_ptr} + 0]",
//             "mov r13, qword ptr [{a_ptr} + 8]",
//             "sbb r13, qword ptr [{b_ptr} + 8]",
//             "mov r14, qword ptr [{a_ptr} + 16]",
//             "sbb r14, qword ptr [{b_ptr} + 16]",
//             "mov r15, qword ptr [{a_ptr} + 24]",
//             "sbb r15, qword ptr [{b_ptr} + 24]",

//             // duplicate (a-b) into [r8, r9, r10, r11]

//             // now make [r12, .., r15] + modulus;
//             // if (a-b) did underflow then (a-b) + modulus > 2^256
//             // so below we will get an overflow

//             "mov r8, r12",
//             "add r12, 1",
//             // "adox r12, qword ptr [rip + {q0_ptr}]",
//             "mov r9, r13",
//             "adc r13, 0",
//             // "adox r13, qword ptr [rip + {q1_ptr}]",
//             "mov r10, r14",
//             "adc r14, 0",
//             // "adox r14, qword ptr [rip + {q2_ptr}]",
//             "mov r11, r15",
//             "adc r15, qword ptr [rip + {q3_ptr}]",
//             // "adox r15, qword ptr [rip + {q3_ptr}]",

//             // if no carry (no overflow) then what we had in [r8, r11] is valid
//             "cmovnc r12, r8",
//             "cmovnc r13, r9",
//             "cmovnc r14, r10",
//             "cmovnc r15, r11", 
//             // end of reduction
//             q3_ptr = sym MODULUS_3_STATIC,
//             a_ptr = in(reg) a.as_ptr(),
//             b_ptr = in(reg) b.as_ptr(),
//             out("r8") _, 
//             out("r9") _, 
//             out("r10") _, 
//             out("r11") _, 
//             out("r12") r0, 
//             out("r13") r1, 
//             out("r14") r2, 
//             out("r15") r3,
//             options(pure, readonly, nostack)
//         );
//     }

//     [r0, r1, r2, r3]
// }

// #[allow(clippy::too_many_lines)]
// #[inline(always)]
// #[cfg(target_arch = "x86_64")]
// pub(crate) fn add_modulus_impl(a: &[u64; 4]) -> [u64; 4] {
//     // we use ADD/ADC only here as it's the same latency, but shorter encoding
    
//     let mut r0: u64;
//     let mut r1: u64;
//     let mut r2: u64;
//     let mut r3: u64;

//     unsafe {
//         asm!(
//             "xor r12d, r12d",
//             "mov r12, qword ptr [{a_ptr} + 0]",
//             "add r12, 1",
//             "mov r13, qword ptr [{a_ptr} + 8]",
//             "adc r13, 0",
//             "mov r14, qword ptr [{a_ptr} + 16]",
//             "adc r14, 0",
//             "mov r15, qword ptr [{a_ptr} + 24]",
//             "adc r15, qword ptr [rip + {q3_ptr}]",
//             a_ptr = in(reg) a.as_ptr(),
//             q3_ptr = sym MODULUS_3_STATIC,
//             out("r12") r0, 
//             out("r13") r1, 
//             out("r14") r2, 
//             out("r15") r3,
//             options(pure, readonly, nostack)
//         );
//     }

//     [r0, r1, r2, r3]
// }


// #[allow(clippy::too_many_lines)]
// #[inline(always)]
// #[cfg(target_arch = "x86_64")]
// pub(crate) fn sub_modulus_impl(a: &[u64; 4]) -> [u64; 4] {
//     // we use SUB/SBB to utilize 0s
    
//     let mut r0: u64;
//     let mut r1: u64;
//     let mut r2: u64;
//     let mut r3: u64;

//     unsafe {
//         asm!(
//             "xor r12d, r12d",
//             "mov r12, qword ptr [{a_ptr} + 0]",
//             "sub r12, 1",
//             "mov r13, qword ptr [{a_ptr} + 8]",
//             "sub r13, 0",
//             "mov r14, qword ptr [{a_ptr} + 16]",
//             "sub r14, 0",
//             "mov r15, qword ptr [{a_ptr} + 24]",
//             "sbb r15, qword ptr [rip + {q3_ptr}]",
//             a_ptr = in(reg) a.as_ptr(),
//             q3_ptr = sym MODULUS_3_STATIC,
//             out("r12") r0, 
//             out("r13") r1, 
//             out("r14") r2, 
//             out("r15") r3,
//             options(pure, readonly, nostack)
//         );
//     }

//     [r0, r1, r2, r3]
// }


// #[allow(clippy::too_many_lines)]
// #[inline(always)]
// #[cfg(target_arch = "x86_64")]
// // #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
// pub(crate) fn reduce_by_modulus_impl(a: &[u64; 4]) -> [u64; 4] {
//     let mut r0: u64;
//     let mut r1: u64;
//     let mut r2: u64;
//     let mut r3: u64;

//     unsafe {
//         asm!(
//             "xor r12d, r12d",
//             "mov r12, qword ptr [{a_ptr} + 0]",
//             "mov r8, r12",
//             "sub r12, 1",
//             "mov r13, qword ptr [{a_ptr} + 8]",
//             "mov r9, r13",
//             "sbb r13, 0",
//             "mov r14, qword ptr [{a_ptr} + 16]",
//             "mov r10, r14",
//             "sbb r14, 0",
//             "mov r15, qword ptr [{a_ptr} + 24]",
//             "mov r11, r15",
//             "sbb r15, qword ptr [rip + {q3_ptr}]",

//             "cmovc r12, r8",
//             "cmovc r13, r9",
//             "cmovc r14, r10",
//             "cmovc r15, r11", 

//             q3_ptr = sym MODULUS_3_STATIC,
//             a_ptr = in(reg) a.as_ptr(),
//             out("r8") _, 
//             out("r9") _, 
//             out("r10") _, 
//             out("r11") _, 
//             out("r12") r0, 
//             out("r13") r1, 
//             out("r14") r2, 
//             out("r15") r3,
//             options(pure, readonly, nostack)
//         );
//     }

//     [r0, r1, r2, r3]
// }


#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
pub(crate) fn mont_mul_with_reduction_impl_by_ref(a: &[u64; 2], b: &[u64; 2]) -> [u64; 2] {
    // we also manually use the facts about INV = -1 and modulus[0] = 1

    let mut r0: u64;
    let mut r1: u64;

    // round 0
    // |    |    | t1 | t0 |
    // |    | t3 | t2 |    |
    // round 1
    // |    | t3 | t1 |    |
    // |    | t2 | t0 |    |
    // | r1 | r0 |    |    |

    unsafe {
        asm!(
            "mov rdx, qword ptr [{a_ptr} + 0]",
            "mulx {t1}, {t0}, qword ptr [{b_ptr} + 0]",
            // a0 * b1
            "mulx {t3}, {t2}, qword ptr [{b_ptr} + 8]",
            // my now we wait to get m = -t0
            // instead of later adding t0 and -t0, we just do set CF"
            "mov rdx, {t0}", // {t0} is reusable
            "neg rdx", // it's m = k * t0
            "xor {in_out_0}, {in_out_0}",
            "mulx {in_out_0}, {t0}, qword ptr [rip + {q1_ptr}]", // second part of first reduction
            "stc",
            "adox {t1}, {t2}", // {t2} is now reusable,
            "adox {t3}, qword ptr [rip + {zero_ptr}]", // complete the carry flag
            "adcx {t1}, {t0}",
            "adcx {t3}, {in_out_0}",

            // round 1
            "mov rdx, qword ptr [{a_ptr} + 8]",
            "mulx {t2}, {t0}, qword ptr [{b_ptr} + 0]",
            "mulx {in_out_1}, {in_out_0}, qword ptr [{b_ptr} + 8]",
            "xor rdx, rdx",
            // t0 + t1
            "adox {t1}, {t0}", // r0
            // carry into t2 without overflow
            "adox {t2}, rdx", // carry
            "mov rdx, {t1}",
            "neg rdx",
            "xor {t0}, {t0}", // reset flags
            "mulx {t1}, {t0}, qword ptr [rip + {q1_ptr}]", // second part of first reduction
            "stc",
            //
            "adox {in_out_0}, {t2}", // carry + low
            "adox {in_out_1}, qword ptr [rip + {zero_ptr}]", // propagate
            "adox {in_out_0}, {t3}", // r1
            "adox {in_out_1}, qword ptr [rip + {zero_ptr}]", // carry
            // m * q1
            "adcx {in_out_0}, {t0}", // t0 + red_carry + r1
            "adcx {in_out_1}, {t1}", // t1 + propagate + carry

            // lowest limb is in 
            "mov {t2}, {in_out_0}",
            "sub {in_out_0}, 1",
            "mov {t3}, {in_out_1}",
            "sbb {in_out_1}, qword ptr [rip + {q1_ptr}]",

            // if CF == 0 then original result was ok (reduction wa not necessary)
            // so if carry (CMOVNQ) then we copy 
            "cmovc {in_out_0}, {t2}",
            "cmovc {in_out_1}, {t3}",
            // end of reduction
            zero_ptr = sym ZERO_U64,
            q1_ptr = sym MODULUS_1_STATIC,
            a_ptr = in(reg) a.as_ptr(),
            b_ptr = in(reg) b.as_ptr(),
            in_out_0 = out(reg) r0,
            in_out_1 = out(reg) r1,
            t0 = out(reg) _,
            t1 = out(reg) _,
            t2 = out(reg) _,
            t3 = out(reg) _,
            out("rdx") _,
            options(pure, readonly, nostack)
        );
    }

    // println!("Tmp = {}", tmp);

    [r0, r1]
}

// #[allow(clippy::too_many_lines)]
// #[inline(always)]
// #[cfg(target_arch = "x86_64")]
// // #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
// pub(crate) fn mont_mul_with_reduction_impl(a: [u64; 2], b: [u64; 2]) -> [u64; 2] {
//     // we also manually use the facts about INV = -1 and modulus[0] = 1

//     let [mut r0, mut r1] = a;
//     let [b0, b1] = b;

//     // let mut tmp: u64;

//     // round 0
//     // |    |    | t1 | t0 |
//     // |    | t3 | t2 |    |
//     // round 1
//     // |    | t3 | t1 |    |
//     // |    | t2 | t0 |    |
//     // | r1 | r0 |    |    |

//     unsafe {
//         asm!(
//             // we have a0 in rdx already
//             // and b0 in {b_0}
//             // it's good for us to run everything in parallel
//             // round a0 * b0
//             "mov rdx, {in_out_0}",
//             "mulx {t1}, {t0}, {b_0}",
//             // a0 * b1
//             "mulx {t3}, {t2}, {b_1}",
//             // my now we wait to get m = -t0
//             "mov rdx, {t0}",
//             "neg rdx", // it's m = k * t0
//             "xor {in_out_0}, {in_out_0}",
//             // carry results of multiplication by a0
//             "adcx {t1}, {t2}", // {t2} is now reusable,
//             "adcx {t3}, {in_out_0}", // complete the carry flag
//             // we do not need to multiply by q0 limb of the modulus, os just add {t0} and m
//             // m * 1 + t0 -> t0 + OF
//             "adox {t0}, rdx", // {t0} is now reusable
//             // m * q1
//             "mulx {t2}, {t0}, qword ptr [rip + {q1_ptr}]", // second part of first reduction
//             // immediately start the next round of multiplications
//             // t3 + (m * q1) + OF
//             "adox {t1}, {t0}",
//             "adox {t3}, {t2}",

//             // round 1
//             "mov rdx, {in_out_1}",
//             "mulx {t2}, {t0}, {b_0}",
//             "mulx {in_out_1}, {in_out_0}, {b_1}",
//             "xor rdx, rdx",
//             // t0 + t1
//             "adcx {t1}, {t0}", // r0
//             // carry into t2 without overflow
//             "adcx {t2}, rdx", // carry
//             "mov rdx, {t1}",
//             "neg rdx",
//             "xor {t0}, {t0}", // reset flags
//             // add into t1
//             "adox {t1}, rdx", // red_carry in OF
//             //
//             "adcx {in_out_0}, {t2}", // carry + low
//             "adcx {in_out_1}, {t0}", // propagate
//             "adcx {in_out_0}, {t3}", // r1
//             "adcx {in_out_1}, {t0}", // carry
//             // m * q1
//             "mulx {t2}, {t0}, qword ptr [rip + {q1_ptr}]", // second part of first reduction
//             "adox {in_out_0}, {t0}", // t0 + red_carry + r1
//             "adox {in_out_1}, {t2}", // t2 + propagate + carry

//             // "mov {t3}, qword ptr [rip + {q1_ptr}]",

//             // // lowest limb is in 
//             // "mov {t0}, {in_out_0}",
//             // "sub {in_out_0}, 1",
//             // "mov {t1}, {in_out_1}",
//             // "sbb {in_out_1}, {t3}",

//             // // if CF == 0 then original result was ok (reduction wa not necessary)
//             // // so if carry (CMOVNQ) then we copy 
//             // "cmovc {in_out_0}, {t0}",
//             // "cmovc {in_out_1}, {t1}",
//             // end of reduction
//             q1_ptr = sym MODULUS_1_STATIC,
//             in_out_0 = inlateout(reg) r0,
//             in_out_1 = inlateout(reg) r1,
//             b_0 = in(reg) b0,
//             b_1 = in(reg) b1,
//             t0 = out(reg) _,
//             t1 = out(reg) _,
//             t2 = out(reg) _,
//             t3 = out(reg) _,
//             // t4 = out(reg) _,
//             // t5 = out(reg) _,
//             // t6 = out(reg) _,
//             // t7 = out(reg) _,
//             // tmp = out(reg) tmp,
//             out("rdx") _,
//             options(pure, readonly, nostack)
//         );
//     }

//     // println!("Tmp = {}", tmp);

//     [r0, r1]
// }


#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
pub(crate) fn mont_mul_with_reduction_impl(a: [u64; 2], b: [u64; 2]) -> [u64; 2] {
    // we also manually use the facts about INV = -1 and modulus[0] = 1

    let [mut r0, mut r1] = a;
    let [b0, b1] = b;

    // let mut tmp: u64;

    // round 0
    // |    |    | t1 | t0 |
    // |    | t3 | t2 |    |
    // round 1
    // |    | t3 | t1 |    |
    // |    | t2 | t0 |    |
    // | r1 | r0 |    |    |

    unsafe {
        asm!(
            "mov rdx, {in_out_0}",
            "mulx {t1}, {t0}, {b_0}",
            // a0 * b1
            "mulx {t3}, {t2}, {b_1}",
            // my now we wait to get m = -t0
            // instead of later adding t0 and -t0, we just do set CF"
            "mov rdx, {t0}", // {t0} is reusable
            "neg rdx", // it's m = k * t0
            "xor {in_out_0}, {in_out_0}",
            "mulx {in_out_0}, {t0}, qword ptr [rip + {q1_ptr}]", // second part of first reduction
            "stc",
            "adox {t1}, {t2}", // {t2} is now reusable,
            "adox {t3}, qword ptr [rip + {zero_ptr}]", // complete the carry flag
            "adcx {t1}, {t0}",
            "adcx {t3}, {in_out_0}",

            // round 1
            "mov rdx, {in_out_1}",
            "mulx {t2}, {t0}, {b_0}",
            "mulx {in_out_1}, {in_out_0}, {b_1}",
            "xor rdx, rdx",
            // t0 + t1
            "adox {t1}, {t0}", // r0
            // carry into t2 without overflow
            "adox {t2}, rdx", // carry
            "mov rdx, {t1}",
            "neg rdx",
            "xor {t0}, {t0}", // reset flags
            "mulx {t1}, {t0}, qword ptr [rip + {q1_ptr}]", // second part of first reduction
            "stc",
            //
            "adox {in_out_0}, {t2}", // carry + low
            "adox {in_out_1}, qword ptr [rip + {zero_ptr}]", // propagate
            "adox {in_out_0}, {t3}", // r1
            "adox {in_out_1}, qword ptr [rip + {zero_ptr}]", // carry
            // m * q1
            "adcx {in_out_0}, {t0}", // t0 + red_carry + r1
            "adcx {in_out_1}, {t1}", // t1 + propagate + carry

            "mov {t3}, qword ptr [rip + {q1_ptr}]",

            // lowest limb is in 
            "mov {t0}, {in_out_0}",
            "sub {in_out_0}, 1",
            "mov {t1}, {in_out_1}",
            "sbb {in_out_1}, {t3}",

            // if CF == 0 then original result was ok (reduction wa not necessary)
            // so if carry (CMOVNQ) then we copy 
            "cmovc {in_out_0}, {t0}",
            "cmovc {in_out_1}, {t1}",
            // end of reduction
            zero_ptr = sym ZERO_U64,
            q1_ptr = sym MODULUS_1_STATIC,
            in_out_0 = inlateout(reg) r0,
            in_out_1 = inlateout(reg) r1,
            b_0 = in(reg) b0,
            b_1 = in(reg) b1,
            t0 = out(reg) _,
            t1 = out(reg) _,
            t2 = out(reg) _,
            t3 = out(reg) _,
            // t4 = out(reg) _,
            // t5 = out(reg) _,
            // t6 = out(reg) _,
            // t7 = out(reg) _,
            // tmp = out(reg) tmp,
            out("rdx") _,
            options(pure, readonly, nostack)
        );
    }

    // println!("Tmp = {}", tmp);

    [r0, r1]
}


#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
pub(crate) fn mont_mul_with_reduction_impl_schoolbook(a: [u64; 2], b: [u64; 2]) -> [u64; 2] {
    // we also manually use the facts about INV = -1 and modulus[0] = 1

    let [mut r0, mut r1] = a;
    let [b0, b1] = b;

    // do it the intel way

    // |    |    | t1 | t0 |
    // |    | t3 | t2 |    |
    // |    | t5 | t4 |    |
    // | t7 | t6 |    |    |
    unsafe {
        asm!(
            "mov rdx, {in_out_0}",
            "xor {in_out_0}, {in_out_0}",
            "mulx {t1}, {t0}, {b_0}",
            "mulx {t3}, {t2}, {b_1}",
            "mov rdx, {in_out_1}",
            "mulx {t5}, {t4}, {b_0}",
            "mulx {t7}, {t6}, {b_1}",

            "mov rdx, {t0}",
            "neg rdx",
            "xor {in_out_0}, {in_out_0}",

            "adox {t1}, {t2}", // limb 1, OF chain
            "adcx {t3}, {t5}", // limb 2, CF chain
            "adox {t3}, {t6}", // limb 2, OF chain
            "adcx {t7}, {in_out_0}", // limb 3, CF chain
            "adox {t7}, {in_out_0}", // limb 3, OF chain

            // // OF chain is for multiplication, CF - for reduction
            // "adox {t1}, {t2}", // limb 1, 
            // "adox {t3}, {t6}", // limb 2, 
            // "adox {t7}, qword ptr [rip + {zero_ptr}]", // limb 3, 
            // "adox {t1}, {t4}", // limb 1,
            // "adox {t3}, {t5}", // limb 2,
            // "adox {t7}, qword ptr [rip + {zero_ptr}]", // limb 3,

            // montgomery reduction

            // round 0

            // -1(INV) * {t0} * 1 (modulus[0]) + t0
            "stc", // virtual one 
            // -1(INV) * {t0} = m; m * modulus[1] + limb 1 + CF
            "mulx {in_out_1}, {in_out_0}, qword ptr [rip + {q1_ptr}]",
            "mov rdx, {t1}",
            "neg rdx",
            "xor {in_out_0}, {in_out_0}",
            "adox {t1}, rdx", // we do it for OF 
            "mulx {in_out_1}, {in_out_0}, qword ptr [rip + {q1_ptr}]",
            "adcx {t1}, {in_out_0}",
            "adcx {t3}, {in_out_1}", // into limb 2
            "adox {t3}, {in_out_0}",
            "adox {t7}, {in_out_1}", // into limb 3

            // result is in [{t3}, {t7}]
            "mov {in_out_0}, {t3}",
            "sub {t3}, 1",
            "mov {in_out_1}, {t7}",
            "sbb {t7}, qword ptr [rip + {q1_ptr}]",

            // // round 0

            // // -1(INV) * {t0} * 1 (modulus[0]) + t0
            // "stc", // virtual one 
            // // -1(INV) * {t0} = m; m * modulus[1] + limb 1 + CF
            // "mulx {in_out_1}, {in_out_0}, qword ptr [rip + {q1_ptr}]",
            // "adcx {t1}, {in_out_0}",
            // "adcx {t3}, {in_out_1}", // into limb 2
            // // CF = 0

            // "mov rdx, {t1}",
            // "neg rdx",
            // "xor {in_out_0}, {in_out_0}",

            // // round 1
            // // -1(INV) * {t1} * 1 (modulus[0]) + t1
            // "stc", // virtual one 
            // // -1(INV) * {t1} = m; m * modulus[1] + limb 2 + CF
            // "mulx {in_out_1}, {in_out_0}, qword ptr [rip + {q1_ptr}]",
            // "adcx {t3}, {in_out_0}",
            // "adcx {t7}, {in_out_1}", // into limb 3
            // // CF = 0

            // // result is in [{t3}, {t7}]
            // "mov {in_out_0}, {t3}",
            // "sub {t3}, 1",
            // "mov {in_out_1}, {t7}",
            // "sbb {t7}, qword ptr [rip + {q1_ptr}]",

            // if CF == 0 then reduced result is ok,
            // so if no carry (CMOVNQ) then we copy 
            "cmovc {in_out_0}, {t3}",
            "cmovc {in_out_1}, {t7}",
            // end of reduction
            // zero_ptr = sym ZERO_U64,
            q1_ptr = sym MODULUS_1_STATIC,
            in_out_0 = inlateout(reg) r0,
            in_out_1 = inlateout(reg) r1,
            b_0 = in(reg) b0,
            b_1 = in(reg) b1,
            t0 = out(reg) _,
            t1 = out(reg) _,
            t2 = out(reg) _,
            t3 = out(reg) _,
            t4 = out(reg) _,
            t5 = out(reg) _,
            t6 = out(reg) _,
            t7 = out(reg) _,
            // tmp = out(reg) tmp,
            out("rdx") _,
            options(pure, readonly, nostack)
        );
    }

    // println!("Tmp = {}", tmp);

    [r0, r1]
}


// #[allow(clippy::too_many_lines)]
// #[inline(always)]
// #[cfg(target_arch = "x86_64")]
// // #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
// pub(crate) fn mont_square_with_reduction_impl(a: &[u64; 4]) -> [u64; 4] {
//     // this is CIOS multiplication when top bit for top word of modulus is not set

//     let mut r0: u64;
//     let mut r1: u64;
//     let mut r2: u64;
//     let mut r3: u64;

//     // let mut tmp: u64;

//     unsafe {
//         asm!(
//             // round 0
//             "mov rdx, qword ptr [{a_ptr} + 0]",
//             "xor r8d, r8d",
//             "mulx r10, r9, qword ptr [{a_ptr} + 8]",
//             "mulx r15, r8, qword ptr [{a_ptr} + 16]",
//             "mulx r12, r11, qword ptr [{a_ptr} + 24]",
//             "adox r10, r8", 
//             "adcx r11, r15", 
//             "mov rdx, qword ptr [{a_ptr} + 8]",
//             "mulx r15, r8, qword ptr [{a_ptr} + 16]",
//             "mulx rcx, rdi, qword ptr [{a_ptr} + 24]",
//             "mov rdx, qword ptr [{a_ptr} + 16]",
//             "mulx r14, r13, qword ptr [{a_ptr} + 24]",
//             "adox r11, r8", 
//             "adcx r12, rdi", 
//             "adox r12, r15", 
//             "adcx r13, rcx",
//             "mov r8, 0",
//             "adox r13, r8",
//             "adcx r14, r8", 

//             // double
//             "adox r9, r9",
//             "adcx r12, r12", 
//             "adox r10, r10",
//             "adcx r13, r13", 
//             "adox r11, r11",
//             "adcx r14, r14", 

//             // square contributions
//             "mov rdx, qword ptr [{a_ptr} + 0]",
//             "mulx rcx, r8, rdx",
//             "mov rdx, qword ptr [{a_ptr} + 16]",
//             "mulx rdi, rdx, rdx",
//             "adox r12, rdx",
//             "adcx r9, rcx",
//             "adox r13, rdi",
//             "mov rdx, qword ptr [{a_ptr} + 24]",
//             "mulx r15, rcx, rdx",
//             "mov rdx, qword ptr [{a_ptr} + 8]",
//             "mulx rdx, rdi, rdx",
//             "adcx r10, rdi",
//             "adox r14, rcx",
//             "mov rdi, 0",
//             "adcx r11, rdx",
//             "adox r15, rdi", // finish the addition chain
//             "adcx r12, rdi",

//             // reduction round 0
//             "mov rdx, r8",
//             "neg rdx",
//             "xor rcx, rcx", // rcx = 0
//             "adox r8, rdx",
//             "mulx rdi, r8, qword ptr [rip + {q3_ptr}]", 
//             "adcx r12, rdi",
//             "adox r9, rcx",
//             "adcx r13, rcx",
//             "adox r10, rcx",
//             "adcx r9, rcx",
//             "adox r11, r8",
//             "adcx r10, rcx",
//             "adcx r11, rcx",
//             "adox r12, rcx",
//             "adcx r12, rcx",

//             // reduction round 1
//             "mov rdx, r9",
//             "neg rdx",
//             "xor rdi, rdi", // rdi = 0
//             "adox r12, rdi",
//             "mulx rcx, r8, qword ptr [rip + {q3_ptr}]", 
//             "adcx r12, r8",
//             "adox r13, rcx",
//             "adcx r13, rdi",
//             "adox r14, rdi",
//             "adcx r9, rdx",
//             "adox r10, rdi",
//             "adcx r10, rdi",
//             "adox r11, rdi",
//             "adcx r11, rdi",

//             // reduction round 2
//             "mov rdx, r10",
//             "neg rdx",
//             "mulx r9, r8, qword ptr [rip + {q3_ptr}]", 
//             "xor rcx, rcx", // rcx = 0
//             "adox r12, rcx",
//             "adcx r12, rcx",
//             "adox r13, rcx",
//             "adcx r13, r8",
//             "adox r14, r9",
//             "adcx r14, rcx",
//             "adox r15, rcx",
//             "adcx r10, rdx",
//             "adox r11, rcx",
//             "adcx r11, rcx",
//             "adox r12, rcx",

//             // reduction round 3
//             "mov rdx, r11",
//             "neg rdx",
//             "xor rcx, rcx", // rcx = 0
//             "adox r11, rdx",
//             "mulx r11, r10, qword ptr [rip + {q3_ptr}]", 
//             "adcx r12, rcx",
//             "adox r12, rcx",
//             "adcx r13, rcx",
//             "adox r13, rcx",
//             "adcx r14, r10",
//             "adox r14, rcx",
//             "adcx r15, r11",
//             "adox r15, rcx",
            
//             // reduction. We use sub/sbb
//             "mov r8, r12",
//             "sub r8, 1",
//             "mov r9, r13",
//             "sbb r9, 0",
//             "mov r10, r14",
//             "sbb r10, 0",
//             "mov r11, r15",
//             "mov rdx, qword ptr [rip + {q3_ptr}]",
//             "sbb r11, rdx",

//             // if CF == 1 then original result was ok (reduction wa not necessary)
//             // so if not carry (CMOVNQ) then we copy 
//             "cmovnc r12, r8",
//             "cmovnc r13, r9",
//             "cmovnc r14, r10",
//             "cmovnc r15, r11",  
//             // end of reduction
//             // inv = const 0xffffffffffffffffu64,
//             q3_ptr = sym MODULUS_3_STATIC,
//             a_ptr = in(reg) a.as_ptr(),
//             out("rcx") _, 
//             out("rdx") _, 
//             out("rdi") _, 
//             out("r8") _, 
//             out("r9") _, 
//             out("r10") _, 
//             out("r11") _, 
//             out("r12") r0, 
//             out("r13") r1, 
//             out("r14") r2, 
//             out("r15") r3,
//             options(pure, readonly, nostack)
//         );
//     }

//     // println!("Tmp = {}", tmp);

//     [r0, r1, r2, r3]
// }