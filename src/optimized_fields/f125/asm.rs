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

#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
pub(crate) fn sub_noborrow_impl(a: [u64; 2], b: [u64; 2]) -> [u64; 2] {
    let [mut r0, mut r1] = a;
    let [b0, b1] = b;

    unsafe {
        asm!(
            "sub {in_out_0}, {b_0}",
            "sbb {in_out_1}, {b_1}",
            in_out_0 = inlateout(reg) r0,
            in_out_1 = inlateout(reg) r1,
            b_0 = in(reg) b0,
            b_1 = in(reg) b1,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1]
}


#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
pub(crate) fn double_nocarry_impl(a: [u64; 2]) -> [u64; 2] {
    // we use ADD/ADC only here as it's the same latency, but shorter encoding
    let [mut r0, mut r1] = a;

    unsafe {
        asm!(
            "add {in_out_0}, {in_out_0}",
            "adc {in_out_1}, {in_out_1}",
            in_out_0 = inlateout(reg) r0,
            in_out_1 = inlateout(reg) r1,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1]
}

#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
pub(crate) fn add_with_reduction_impl(a: [u64; 2], b: [u64; 2]) -> [u64; 2] {
    let [mut r0, mut r1] = a;
    let [b0, b1] = b;

    unsafe {
        asm!(
            "adox {in_out_0}, {b_0}",
            "mov {t0}, {in_out_0}",
            "adcx {in_out_0}, qword ptr [rip + {q0_ptr}]",
            "adox {in_out_1}, {b_1}",
            "mov {t1}, {in_out_1}",
            "adcx {in_out_1}, qword ptr [rip + {q1_ptr}]",
            "cmovnc {in_out_0}, {t0}",
            "cmovnc {in_out_1}, {t1}",
            q0_ptr = sym MODULUS_NEGATED_STATIC_0,
            q1_ptr = sym MODULUS_NEGATED_STATIC_1,
            in_out_0 = inlateout(reg) r0,
            in_out_1 = inlateout(reg) r1,
            b_0 = in(reg) b0,
            b_1 = in(reg) b1,
            t0 = out(reg) _,
            t1 = out(reg) _,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1]
}

#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
pub(crate) fn double_with_reduction_impl(a: [u64; 2]) -> [u64; 2] {
    let [mut r0, mut r1] = a;

    unsafe {
        asm!(
            "adox {in_out_0}, {in_out_0}",
            "mov {t0}, {in_out_0}",
            "adcx {in_out_0}, qword ptr [rip + {q0_ptr}]",
            "adox {in_out_1}, {in_out_1}",
            "mov {t1}, {in_out_1}",
            "adcx {in_out_1}, qword ptr [rip + {q1_ptr}]",
            "cmovnc {in_out_0}, {t0}",
            "cmovnc {in_out_1}, {t1}",
            q0_ptr = sym MODULUS_NEGATED_STATIC_0,
            q1_ptr = sym MODULUS_NEGATED_STATIC_1,
            in_out_0 = inlateout(reg) r0,
            in_out_1 = inlateout(reg) r1,
            t0 = out(reg) _,
            t1 = out(reg) _,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1]
}

#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
pub(crate) fn sub_with_reduction_impl(a: [u64; 2], b: [u64; 2]) -> [u64; 2] {
    let [mut r0, mut r1] = a;
    let [b0, b1] = b;

    unsafe {
        asm!(
            "sub {in_out_0}, {b_0}",
            "sbb {in_out_1}, {b_1}",
            "mov {t0}, {in_out_0}",
            "add {in_out_0}, 1",
            "mov {t1}, {in_out_1}",
            "adc {in_out_1}, qword ptr [rip + {q1_ptr}]",
            "cmovnc {in_out_0}, {t0}",
            "cmovnc {in_out_1}, {t1}",
            q1_ptr = sym MODULUS_1_STATIC,
            in_out_0 = inlateout(reg) r0,
            in_out_1 = inlateout(reg) r1,
            b_0 = in(reg) b0,
            b_1 = in(reg) b1,
            t0 = out(reg) _,
            t1 = out(reg) _,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1]
}

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


#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
pub(crate) fn reduce_by_modulus_impl(a: [u64; 2]) -> [u64; 2] {
    let [mut r0, mut r1] = a;
    // modulus : higher-limb: 0x3000000300000001  lower-limb: 1
    // reduce number when it is greater than modulus
    // in order to check this test
    // we substract 1 from lower limb.
    // if it borrows(CF) then it means number is higher than modulus
    // we substract higer limb of modulus from higher limb of number
    // if it generates and CF then it means number is higher than modulus
    // and reduction will happen
    // same for lower limb.
    unsafe {
        asm!(
            "mov {t0}, {in_out_0}",
            "sub {in_out_0}, 1",
            "mov {t1}, {in_out_1}",
            "sbb {in_out_1}, {q1_const}",
            // "sbb {in_out_1}, qword ptr [rip + {q1_ptr}]",
            "cmovc {in_out_0}, {t0}",
            "cmovc {in_out_1}, {t1}",
            q1_const = const MODULUS_LIMB_1,
            // q1_ptr = sym MODULUS_1_STATIC,
            in_out_0 = inlateout(reg) r0,
            in_out_1 = inlateout(reg) r1,
            t0 = out(reg) _,
            t1 = out(reg) _,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1]
}

#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
pub(crate) fn mont_mul_with_reduction_impl(a: [u64; 2], b: [u64; 2]) -> [u64; 2] {
    // we also manually use the facts about INV = -1 and modulus[0] = 1

    let [mut r0, mut r1] = a;
    let [b0, b1] = b;

    // round 0
    // |    |    | t1 | t0 |
    // |    | t3 | t2 |    |
    // reduction 0
    // |    |    | 1  |  0 |
    // |    | i0 | t0 |    |
    // round 1
    // |    | t3 | t1 |    |
    // |    | t2 | t0 |    |
    // | i1 | i0 |    |    |
    // reduction 1
    // |    | 1  | 0  |    |
    // | i1 | i0 |    |    |

    unsafe {
        asm!(
            "mov rdx, {in_out_0}",
            "mulx {t1}, {t0}, {b_0}",
            // a0 * b1
            "mulx {t3}, {t2}, {b_1}",
            // my now we wait to get m = -t0
            // instead of later adding t0 and -t0, we just do set CF
            "mov rdx, {t0}", // {t0} is reusable
            "neg rdx",
            "xor {t0}, {t0}",
            "mulx {in_out_0}, {t0}, qword ptr [rip + {q1_ptr}]", // second part of first reduction
            "adox {t1}, {t2}", // {t2} is now reusable,
            "adox {t3}, qword ptr [rip + {zero_ptr}]", // complete the carry flag
            // OF = 0
            "stc",
            "adcx {t1}, {t0}",
            "adcx {t3}, {in_out_0}",
            // OF = CF = 0

            // round 1
            "mov rdx, {in_out_1}",
            "mulx {t2}, {t0}, {b_0}",
            "mulx {in_out_1}, {in_out_0}, {b_1}",
            "adox {t1}, {t0}", // r0
            "adox {t2}, qword ptr [rip + {zero_ptr}]", // carry
            // OF = 0
            "mov rdx, {t1}",
            "neg rdx",
            "xor {t0}, {t0}",
            // "xor rdx, qword ptr [rip + {u64_max_ptr}]",
            "mulx {t1}, {t0}, qword ptr [rip + {q1_ptr}]", // second part of first reduction
            "stc",
            "adcx {in_out_0}, {t2}", // carry + low
            "adcx {in_out_1}, qword ptr [rip + {zero_ptr}]", // propagate
            "adox {in_out_0}, {t3}", // r1
            "adox {in_out_1}, qword ptr [rip + {zero_ptr}]", // carry
            // m * q1
            "adcx {in_out_0}, {t0}", // t0 + red_carry + r1
            "adcx {in_out_1}, {t1}", // t1 + propagate + carry

            "mov {t0}, {in_out_0}",
            "sub {in_out_0}, 1",
            "mov {t1}, {in_out_1}",
            "sbb {in_out_1}, qword ptr [rip + {q1_ptr}]",

            // if CF == 0 then original result was ok (reduction wa not necessary)
            "cmovc {in_out_0}, {t0}",
            "cmovc {in_out_1}, {t1}",
            // end of reduction
            zero_ptr = sym ZERO_U64,
            // u64_max_ptr = sym U64_MAX,
            q1_ptr = sym MODULUS_1_STATIC,
            in_out_0 = inlateout(reg) r0,
            in_out_1 = inlateout(reg) r1,
            b_0 = in(reg) b0,
            b_1 = in(reg) b1,
            t0 = out(reg) _,
            t1 = out(reg) _,
            t2 = out(reg) _,
            t3 = out(reg) _,
            out("rdx") _,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1]
}


// TODO: not yet valid, need to take care when r0/r1 == 0, so we do not 
// need to "stc"

#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
pub(crate) fn mont_mul_with_partial_reduction_impl(a: [u64; 2], b: [u64; 2]) -> [u64; 2] {
    // we also manually use the facts about INV = -1 and modulus[0] = 1

    let [mut r0, mut r1] = a;
    let [b0, b1] = b;

    // round 0
    // |    |    | t1 | t0 |
    // |    | t3 | t2 |    |
    // reduction 0
    // |    |    | 1  |  0 |
    // |    | i0 | t0 |    |
    // round 1
    // |    | t3 | t1 |    |
    // |    | t2 | t0 |    |
    // | i1 | i0 |    |    |
    // reduction 1
    // |    | 1  | 0  |    |
    // | i1 | i0 |    |    |

    unsafe {
        asm!(
            "mov rdx, {in_out_0}",
            "mulx {t1}, {t0}, {b_0}",
            // a0 * b1
            "mulx {t3}, {t2}, {b_1}",
            // my now we wait to get m = -t0
            // instead of later adding t0 and -t0, we just do set CF
            "mov rdx, {t0}", // {t0} is reusable
            "neg rdx",
            "xor {t0}, {t0}",
            "mulx {in_out_0}, {t0}, qword ptr [rip + {q1_ptr}]", // second part of first reduction
            "adox {t1}, {t2}", // {t2} is now reusable,
            "adox {t3}, qword ptr [rip + {zero_ptr}]", // complete the carry flag
            // OF = 0
            "stc",
            "adcx {t1}, {t0}",
            "adcx {t3}, {in_out_0}",
            // OF = CF = 0

            // round 1
            "mov rdx, {in_out_1}",
            "mulx {t2}, {t0}, {b_0}",
            "mulx {in_out_1}, {in_out_0}, {b_1}",
            "adox {t1}, {t0}", // r0
            "adox {t2}, qword ptr [rip + {zero_ptr}]", // carry
            // OF = 0
            "mov rdx, {t1}",
            "neg rdx",
            "xor {t0}, {t0}",
            // "xor rdx, qword ptr [rip + {u64_max_ptr}]",
            "mulx {t1}, {t0}, qword ptr [rip + {q1_ptr}]", // second part of first reduction
            "stc",
            "adcx {in_out_0}, {t2}", // carry + low
            "adcx {in_out_1}, qword ptr [rip + {zero_ptr}]", // propagate
            "adox {in_out_0}, {t3}", // r1
            "adox {in_out_1}, qword ptr [rip + {zero_ptr}]", // carry
            // m * q1
            "adcx {in_out_0}, {t0}", // t0 + red_carry + r1
            "adcx {in_out_1}, {t1}", // t1 + propagate + carry
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
            out("rdx") _,
            options(pure, readonly, nostack)
        );
    }

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