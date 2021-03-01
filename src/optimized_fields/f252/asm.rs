use super::*;

#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
pub(crate) fn add_nocarry_impl(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
    // we use ADD/ADC only here as it's the same latency, but shorter encoding
    
    let mut r0: u64;
    let mut r1: u64;
    let mut r2: u64;
    let mut r3: u64;

    unsafe {
        asm!(
            "xor r12d, r12d",
            "mov r12, qword ptr [{a_ptr} + 0]",
            "add r12, qword ptr [{b_ptr} + 0]",
            "mov r13, qword ptr [{a_ptr} + 8]",
            "adc r13, qword ptr [{b_ptr} + 8]",
            "mov r14, qword ptr [{a_ptr} + 16]",
            "adc r14, qword ptr [{b_ptr} + 16]",
            "mov r15, qword ptr [{a_ptr} + 24]",
            "adc r15, qword ptr [{b_ptr} + 24]",
            a_ptr = in(reg) a.as_ptr(),
            b_ptr = in(reg) b.as_ptr(),
            out("r12") r0, 
            out("r13") r1, 
            out("r14") r2, 
            out("r15") r3,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1, r2, r3]
}

#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
pub(crate) fn sub_noborrow_impl(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
    let mut r0: u64;
    let mut r1: u64;
    let mut r2: u64;
    let mut r3: u64;

    unsafe {
        asm!(
            "xor r12d, r12d",
            "mov r12, qword ptr [{a_ptr} + 0]",
            "sub r12, qword ptr [{b_ptr} + 0]",
            "mov r13, qword ptr [{a_ptr} + 8]",
            "sbb r13, qword ptr [{b_ptr} + 8]",
            "mov r14, qword ptr [{a_ptr} + 16]",
            "sbb r14, qword ptr [{b_ptr} + 16]",
            "mov r15, qword ptr [{a_ptr} + 24]",
            "sbb r15, qword ptr [{b_ptr} + 24]",
            a_ptr = in(reg) a.as_ptr(),
            b_ptr = in(reg) b.as_ptr(),
            out("r12") r0, 
            out("r13") r1, 
            out("r14") r2, 
            out("r15") r3,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1, r2, r3]
}


#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
pub(crate) fn double_nocarry_impl(a: &[u64; 4]) -> [u64; 4] {
    // we use ADD/ADC only here as it's the same latency, but shorter encoding
    
    let mut r0: u64;
    let mut r1: u64;
    let mut r2: u64;
    let mut r3: u64;

    unsafe {
        asm!(
            "xor r12d, r12d",
            "mov r12, qword ptr [{a_ptr} + 0]",
            "add r12, r12",
            "mov r13, qword ptr [{a_ptr} + 8]",
            "adc r13, r13",
            "mov r14, qword ptr [{a_ptr} + 16]",
            "adc r14, r14",
            "mov r15, qword ptr [{a_ptr} + 24]",
            "adc r15, r15",
            a_ptr = in(reg) a.as_ptr(),
            out("r12") r0, 
            out("r13") r1, 
            out("r14") r2, 
            out("r15") r3,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1, r2, r3]
}

#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
pub(crate) fn add_with_reduction_impl(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
    let mut r0: u64;
    let mut r1: u64;
    let mut r2: u64;
    let mut r3: u64;

    unsafe {
        asm!(
            // we sum (a+b) using addition chain with OF
            // and sum (a+b) - p using addition chain with CF
            // if (a+b) does not overflow the modulus
            // then sum (a+b) - p will produce CF
            "xor r12d, r12d",
            "mov r12, qword ptr [{a_ptr} + 0]",
            "mov r13, qword ptr [{a_ptr} + 8]",
            "mov r14, qword ptr [{a_ptr} + 16]",
            "mov r15, qword ptr [{a_ptr} + 24]",
            "adox r12, qword ptr [{b_ptr} + 0]",
            "mov r8, r12",
            "adcx r8, qword ptr [rip + {q0_ptr}]",
            "adox r13, qword ptr [{b_ptr} + 8]",
            "mov r9, r13",
            "adcx r9, qword ptr [rip + {q1_ptr}]",
            "adox r14, qword ptr [{b_ptr} + 16]",
            "mov r10, r14",
            "adcx r10, qword ptr [rip + {q2_ptr}]",
            "adox r15, qword ptr [{b_ptr} + 24]",
            "mov r11, r15",
            "adcx r11, qword ptr [rip + {q3_ptr}]",

            // if CF = 0 then take value (a+b) from [r12, .., r15]
            // otherwise take (a+b) - p

            "cmovc r12, r8",
            "cmovc r13, r9",
            "cmovc r14, r10",
            "cmovc r15, r11",  

            q0_ptr = sym MODULUS_NEGATED_STATIC_0,
            q1_ptr = sym MODULUS_NEGATED_STATIC_1,
            q2_ptr = sym MODULUS_NEGATED_STATIC_2,
            q3_ptr = sym MODULUS_NEGATED_STATIC_3,
            // end of reduction
            a_ptr = in(reg) a.as_ptr(),
            b_ptr = in(reg) b.as_ptr(),
            out("r8") _, 
            out("r9") _, 
            out("r10") _, 
            out("r11") _, 
            out("r12") r0, 
            out("r13") r1, 
            out("r14") r2, 
            out("r15") r3,
            options(pure, readonly, nostack)
        );
    }

    // unsafe {
    //     asm!(
    //         "xor r12d, r12d",
    //         "mov r12, qword ptr [{a_ptr} + 0]",
    //         "mov r13, qword ptr [{a_ptr} + 8]",
    //         "mov r14, qword ptr [{a_ptr} + 16]",
    //         "mov r15, qword ptr [{a_ptr} + 24]",
    //         "add r12, qword ptr [{b_ptr} + 0]",
    //         "adc r13, qword ptr [{b_ptr} + 8]",
    //         "adc r14, qword ptr [{b_ptr} + 16]",
    //         "adc r15, qword ptr [{b_ptr} + 24]",

    //         "mov r8, r12",
    //         "mov rdx, qword ptr [rip + {q0_ptr}]",
    //         "sub r8, rdx",
    //         "mov r9, r13",
    //         "mov rdx, qword ptr [rip + {q1_ptr}]",
    //         "sbb r9, rdx",
    //         "mov r10, r14",
    //         "mov rdx, qword ptr [rip + {q2_ptr}]",
    //         "sbb r10, rdx",
    //         "mov r11, r15",
    //         "mov rdx, qword ptr [rip + {q3_ptr}]",
    //         "sbb r11, rdx",

    //         "cmovnc r12, r8",
    //         "cmovnc r13, r9",
    //         "cmovnc r14, r10",
    //         "cmovnc r15, r11",  

    //         q0_ptr = sym #m0,
    //         q1_ptr = sym #m1,
    //         q2_ptr = sym #m2,
    //         q3_ptr = sym #m3,
    //         // end of reduction
    //         a_ptr = in(reg) a.as_ptr(),
    //         b_ptr = in(reg) b.as_ptr(),
    //         out("rdx") _, 
    //         out("r8") _, 
    //         out("r9") _, 
    //         out("r10") _, 
    //         out("r11") _, 
    //         out("r12") r0, 
    //         out("r13") r1, 
    //         out("r14") r2, 
    //         out("r15") r3,
    //         options(pure, readonly, nostack)
    //     );
    // }

    [r0, r1, r2, r3]
}

#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
pub(crate) fn double_with_reduction_impl(a: &[u64; 4]) -> [u64; 4] {
    let mut r0: u64;
    let mut r1: u64;
    let mut r2: u64;
    let mut r3: u64;

    unsafe {
        asm!(
            // we sum (a+b) using addition chain with OF
            // and sum (a+b) - p using addition chain with CF
            // if (a+b) does not overflow the modulus
            // then sum (a+b) will produce CF
            "xor r12d, r12d",
            "mov r12, qword ptr [{a_ptr} + 0]",
            "mov r13, qword ptr [{a_ptr} + 8]",
            "mov r14, qword ptr [{a_ptr} + 16]",
            "mov r15, qword ptr [{a_ptr} + 24]",
            "adox r12, r12",
            "mov r8, r12",
            "adcx r8, qword ptr [rip + {q0_ptr}]",
            "adox r13, r13",
            "mov r9, r13",
            "adcx r9, qword ptr [rip + {q1_ptr}]",
            "adox r14, r14",
            "mov r10, r14",
            "adcx r10, qword ptr [rip + {q2_ptr}]",
            "adox r15, r15",
            "mov r11, r15",
            "adcx r11, qword ptr [rip + {q3_ptr}]",

            // if CF = 0 then take value (a+b) from [r12, .., r15]
            // otherwise take (a+b) - p

            "cmovc r12, r8",
            "cmovc r13, r9",
            "cmovc r14, r10",
            "cmovc r15, r11",  

            q0_ptr = sym MODULUS_NEGATED_STATIC_0,
            q1_ptr = sym MODULUS_NEGATED_STATIC_1,
            q2_ptr = sym MODULUS_NEGATED_STATIC_2,
            q3_ptr = sym MODULUS_NEGATED_STATIC_3,
            // end of reduction
            a_ptr = in(reg) a.as_ptr(),
            out("r8") _, 
            out("r9") _, 
            out("r10") _, 
            out("r11") _, 
            out("r12") r0, 
            out("r13") r1, 
            out("r14") r2, 
            out("r15") r3,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1, r2, r3]
}

#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
pub(crate) fn sub_with_reduction_impl(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
    let mut r0: u64;
    let mut r1: u64;
    let mut r2: u64;
    let mut r3: u64;

    unsafe {
        asm!(
            "xor r12d, r12d",
            "mov r12, qword ptr [{a_ptr} + 0]",
            "sub r12, qword ptr [{b_ptr} + 0]",
            "mov r13, qword ptr [{a_ptr} + 8]",
            "sbb r13, qword ptr [{b_ptr} + 8]",
            "mov r14, qword ptr [{a_ptr} + 16]",
            "sbb r14, qword ptr [{b_ptr} + 16]",
            "mov r15, qword ptr [{a_ptr} + 24]",
            "sbb r15, qword ptr [{b_ptr} + 24]",

            // duplicate (a-b) into [r8, r9, r10, r11]

            // now make [r12, .., r15] + modulus;
            // if (a-b) did underflow then (a-b) + modulus > 2^256
            // so below we will get an overflow

            "mov r8, r12",
            "add r12, 1",
            // "adox r12, qword ptr [rip + {q0_ptr}]",
            "mov r9, r13",
            "adc r13, 0",
            // "adox r13, qword ptr [rip + {q1_ptr}]",
            "mov r10, r14",
            "adc r14, 0",
            // "adox r14, qword ptr [rip + {q2_ptr}]",
            "mov r11, r15",
            "adc r15, qword ptr [rip + {q3_ptr}]",
            // "adox r15, qword ptr [rip + {q3_ptr}]",

            // if no carry (no overflow) then what we had in [r8, r11] is valid
            "cmovnc r12, r8",
            "cmovnc r13, r9",
            "cmovnc r14, r10",
            "cmovnc r15, r11", 
            // end of reduction
            q3_ptr = sym MODULUS_3_STATIC,
            a_ptr = in(reg) a.as_ptr(),
            b_ptr = in(reg) b.as_ptr(),
            out("r8") _, 
            out("r9") _, 
            out("r10") _, 
            out("r11") _, 
            out("r12") r0, 
            out("r13") r1, 
            out("r14") r2, 
            out("r15") r3,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1, r2, r3]
}

#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
pub(crate) fn add_modulus_impl(a: &[u64; 4]) -> [u64; 4] {
    // we use ADD/ADC only here as it's the same latency, but shorter encoding
    
    let mut r0: u64;
    let mut r1: u64;
    let mut r2: u64;
    let mut r3: u64;

    unsafe {
        asm!(
            "xor r12d, r12d",
            "mov r12, qword ptr [{a_ptr} + 0]",
            "add r12, 1",
            "mov r13, qword ptr [{a_ptr} + 8]",
            "adc r13, 0",
            "mov r14, qword ptr [{a_ptr} + 16]",
            "adc r14, 0",
            "mov r15, qword ptr [{a_ptr} + 24]",
            "adc r15, qword ptr [rip + {q3_ptr}]",
            a_ptr = in(reg) a.as_ptr(),
            q3_ptr = sym MODULUS_3_STATIC,
            out("r12") r0, 
            out("r13") r1, 
            out("r14") r2, 
            out("r15") r3,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1, r2, r3]
}


#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
pub(crate) fn sub_modulus_impl(a: &[u64; 4]) -> [u64; 4] {
    // we use SUB/SBB to utilize 0s
    
    let mut r0: u64;
    let mut r1: u64;
    let mut r2: u64;
    let mut r3: u64;

    unsafe {
        asm!(
            "xor r12d, r12d",
            "mov r12, qword ptr [{a_ptr} + 0]",
            "sub r12, 1",
            "mov r13, qword ptr [{a_ptr} + 8]",
            "sub r13, 0",
            "mov r14, qword ptr [{a_ptr} + 16]",
            "sub r14, 0",
            "mov r15, qword ptr [{a_ptr} + 24]",
            "sbb r15, qword ptr [rip + {q3_ptr}]",
            a_ptr = in(reg) a.as_ptr(),
            q3_ptr = sym MODULUS_3_STATIC,
            out("r12") r0, 
            out("r13") r1, 
            out("r14") r2, 
            out("r15") r3,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1, r2, r3]
}


#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
pub(crate) fn reduce_by_modulus_impl(a: &[u64; 4]) -> [u64; 4] {
    let mut r0: u64;
    let mut r1: u64;
    let mut r2: u64;
    let mut r3: u64;

    unsafe {
        asm!(
            "xor r12d, r12d",
            "mov r12, qword ptr [{a_ptr} + 0]",
            "mov r8, r12",
            "sub r12, 1",
            "mov r13, qword ptr [{a_ptr} + 8]",
            "mov r9, r13",
            "sbb r13, 0",
            "mov r14, qword ptr [{a_ptr} + 16]",
            "mov r10, r14",
            "sbb r14, 0",
            "mov r15, qword ptr [{a_ptr} + 24]",
            "mov r11, r15",
            "sbb r15, qword ptr [rip + {q3_ptr}]",

            "cmovc r12, r8",
            "cmovc r13, r9",
            "cmovc r14, r10",
            "cmovc r15, r11", 

            q3_ptr = sym MODULUS_3_STATIC,
            a_ptr = in(reg) a.as_ptr(),
            out("r8") _, 
            out("r9") _, 
            out("r10") _, 
            out("r11") _, 
            out("r12") r0, 
            out("r13") r1, 
            out("r14") r2, 
            out("r15") r3,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1, r2, r3]
}

#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
pub(crate) fn mont_mul_with_reduction_impl(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
    // we also manually use the facts about INV = -1 and modulus[0] = 1, modulus[1] = modulus[2] = 0

    let mut r0: u64;
    let mut r1: u64;
    let mut r2: u64;
    let mut r3: u64;

    unsafe {
        asm!(
            // round 0
            "mov rdx, qword ptr [{a_ptr} + 0]",
            "xor r8d, r8d",
            "mulx r14, r13, qword ptr [{b_ptr} + 0]",
            "mulx r9, r8, qword ptr [{b_ptr} + 8]",
            "mulx r10, r15, qword ptr [{b_ptr} + 16]",
            "mulx r12, rdi, qword ptr [{b_ptr} + 24]",
            "mov rdx, r13", 
            // here we multiply rdx * inv -> k0 (we ignore r11), so just negate
            "neg rdx",
            "xor r11d, r11d", // clear the flags, and from now on r11 == 0
            "adcx r14, r8", 
            "adox r10, rdi", 
            // start calculating immediately
            "mulx r8, rdi, qword ptr [rip + {q3_ptr}]", 
            "adcx r15, r9",
            "adox r12, r11",
            "adcx r10, r11",
            "adox r13, rdx", 
            "adcx r14, r11",
            "adox r14, r11", 
            "adcx r15, r11", 
            "adox r15, r11", 
            "adcx r10, rdi",
            "adox r10, r11", 
            "adcx r12, r8", 
            "adox r12, r11",

            // round 1
            "mov rdx, qword ptr [{a_ptr} + 8]",
            "mulx r9, r8, qword ptr [{b_ptr} + 0]",
            "mulx r11, rdi, qword ptr [{b_ptr} + 8]",
            "adcx r14, r8",
            "adox r15, r9",
            "mulx r9, r8, qword ptr [{b_ptr} + 16]",
            "adcx r15, rdi",
            "adox r10, r11",
            "mulx r13, rdi, qword ptr [{b_ptr} + 24]",
            "adcx r10, r8",
            "adox r12, rdi",
            "adcx r12, r9",
            "mov rdi, 0",
            "adox r13, rdi",
            "adcx r13, rdi",
            "mov rdx, r14",
            "neg rdx",
            "mulx r11, rdi, qword ptr [rip + {q3_ptr}]",
            "xor r8d, r8d", // clear the flags!, r8 == 0
            "adox r14, rdx",
            "adcx r15, r8",
            "adox r15, r8",
            "adcx r10, r8",
            "adox r10, r8",
            "adcx r12, r8",
            "adox r12, rdi",
            "adcx r13, r11",
            "adox r13, r8",

            // round 2
            "mov rdx, qword ptr [{a_ptr} + 16]",
            "mulx r9, r8, qword ptr [{b_ptr} + 0]",
            "mulx r11, rdi, qword ptr [{b_ptr} + 8]",
            "adcx r15, r8",
            "adox r10, r9",
            "mulx r9, r8, qword ptr [{b_ptr} + 16]",
            "adcx r10, rdi",
            "adox r12, r11",
            "mulx r14, rdi, qword ptr [{b_ptr} + 24]",
            "adcx r12, r8",
            "adox r13, r9",
            "adcx r13, rdi",
            "mov r9, 0",
            "adox r14, r9",
            "adcx r14, r9",
            "mov rdx, r15",
            "neg rdx",
            "mulx r11, rdi, qword ptr [rip + {q3_ptr}]",
            "xor r8d, r8d", // clear the flags!, r8 == 0
            "adox r15, rdx",
            "adcx r10, r8",
            "adox r10, r8",
            "adcx r12, r8",
            "adox r12, r8",
            "adcx r13, r8",
            "adox r13, rdi",
            "adcx r14, r11",
            "adox r14, r8",

            // round 3
            "mov rdx, qword ptr [{a_ptr} + 24]",
            "mulx r9, r8, qword ptr [{b_ptr} + 0]",
            "mulx r11, rdi, qword ptr [{b_ptr} + 8]",
            "adcx r10, r8",
            "adox r12, r9",
            "mulx r9, r8, qword ptr [{b_ptr} + 16]",
            "adcx r12, rdi",
            "adox r13, r11",
            "mulx r15, rdi, qword ptr [{b_ptr} + 24]",
            "adcx r13, r8",
            "adox r14, r9",
            "adcx r14, rdi",
            "mov r9, 0",
            "adox r15, r9",
            "adcx r15, r9",
            "mov rdx, r10",
            "neg rdx",
            "mulx r9, rdi, qword ptr [rip + {q3_ptr}]",
            "xor r8d, r8d", // clear the flags! r8 == 0
            "adox r10, rdx",
            "adcx r12, r8",
            "adox r12, r8",
            "adcx r13, r8",
            "adox r13, r8",
            "adcx r14, r8",
            "adox r14, rdi",
            "adcx r15, r9",
            "adox r15, r8",

            "mov r8, r12",
            "sub r8, 1",
            "mov r9, r13",
            "sbb r9, 0",
            "mov r10, r14",
            "sbb r10, 0",
            "mov r11, r15",
            "mov rdx, qword ptr [rip + {q3_ptr}]",
            "sbb r11, rdx",

            // if CF == 1 then original result was ok (reduction wa not necessary)
            // so if not carry (CMOVNQ) then we copy 
            "cmovnc r12, r8",
            "cmovnc r13, r9",
            "cmovnc r14, r10",
            "cmovnc r15, r11",  
            // end of reduction
            q3_ptr = sym MODULUS_3_STATIC,
            a_ptr = in(reg) a.as_ptr(),
            b_ptr = in(reg) b.as_ptr(),
            out("rdx") _, 
            out("rdi") _, 
            out("r8") _, 
            out("r9") _, 
            out("r10") _, 
            out("r11") _, 
            out("r12") r0, 
            out("r13") r1, 
            out("r14") r2, 
            out("r15") r3,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1, r2, r3]
}

fn mont_mul_new(a: &[u64;4], b: &[u64;4]) -> [u64;4]{

    let mut r0: u64;
    let mut r1: u64;
    let mut r2: u64;
    let mut r3: u64;
    
    unsafe{
        asm!{
            "mov rdx, qword ptr[{a_ptr} + 0]", // a0 @ rdx
            "xor rax, rax",  // clear r10 register, we use this when we need 0 
            "mulx r14, r13, qword ptr [{b_ptr} + 0]", // a0 * b0
            "mulx r9, r8, qword ptr [{b_ptr} + 8]", // a0 * b1
            
            "mulx r8, rcx, qword ptr [{b_ptr} + 0]", // a0 * b0
            "mulx r9, rax, qword ptr [{b_ptr} + 8]", // a0 * b1

            // "adcx r8, rax", 
            // "mulx r10, rax, qword ptr[{b_ptr} + 16]", // a0 * b2
            // "adcx r9, rax", 
            // "mulx r11, rax, qword ptr[{b_ptr} + 24]", // a0 * b2
            // "adcx r10, rax", 
            
            "mov r12, r8",
            "mov r13, r9",
            "mov r14, r10",
            "mov r15, r11",
            // "mulx r14 r13 qword ptr[{b_ptr} + 0]",  // (r[0], r[1]) <- a[0] * b[0]
            // "mulx r9 r8 qword ptr[{b_ptr} + 8]",    // (t[0], t[1]) <- a[0] * b[1]
            // "mulx r15 r8 qword ptr[{b_ptr} + 16]",    // (r[2] , r[3]) <- a[0] * b[2]
            // "mulx r9 r8 qword ptr[{b_ptr} + 8]",    // (t[2], r[4]) <- a[0] * b[3] (overwrite a[0]) 
            a_ptr = in(reg) a.as_ptr(),
            b_ptr = in(reg) b.as_ptr(),
            out("rdx") _,
            out("rdi") _,
            out("r8") _,
            out("r9") _,
            out("r10") _,
            out("r11") _,
            out("r12") r0,
            out("r13") r1,
            out("r14") r2,
            out("r15") r3,
        }
    };

    [r0,r1,r2,r3]
}

#[test]
fn test_new_mont_mul(){
    // use rand::{XorShiftRng, SeedableRng, Rand};

    let a: [u64;4] = [1,2,3,4]; 
    let b: [u64;4] = [5,6,7,8]; 

    let c = mont_mul_new(&a, &b);

    println!("c {:?}", c);
}


#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
pub(crate) fn mont_mul_with_partial_reduction_impl(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
    // we also manually use the facts about INV = -1 and modulus[0] = 1, modulus[1] = modulus[2] = 0

    let mut r0: u64;
    let mut r1: u64;
    let mut r2: u64;
    let mut r3: u64;

    unsafe {
        asm!(
            // round 0
            "mov rdx, qword ptr [{a_ptr} + 0]",
            "xor r8d, r8d",
            "mulx r14, r13, qword ptr [{b_ptr} + 0]",
            "mulx r9, r8, qword ptr [{b_ptr} + 8]",
            "mulx r10, r15, qword ptr [{b_ptr} + 16]",
            "mulx r12, rdi, qword ptr [{b_ptr} + 24]",
            "mov rdx, r13", 
            // here we multiply rdx * inv -> k0 (we ignore r11), so just negate
            "neg rdx",
            "xor r11d, r11d", // clear the flags, and from now on r11 == 0
            "adcx r14, r8", 
            "adox r10, rdi", 
            // start calculating immediately
            "mulx r8, rdi, qword ptr [rip + {q3_ptr}]", 
            "adcx r15, r9",
            "adox r12, r11",
            "adcx r10, r11",
            "adox r13, rdx", 
            "adcx r14, r11",
            "adox r14, r11", 
            "adcx r15, r11", 
            "adox r15, r11", 
            "adcx r10, rdi",
            "adox r10, r11", 
            "adcx r12, r8", 
            "adox r12, r11",

            // round 1
            "mov rdx, qword ptr [{a_ptr} + 8]",
            "mulx r9, r8, qword ptr [{b_ptr} + 0]",
            "mulx r11, rdi, qword ptr [{b_ptr} + 8]",
            "adcx r14, r8",
            "adox r15, r9",
            "mulx r9, r8, qword ptr [{b_ptr} + 16]",
            "adcx r15, rdi",
            "adox r10, r11",
            "mulx r13, rdi, qword ptr [{b_ptr} + 24]",
            "adcx r10, r8",
            "adox r12, rdi",
            "adcx r12, r9",
            "mov rdi, 0",
            "adox r13, rdi",
            "adcx r13, rdi",
            "mov rdx, r14",
            "neg rdx",
            "mulx r11, rdi, qword ptr [rip + {q3_ptr}]",
            "xor r8d, r8d", // clear the flags!, r8 == 0
            "adox r14, rdx",
            "adcx r15, r8",
            "adox r15, r8",
            "adcx r10, r8",
            "adox r10, r8",
            "adcx r12, r8",
            "adox r12, rdi",
            "adcx r13, r11",
            "adox r13, r8",

            // round 2
            "mov rdx, qword ptr [{a_ptr} + 16]",
            "mulx r9, r8, qword ptr [{b_ptr} + 0]",
            "mulx r11, rdi, qword ptr [{b_ptr} + 8]",
            "adcx r15, r8",
            "adox r10, r9",
            "mulx r9, r8, qword ptr [{b_ptr} + 16]",
            "adcx r10, rdi",
            "adox r12, r11",
            "mulx r14, rdi, qword ptr [{b_ptr} + 24]",
            "adcx r12, r8",
            "adox r13, r9",
            "adcx r13, rdi",
            "mov r9, 0",
            "adox r14, r9",
            "adcx r14, r9",
            "mov rdx, r15",
            "neg rdx",
            "mulx r11, rdi, qword ptr [rip + {q3_ptr}]",
            "xor r8d, r8d", // clear the flags!, r8 == 0
            "adox r15, rdx",
            "adcx r10, r8",
            "adox r10, r8",
            "adcx r12, r8",
            "adox r12, r8",
            "adcx r13, r8",
            "adox r13, rdi",
            "adcx r14, r11",
            "adox r14, r8",

            // round 3
            "mov rdx, qword ptr [{a_ptr} + 24]",
            "mulx r9, r8, qword ptr [{b_ptr} + 0]",
            "mulx r11, rdi, qword ptr [{b_ptr} + 8]",
            "adcx r10, r8",
            "adox r12, r9",
            "mulx r9, r8, qword ptr [{b_ptr} + 16]",
            "adcx r12, rdi",
            "adox r13, r11",
            "mulx r15, rdi, qword ptr [{b_ptr} + 24]",
            "adcx r13, r8",
            "adox r14, r9",
            "adcx r14, rdi",
            "mov r9, 0",
            "adox r15, r9",
            "adcx r15, r9",
            "mov rdx, r10",
            "neg rdx",
            "mulx r9, rdi, qword ptr [rip + {q3_ptr}]",
            "xor r8d, r8d", // clear the flags! r8 == 0
            "adox r10, rdx",
            "adcx r12, r8",
            "adox r12, r8",
            "adcx r13, r8",
            "adox r13, r8",
            "adcx r14, r8",
            "adox r14, rdi",
            "adcx r15, r9",
            "adox r15, r8",

            "mov r8, r12",
            "sub r8, 1",
            "mov r9, r13",
            "sbb r9, 0",
            "mov r10, r14",
            "sbb r10, 0",
            "mov r11, r15",
            "mov rdx, qword ptr [rip + {q3_ptr}]",
            "sbb r11, rdx",

            // if CF == 1 then original result was ok (reduction wa not necessary)
            // so if not carry (CMOVNQ) then we copy 
            "cmovnc r12, r8",
            "cmovnc r13, r9",
            "cmovnc r14, r10",
            "cmovnc r15, r11",  
            // end of reduction
            q3_ptr = sym MODULUS_3_STATIC,
            a_ptr = in(reg) a.as_ptr(),
            b_ptr = in(reg) b.as_ptr(),
            out("rdx") _, 
            out("rdi") _, 
            out("r8") _, 
            out("r9") _, 
            out("r10") _, 
            out("r11") _, 
            out("r12") r0, 
            out("r13") r1, 
            out("r14") r2, 
            out("r15") r3,
            options(pure, readonly, nostack)
        );
    }

    [r0, r1, r2, r3]
}

#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(target_arch = "x86_64")]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
pub(crate) fn mont_square_with_reduction_impl(a: &[u64; 4]) -> [u64; 4] {
    // this is CIOS multiplication when top bit for top word of modulus is not set

    let mut r0: u64;
    let mut r1: u64;
    let mut r2: u64;
    let mut r3: u64;

    // let mut tmp: u64;

    unsafe {
        asm!(
            // round 0
            "mov rdx, qword ptr [{a_ptr} + 0]",
            "xor r8d, r8d",
            "mulx r10, r9, qword ptr [{a_ptr} + 8]",
            "mulx r15, r8, qword ptr [{a_ptr} + 16]",
            "mulx r12, r11, qword ptr [{a_ptr} + 24]",
            "adox r10, r8", 
            "adcx r11, r15", 
            "mov rdx, qword ptr [{a_ptr} + 8]",
            "mulx r15, r8, qword ptr [{a_ptr} + 16]",
            "mulx rcx, rdi, qword ptr [{a_ptr} + 24]",
            "mov rdx, qword ptr [{a_ptr} + 16]",
            "mulx r14, r13, qword ptr [{a_ptr} + 24]",
            "adox r11, r8", 
            "adcx r12, rdi", 
            "adox r12, r15", 
            "adcx r13, rcx",
            "mov r8, 0",
            "adox r13, r8",
            "adcx r14, r8", 

            // double
            "adox r9, r9",
            "adcx r12, r12", 
            "adox r10, r10",
            "adcx r13, r13", 
            "adox r11, r11",
            "adcx r14, r14", 

            // square contributions
            "mov rdx, qword ptr [{a_ptr} + 0]",
            "mulx rcx, r8, rdx",
            "mov rdx, qword ptr [{a_ptr} + 16]",
            "mulx rdi, rdx, rdx",
            "adox r12, rdx",
            "adcx r9, rcx",
            "adox r13, rdi",
            "mov rdx, qword ptr [{a_ptr} + 24]",
            "mulx r15, rcx, rdx",
            "mov rdx, qword ptr [{a_ptr} + 8]",
            "mulx rdx, rdi, rdx",
            "adcx r10, rdi",
            "adox r14, rcx",
            "mov rdi, 0",
            "adcx r11, rdx",
            "adox r15, rdi", // finish the addition chain
            "adcx r12, rdi",

            // reduction round 0
            "mov rdx, r8",
            "neg rdx",
            "xor rcx, rcx", // rcx = 0
            "adox r8, rdx",
            "mulx rdi, r8, qword ptr [rip + {q3_ptr}]", 
            "adcx r12, rdi",
            "adox r9, rcx",
            "adcx r13, rcx",
            "adox r10, rcx",
            "adcx r9, rcx",
            "adox r11, r8",
            "adcx r10, rcx",
            "adcx r11, rcx",
            "adox r12, rcx",
            "adcx r12, rcx",

            // reduction round 1
            "mov rdx, r9",
            "neg rdx",
            "xor rdi, rdi", // rdi = 0
            "adox r12, rdi",
            "mulx rcx, r8, qword ptr [rip + {q3_ptr}]", 
            "adcx r12, r8",
            "adox r13, rcx",
            "adcx r13, rdi",
            "adox r14, rdi",
            "adcx r9, rdx",
            "adox r10, rdi",
            "adcx r10, rdi",
            "adox r11, rdi",
            "adcx r11, rdi",

            // reduction round 2
            "mov rdx, r10",
            "neg rdx",
            "mulx r9, r8, qword ptr [rip + {q3_ptr}]", 
            "xor rcx, rcx", // rcx = 0
            "adox r12, rcx",
            "adcx r12, rcx",
            "adox r13, rcx",
            "adcx r13, r8",
            "adox r14, r9",
            "adcx r14, rcx",
            "adox r15, rcx",
            "adcx r10, rdx",
            "adox r11, rcx",
            "adcx r11, rcx",
            "adox r12, rcx",

            // reduction round 3
            "mov rdx, r11",
            "neg rdx",
            "xor rcx, rcx", // rcx = 0
            "adox r11, rdx",
            "mulx r11, r10, qword ptr [rip + {q3_ptr}]", 
            "adcx r12, rcx",
            "adox r12, rcx",
            "adcx r13, rcx",
            "adox r13, rcx",
            "adcx r14, r10",
            "adox r14, rcx",
            "adcx r15, r11",
            "adox r15, rcx",
            
            // reduction. We use sub/sbb
            "mov r8, r12",
            "sub r8, 1",
            "mov r9, r13",
            "sbb r9, 0",
            "mov r10, r14",
            "sbb r10, 0",
            "mov r11, r15",
            "mov rdx, qword ptr [rip + {q3_ptr}]",
            "sbb r11, rdx",

            // if CF == 1 then original result was ok (reduction wa not necessary)
            // so if not carry (CMOVNQ) then we copy 
            "cmovnc r12, r8",
            "cmovnc r13, r9",
            "cmovnc r14, r10",
            "cmovnc r15, r11",  
            // end of reduction
            // inv = const 0xffffffffffffffffu64,
            q3_ptr = sym MODULUS_3_STATIC,
            a_ptr = in(reg) a.as_ptr(),
            out("rcx") _, 
            out("rdx") _, 
            out("rdi") _, 
            out("r8") _, 
            out("r9") _, 
            out("r10") _, 
            out("r11") _, 
            out("r12") r0, 
            out("r13") r1, 
            out("r14") r2, 
            out("r15") r3,
            options(pure, readonly, nostack)
        );
    }

    // println!("Tmp = {}", tmp);

    [r0, r1, r2, r3]
}