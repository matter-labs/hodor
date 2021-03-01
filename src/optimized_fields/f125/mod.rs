use ff::*;

mod asm;
use asm::*;

#[derive(Copy, Clone, PartialEq, Eq, Default, Debug, Hash)]
pub struct FrRepr([u64; 2]);

const MODULUS_RAW: [u64; 2] = [1, MODULUS_LIMB_1];
const MODULUS: FrRepr = FrRepr(MODULUS_RAW);

// 3
const GENERATOR: FrRepr = FrRepr([0x720b1b19d49ea8f1, 0xbf4aa36101f13a58]);

// 2^S * t = MODULUS - 1 with t odd
const S: u32 = 64;

// 2^S root of unity computed by GENERATOR^t
const ROOT_OF_UNITY: FrRepr = FrRepr([0xaa9f02ab1d6124de, 0xb3524a6466112932]);

// -((2**256) mod s) mod s
const NEGATIVE_ONE: Fr = Fr(FrRepr([0xaa9f02ab1d6124de, 0xb3524a6466112932]));

const MODULUS_BITS: u32 = 125;
const REPR_SHAVE_BITS: u32 = 5;

static MONT_INV: u64 = 0xffffffffffffffff;
const R: FrRepr = FrRepr([0xffffffffffffffe1, 0xffffffffffffffff]);
const R2: FrRepr = FrRepr([0xfffffd737e000401, 0x00000001330fffff]);

const MODULUS_LIMB_1: u64 = 0x3000000300000001;

static ZERO_U64: u64 = 0;
static ONE_U64: u64 = 1;
static U64_MAX: u64 = 0xffffffffffffffff;

// static MODULUS_0_STATIC: u64 = 1;
static MODULUS_1_STATIC: u64 = 0x3000000300000001;

// 2^256 - MODULUS
static MODULUS_NEGATED_STATIC_0: u64 = 0xffffffffffffffff;
static MODULUS_NEGATED_STATIC_1: u64 = 0xcffffffcfffffffe;

impl FrRepr {

}

impl ::rand::Rand for FrRepr {
    #[inline(always)]
    fn rand<R: ::rand::Rng>(rng: &mut R) -> Self {
        FrRepr(rng.gen())
    }
}


impl ::std::fmt::Display for FrRepr
{
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "0x")?;
        for i in self.0.iter().rev() {
            write!(f, "{:016x}", *i)?;
        }

        Ok(())
    }
}

impl AsRef<[u64]> for FrRepr {
    #[inline(always)]
    fn as_ref(&self) -> &[u64] {
        &self.0
    }
}

impl AsMut<[u64]> for FrRepr {
    #[inline(always)]
    fn as_mut(&mut self) -> &mut [u64] {
        &mut self.0
    }
}

impl From<u64> for FrRepr {
    #[inline(always)]
    fn from(val: u64) -> FrRepr {
        let mut repr = Self::default();
        repr.0[0] = val;
        repr
    }
}


impl Ord for FrRepr {
    #[inline(always)]
    fn cmp(&self, other: &FrRepr) -> ::std::cmp::Ordering {
        for (a, b) in self.0.iter().rev().zip(other.0.iter().rev()) {
            if a < b {
                return ::std::cmp::Ordering::Less
            } else if a > b {
                return ::std::cmp::Ordering::Greater
            }
        }

        ::std::cmp::Ordering::Equal
    }
}

impl PartialOrd for FrRepr {
    #[inline(always)]
    fn partial_cmp(&self, other: &FrRepr) -> Option<::std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl PrimeFieldRepr for FrRepr {
    #[inline(always)]
    fn is_odd(&self) -> bool {
        self.0[0] & 1 == 1
    }

    #[inline(always)]
    fn is_even(&self) -> bool {
        !self.is_odd()
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.0[0] == 0 && self.0[1] == 0
    }

    #[inline(always)]
    fn shr(&mut self, mut n: u32) {
        if n >= 64 * 4 {
            *self = Self::from(0);
            return;
        }

        while n >= 64 {
            let mut t = 0;
            for i in self.0.iter_mut().rev() {
                ::std::mem::swap(&mut t, i);
            }
            n -= 64;
        }

        if n > 0 {
            let mut t = 0;
            for i in self.0.iter_mut().rev() {
                let t2 = *i << (64 - n);
                *i >>= n;
                *i |= t;
                t = t2;
            }
        }
    }

    #[inline(always)]
    fn div2(&mut self) {
        let mut t = 0;
        for i in self.0.iter_mut().rev() {
            let t2 = *i << 63;
            *i >>= 1;
            *i |= t;
            t = t2;
        }
    }

    #[inline(always)]
    fn mul2(&mut self) {
        *self = FrRepr(self::asm::double_nocarry_impl(self.0));
    }

    #[inline(always)]
    fn shl(&mut self, mut n: u32) {
        if n >= 64 * 4 {
            *self = Self::from(0);
            return;
        }

        while n >= 64 {
            let mut t = 0;
            for i in &mut self.0 {
                ::std::mem::swap(&mut t, i);
            }
            n -= 64;
        }

        if n > 0 {
            let mut t = 0;
            for i in &mut self.0 {
                let t2 = *i >> (64 - n);
                *i <<= n;
                *i |= t;
                t = t2;
            }
        }
    }

    #[inline(always)]
    fn num_bits(&self) -> u32 {
        let mut ret = (4 as u32) * 64;
        for i in self.0.iter().rev() {
            let leading = i.leading_zeros();
            ret -= leading;
            if leading != 64 {
                break;
            }
        }

        ret
    }

    #[inline(always)]
    fn add_nocarry(&mut self, other: &FrRepr) {
        *self = FrRepr(self::asm::add_nocarry_impl(self.0, other.0));
    }

    #[inline(always)]
    fn sub_noborrow(&mut self, other: &FrRepr) {
        *self = FrRepr(self::asm::sub_noborrow_impl(self.0, other.0));
    }
}

#[derive(Copy, Clone, PartialEq, Eq, Default, Debug, Hash)]
pub struct Fr(FrRepr);

impl ::std::fmt::Display for Fr
{
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Fr({})", self.into_repr())
    }
}

impl ::rand::Rand for Fr {
    fn rand<R: ::rand::Rng>(rng: &mut R) -> Self {
        loop {
            let mut tmp = Fr(FrRepr::rand(rng));

            // Mask away the unused bits at the beginning.
            tmp.0.as_mut()[1] &= 0xffffffffffffffff >> REPR_SHAVE_BITS;

            if tmp.is_valid() {
                return tmp
            }
        }
    }
}

impl From<Fr> for FrRepr {
    fn from(e: Fr) -> FrRepr {
        e.into_repr()
    }
}

impl PrimeField for Fr {
    type Repr = FrRepr;

    fn from_repr(r: FrRepr) -> Result<Fr, PrimeFieldDecodingError> {
        let mut r = Fr(r);
        if r.is_valid() {
            r.mul_assign(&Fr(R2));

            Ok(r)
        } else {
            Err(PrimeFieldDecodingError::NotInField(format!("{}", r.0)))
        }
    }

    fn from_raw_repr(r: FrRepr) -> Result<Fr, PrimeFieldDecodingError> {
        let r = Fr(r);
        if r.is_valid() {
            Ok(r)
        } else {
            Err(PrimeFieldDecodingError::NotInField(format!("{}", r.0)))
        }
    }

    fn into_repr(&self) -> FrRepr {
        const OTHER: Fr = Fr(FrRepr([1,0]));
        let mut tmp = *self;
        tmp.mul_assign(&OTHER);

        tmp.into_raw_repr()
    }

    fn into_raw_repr(&self) -> FrRepr {
        let r = *self;
        r.0
    }

    fn char() -> FrRepr {
        MODULUS
    }

    const NUM_BITS: u32 = MODULUS_BITS;

    const CAPACITY: u32 = Self::NUM_BITS - 1;

    fn multiplicative_generator() -> Self {
        Fr(GENERATOR)
    }

    const S: u32 = S;

    fn root_of_unity() -> Self {
        Fr(ROOT_OF_UNITY)
    }
}

impl Field for Fr {
    #[inline]
    fn zero() -> Self {
        Fr(FrRepr::from(0))
    }

    #[inline]
    fn one() -> Self {
        Fr(R)
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    #[inline]
    fn add_assign(&mut self, other: &Fr) {
        *self = Fr(FrRepr(add_with_reduction_impl(self.0.0, other.0.0)));
    }

    #[inline]
    fn double(&mut self) {
        *self = Fr(FrRepr(double_with_reduction_impl(self.0.0)));
    }

    #[inline]
    fn sub_assign(&mut self, other: &Fr) {
        *self = Fr(FrRepr(sub_with_reduction_impl(self.0.0, other.0.0)));
    }

    #[inline]
    fn negate(&mut self) {
        if !self.is_zero() {
            // unimplemented!()
            *self = Fr(FrRepr(sub_noborrow_impl(MODULUS.0, self.0.0)));
        }
    }

    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            // Guajardo Kumar Paar Pelzl
            // Efficient Software-Implementation of Finite Fields with Applications to Cryptography
            // Algorithm 16 (BEA for Inversion in Fp)

            let one = FrRepr::from(1);

            let mut u = self.0;
            let mut v = MODULUS;
            let mut b = Fr(R2); // Avoids unnecessary reduction step.
            let mut c = Self::zero();

            while u != one && v != one {
                while u.is_even() {
                    u.div2();

                    if b.0.is_even() {
                        b.0.div2();
                    } else {
                        b.0.add_nocarry(&MODULUS);
                        b.0.div2();
                    }
                }

                while v.is_even() {
                    v.div2();

                    if c.0.is_even() {
                        c.0.div2();
                    } else {
                        c.0.add_nocarry(&MODULUS);
                        c.0.div2();
                    }
                }

                if v < u {
                    u.sub_noborrow(&v);
                    b.sub_assign(&c);
                } else {
                    v.sub_noborrow(&u);
                    c.sub_assign(&b);
                }
            }

            if u == one {
                Some(b)
            } else {
                Some(c)
            }
        }
    }

    #[inline(always)]
    fn frobenius_map(&mut self, _: usize) {
        // This has no effect in a prime field.
    }

    #[inline]
    fn mul_assign(&mut self, other: &Fr)
    {
        *self = Fr(FrRepr(mont_mul_with_reduction_impl(self.0.0, other.0.0)));
    }

    #[inline]
    fn square(&mut self)
    {
        // unimplemented!()
        // *self = Fr(FrRepr(mont_square_with_reduction_impl(self.0.0)));
        self.mul_assign(&self.clone());
    }
}

impl Fr {
    #[inline(always)]
    fn is_valid(&self) -> bool {
        self.0 < MODULUS
    }

    #[inline(always)]
    fn reduce(&mut self) {
        unimplemented!()
        // *self = Fr(FrRepr(reduce_by_modulus_impl(self.0.0)));
    }
}

pub fn reg_pass(a: [u64; 2], b: [u64; 2], c: [u64; 2]) -> [u64; 2] {
    let tmp = add_nocarry_impl(a, b);
    add_nocarry_impl(tmp, c)
}

#[inline(always)]
pub fn cios_debug(a: &[u64; 2], b: &[u64; 2]) -> [u64; 2] {
    let mut m;

    const INV: u64 = 0xffffffffffffffffu64;

    // round 0 
    let b0 = b[0];
    let (r0, carry) = full_width_mul(a[0], b0);
    m = r0.wrapping_mul(INV);
    // everywhere semantic is arg0 + (arg1 * arg2)
    let red_carry = mac_by_value_return_carry_only(r0, m, 1);

    // loop over the rest
    let (r1, carry) = mac_by_value(carry, a[1], b0);
    let (r0, red_carry) = mac_with_carry_by_value(r1, m, MODULUS_LIMB_1, red_carry);

    // this will check overflow in debug
    let r1 = red_carry + carry;
    drop(b0);

    // round 1
    let b1 = b[1];
    let (r0, carry) = mac_by_value(r0, a[0], b1);
    m = r0.wrapping_mul(INV);
    // everywhere semantic is arg0 + (arg1 * arg2)
    let red_carry = mac_by_value_return_carry_only(r0, m, 1);

    // loop over the rest
    let (r1, carry) = mac_with_carry_by_value(r1, a[1], b1, carry);
    let (r0, red_carry) = mac_with_carry_by_value(r1, m, MODULUS_LIMB_1, red_carry);

    // this will check overflow in debug
    let r1 = red_carry + carry;
    drop(b1);

    return [r0, r1];
}

#[cfg(test)]
mod test {

    // #[test]
    // fn test_asm_steps() {
    //     let modulus: u128 = (0x3000000300000001 as u128) << 64 + 1;
    //     use super::*;
    //     use std::time::Instant;
    //     use rand::{XorShiftRng, SeedableRng, Rng};
    //     let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    //     let a = [0xfffffffffffffffc, 0x3000000300000000];
    //     let b = [0xfffffffffffffff0, 0x3000000300000000];

    //     let c_asm = super::asm::mont_mul_with_reduction_impl(a, b);
    //     let c_native = super::cios_debug(b, a);

    //     dbg!(c_asm);
    //     dbg!(c_native);
    // }


    // #[test]
    // fn test_asm_f125() {
    //     let modulus: u128 = (0x3000000300000001 as u128) << 64 + 1u128;
    //     use super::*;
    //     use std::time::Instant;
    //     use rand::{XorShiftRng, SeedableRng, Rng};
    //     let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    //     for i in 0..100 {
    //         let a: u128 = ((rng.gen::<u64>() as u128) << 64) + rng.gen::<u64>() as u128;
    //         let b: u128 = ((rng.gen::<u64>() as u128) << 64) + rng.gen::<u64>() as u128;

    //         let a = a % modulus;
    //         let b = b % modulus;

    //         let a = [a as u64, (a >> 64) as u64];
    //         let b = [b as u64, (b >> 64) as u64];

    //         let c_asm = super::asm::mont_mul_with_reduction_impl(a, b);
    //         let c_native = super::cios_debug(b, a);

    //         assert_eq!(c_asm, c_native, "failed at attempt {}", i);
    //     }
    // }

    // #[test]
    // fn try_pass_through_registers() {
    //     use super::*;
    //     use std::time::Instant;
    //     use rand::{XorShiftRng, SeedableRng, Rng};
    //     let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    //     const SIZE: usize = 1 << 24;
    //     let mut a: Fr = rng.gen();
    //     let b: Fr = rng.gen();
    //     let start = Instant::now();
    //     for _ in 0..SIZE {
    //         a.mul_assign(&b);
    //     }
    //     let ns_spent = start.elapsed().as_nanos() as f64;
    //     println!("Spent {} ns per multiplication", ns_spent / (SIZE as f64));
    // }

    // #[test]
    // fn test_asm_repr_correctness() {
    //     use super::*;
    //     use rand::{XorShiftRng, SeedableRng, Rng};
    //     use super::super::naive_f252::FrRepr as FrReprNaive;
    //     let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    //     const SIZE: usize = 1 << 10;
    //     for i in 0..SIZE {
    //         let a: Fr = rng.gen();
    //         let b: Fr = rng.gen();
    
    //         let mut a = a.into_raw_repr();
    //         let b = b.into_raw_repr();

    //         let mut a_n: FrReprNaive = unsafe {std::mem::transmute_copy(&a)};
    //         let b_n: FrReprNaive = unsafe {std::mem::transmute_copy(&b)};

    //         a.add_nocarry(&b);
    //         a_n.add_nocarry(&b_n);

    //         assert_eq!(unsafe { std::mem::transmute_copy::<_, FrReprNaive>(&a) }, a_n, "failed on attempt {}", i);

    //         a.sub_noborrow(&b);
    //         a_n.sub_noborrow(&b_n);

    //         assert_eq!(unsafe { std::mem::transmute_copy::<_, FrReprNaive>(&a) }, a_n, "failed on attempt {}", i);
    //     }
    // }

    #[test]
    fn test_asm_correctness() {
        use super::*;
        use rand::{XorShiftRng, SeedableRng, Rng};
        use super::super::naive_f125::Fr as FrNaive;
        let rng = &mut XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        const SIZE: usize = 1 << 20;
        let mut a: Fr = rng.gen();
        let b: Fr = rng.gen();

        let mut a_n: FrNaive = unsafe {std::mem::transmute(a)};
        let b_n: FrNaive = unsafe {std::mem::transmute(b)};
        for i in 0..SIZE {
            a.sub_assign(&b);
            a_n.sub_assign(&b_n);

            a.double();
            a_n.double();

            a.add_assign(&b);
            a_n.add_assign(&b_n);

            // a.square();
            // a_n.square();

            a.mul_assign(&b);
            a_n.mul_assign(&b_n);

            assert_eq!(unsafe { std::mem::transmute::<_, FrNaive>(a).into_raw_repr() }, a_n.into_raw_repr(), "failed on attempt {}", i);
            assert_eq!(unsafe { std::mem::transmute::<_, <FrNaive as PrimeField>::Repr>(a.into_repr()) }, a_n.into_repr(), "failed on attempt {}", i);
        }
    }
}