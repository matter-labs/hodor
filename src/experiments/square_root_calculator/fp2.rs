// use super::fq::{FROBENIUS_COEFF_FQ2_C1, Fq, NEGATIVE_ONE};
use ff::{Field, SqrtField};
use rand::{Rand, Rng};

use super::super::Fr as Fq;
use super::super::FrRepr as FqRepr;

use std::cmp::Ordering;

const NON_RESIDUE: Fq = Fq(FqRepr([
    0xffffffffffffffa1,
    0xffffffffffffffff,
    0xffffffffffffffff,
    0x07fffffffffff9b0
]));

const MINUS_ONE: Fq = Fq(FqRepr([
    0x0000000000000020,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000220
]));

// 0x20 0000 0000 0000 8800 0000 0000 0090 8000 0000 0000 0000 0000 0000 0000 0000 0800 0000 0000 0011 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000

const Q_SQUARED_MINUT_ONE_BY_TWO: [u64; 8] = [
    0x0,
    0x0,
    0x0,
    0x0800000000000011,
    0x0,
    0x8000000000000000,
    0x8800000000000090,
    0x20000000000000
];

const Q_MINUS_ONE_BY_TWO: [u64; 4] = [
    0x0,
    0x0,
    0x8000000000000000,
    0x400000000000008 
];

const Q_MINUS_ONE_BY_FOUR: [u64; 4] = [
    0x0,
    0x0,
    0x4000000000000000,
    0x200000000000004 
];

const E_PRECOMPUTED: Fq2 = Fq2 {
    c0: Fq(FqRepr([
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000
    ])),
    
    c1: Fq(FqRepr([
        0xb11079dbfdb6981f,
        0x701bb8e5e5b53751,
        0x6fe88e46d707bcab,
        0x02365bb6d67e6298
    ]))
};

const F_PRECOMPUTED: Fq2 = Fq2 {
    c0: Fq(FqRepr([
        0xffffffffffffff41,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0x07fffffffffff350
    ])),
    
    c1: Fq(FqRepr([
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000
    ]))
};

const MINUS_ONE_FP2: Fq2 = Fq2 {
    c0: MINUS_ONE,
    
    c1: Fq(FqRepr([
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000
    ]))
};

/// An element of Fq2, represented by c0 + c1 * u.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct Fq2 {
    pub c0: Fq,
    pub c1: Fq,
}

impl ::std::fmt::Display for Fq2 {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Fq2({} + {} * u)", self.c0, self.c1)
    }
}

/// `Fq2` elements are ordered lexicographically.
impl Ord for Fq2 {
    #[inline(always)]
    fn cmp(&self, other: &Fq2) -> Ordering {
        match self.c1.cmp(&other.c1) {
            Ordering::Greater => Ordering::Greater,
            Ordering::Less => Ordering::Less,
            Ordering::Equal => self.c0.cmp(&other.c0),
        }
    }
}

impl PartialOrd for Fq2 {
    #[inline(always)]
    fn partial_cmp(&self, other: &Fq2) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Fq2 {
    /// Norm of Fq2 as extension field in i over Fq
    pub fn norm(&self) -> Fq {
        let mut t0 = self.c0;
        let mut t1 = self.c1;
        t0.square();
        t1.square();
        t1.mul_assign(&NON_RESIDUE);
        t1.add_assign(&t0);

        t1
    }
}

impl Rand for Fq2 {
    fn rand<R: Rng>(rng: &mut R) -> Self {
        Fq2 {
            c0: rng.gen(),
            c1: rng.gen(),
        }
    }
}

impl Field for Fq2 {
    fn zero() -> Self {
        Fq2 {
            c0: Fq::zero(),
            c1: Fq::zero(),
        }
    }

    fn one() -> Self {
        Fq2 {
            c0: Fq::one(),
            c1: Fq::zero(),
        }
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn square(&mut self) {
        // v0 = c0 - c1
        let mut v0 = self.c0;
        v0.sub_assign(&self.c1);
        // v3 = c0 - beta * c1
        let mut v3 = self.c0;
        let mut t0 = self.c1;
        t0.mul_assign(&NON_RESIDUE);
        v3.sub_assign(&t0);
        // v2 = c0 * c1
        let mut v2 = self.c0;
        v2.mul_assign(&self.c1);

        // v0 = (v0 * v3) + v2
        v0.mul_assign(&v3);
        v0.add_assign(&v2);

        self.c1 = v2.clone();
        self.c1.double();
        self.c0 = v0;
        v2.mul_assign(&NON_RESIDUE);
        self.c0.add_assign(&v2);
    }

    fn double(&mut self) {
        self.c0.double();
        self.c1.double();
    }

    fn negate(&mut self) {
        self.c0.negate();
        self.c1.negate();
    }

    fn add_assign(&mut self, other: &Self) {
        self.c0.add_assign(&other.c0);
        self.c1.add_assign(&other.c1);
    }

    fn sub_assign(&mut self, other: &Self) {
        self.c0.sub_assign(&other.c0);
        self.c1.sub_assign(&other.c1);
    }

    fn mul_assign(&mut self, other: &Self) {
        let mut aa = self.c0;
        aa.mul_assign(&other.c0);
        let mut bb = self.c1;
        bb.mul_assign(&other.c1);

        let mut o = other.c0;
        o.add_assign(&other.c1);

        self.c1.add_assign(&self.c0);
        self.c1.mul_assign(&o);
        self.c1.sub_assign(&aa);
        self.c1.sub_assign(&bb);
        self.c0 = aa;
        bb.mul_assign(&NON_RESIDUE);
        self.c0.add_assign(&bb);
    }

    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        let mut v0 = self.c0;
        v0.square();
        let mut v1 = self.c1;
        v1.square();

        // v0 = v0 - beta * v1
        let mut v1_by_nonresidue = v1;
        v1_by_nonresidue.mul_assign(&NON_RESIDUE);
        v0.sub_assign(&v1_by_nonresidue);
        v0.inverse().map(|v1| {
            let mut c0 = self.c0;
            c0.mul_assign(&v1);
            let mut c1 = self.c1;
            c1.mul_assign(&v1);
            c1.negate();

            Self {
                c0: c0, 
                c1: c1,
            }
        })
    }

    fn frobenius_map(&mut self, power: usize) {
        unimplemented!();
        // self.c1.mul_assign(&FROBENIUS_COEFF_FQ2_C1[power % 2]);
    }
}

impl Fq2 {
    fn xi(&self) -> Self {
        self.pow(&Q_SQUARED_MINUT_ONE_BY_TWO)
    }

    pub fn mul_by_0(&mut self, fq: &Fq) {
        self.c0.mul_assign(&fq);
        self.c1.mul_assign(&fq);
    }

    pub fn mul_by_1(&mut self, fq: &Fq) {
        let aa = Fq::zero();
        let mut bb = self.c1;
        bb.mul_assign(&fq);

        self.c1.add_assign(&self.c0);
        self.c1.mul_assign(&fq);
        self.c1.sub_assign(&bb);
        self.c0 = aa;
        bb.mul_assign(&NON_RESIDUE);
        self.c0 = bb;
    }
}

impl SqrtField for Fq2 {
    fn legendre(&self) -> ::ff::LegendreSymbol {
        self.norm().legendre()
    }

    fn sqrt(&self) -> Option<Self> {
        // Algorithm 9, https://eprint.iacr.org/2012/685.pdf

        if self.is_zero() {
            Some(Self::zero())
        } else {
            use ff::PrimeField;

            let b = self.pow(&Q_MINUS_ONE_BY_FOUR);

            let mut b_squared = b;
            b_squared.square();

            // TODO: Implement Frobenius map
            let b_in_q = b.pow(Fq::char());

            let mut b_in_q_by_b = b_in_q;
            b_in_q_by_b.mul_assign(&b);

            let mut a0 = b_in_q_by_b;
            a0.square();

            if a0 == MINUS_ONE_FP2 {
                return None;
            }

            if b_in_q_by_b == Fq2::one() {
                let mut b_squared_by_a = b_squared;
                b_squared_by_a.mul_assign(&self);
                debug_assert!(b_squared_by_a.c1.is_zero());
                let x0 = b_squared_by_a.c0.sqrt();
                if x0.is_none() {
                    return None;
                }

                let mut x = b_in_q;
                x.mul_by_0(&x0.expect("is some"));

                return Some(x);
            } else {
                let mut b_squared_by_a_by_f = b_squared;
                b_squared_by_a_by_f.mul_assign(&self);
                debug_assert!(F_PRECOMPUTED.c1.is_zero());
                b_squared_by_a_by_f.mul_by_0(&F_PRECOMPUTED.c0);
                debug_assert!(b_squared_by_a_by_f.c1.is_zero());

                let x0 = b_squared_by_a_by_f.c0.sqrt();
                if x0.is_none() {
                    return None;
                }

                let mut x = b_in_q;
                x.mul_by_0(&x0.expect("is some"));
                debug_assert!(E_PRECOMPUTED.c0.is_zero());
                x.mul_by_1(&E_PRECOMPUTED.c1);

                return Some(x);
            }
        }
    }
}

#[test]
fn find_c() {
    // Precomputation part of Alg 10 from https://eprint.iacr.org/2012/685.pdf
    let one_fp2 = Fq2::one();

    let mut c = Fq2 {
        c0: Fq::one(),
        c1: Fq::one(),
    };
    loop {
        let xi = c.xi();
        if xi != one_fp2 {
            break;
        } else {
            c.add_assign(&one_fp2);
        }
    }
    println!("C = {}", c);

    let d = c.pow(&Q_MINUS_ONE_BY_TWO);
    println!("D = {}", d);

    let mut dc = d;
    dc.mul_assign(&c);

    let e = dc.inverse().unwrap();

    let mut f = dc;
    f.square();


    let mut may_be_one = e;
    may_be_one.mul_assign(&dc);

    assert!(may_be_one == one_fp2);

    println!("E = {}", e);
    println!("F = {}", f);

    use ff::PrimeField;

    let repr = e.c0.into_raw_repr();
    for el in repr.as_ref().iter() {
        println!("0x{:016x}", el);
    }
    println!("");
    let repr = e.c1.into_raw_repr();
    for el in repr.as_ref().iter() {
        println!("0x{:016x}", el);
    }
    println!("");
    let repr = f.c0.into_raw_repr();
    for el in repr.as_ref().iter() {
        println!("0x{:016x}", el);
    }
    println!("");
    let repr = f.c1.into_raw_repr();
    for el in repr.as_ref().iter() {
        println!("0x{:016x}", el);
    }
    println!("");
    let mut minus_one = Fq::one();
    minus_one.negate();
    let repr = minus_one.into_raw_repr();
    for el in repr.as_ref().iter() {
        println!("0x{:016x}", el);
    }
}

#[test]
fn test_fq2_squaring() {
    use rand::{XorShiftRng, SeedableRng};
    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    for _ in 0..1000 {
        let mut a = Fq2::rand(&mut rng);
        let mut b = a;
        b.mul_assign(&a);
        a.square();

        assert_eq!(a, b);
    }
}

#[test]
fn test_sparse_mul() {
    use crate::ff::PrimeField;

    let a = Fq2 {
        c0: Fq::one(),
        c1: Fq::one(),
    }; // u + 1

    let seven = Fq::from_str("7").unwrap();

    let b = Fq2 {
        c0: Fq::zero(),
        c1: seven,
    };

    let c = Fq2 {
        c0: seven,
        c1: Fq::zero(),
    };

    let mut naive_mul_0 = a;
    naive_mul_0.mul_assign(&c);

    let mut fast_mul_0 = a;
    fast_mul_0.mul_by_0(&seven);

    assert_eq!(naive_mul_0, fast_mul_0);

    let mut naive_mul_1 = a;
    naive_mul_1.mul_assign(&b);

    let mut fast_mul_1 = a;
    fast_mul_1.mul_by_1(&seven);

    assert_eq!(naive_mul_1, fast_mul_1);
}

#[test]
fn test_square_root() {
    use rand::{XorShiftRng, SeedableRng};
    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    let mut cnt = 0;
    for i in 0..1000 {
        let a = Fq2::rand(&mut rng);
        if let Some(sqrt_a) = a.sqrt() {
            let mut may_be_a = sqrt_a;
            may_be_a.square();
            assert_eq!(may_be_a, a, "Failed at attempt {}, counter {}", i, cnt);
            cnt += 1;
        }
    }

    assert!(cnt > 0);
}

#[test]
fn test_bench_square_root() {
    use ff::PrimeField;
    use rand::{XorShiftRng, SeedableRng};
    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    // const NUM_ROUNDS: usize = (1 << 20) - 1;
    const NUM_ROUNDS: usize = (1 << 18) - 1;

    let k0 = Fq::from_str("3").unwrap();
    let k1 = Fq::from_str("5").unwrap();

    let mut qr = Fq2::zero();
    for _ in 0..1000 {
        let a = Fq2::rand(&mut rng);
        if let Some(_) = a.sqrt() {
            qr = a;
            break;
        }
    }

    let now = std::time::Instant::now();

    for _ in 0..NUM_ROUNDS {
        let mut sqrt = qr.sqrt();
        if sqrt.is_none() {
            qr.negate();
            sqrt = qr.sqrt();
        }
        assert!(sqrt.is_some());
        let sqrt = sqrt.unwrap();
        qr = sqrt;

        // qr.c0 = sqrt.c1;
        // qr.c0.add_assign(&k0);

        // qr.c1 = sqrt.c0;
        // qr.c1.add_assign(&k1);
    }

    println!("{} square root calculations taken {}ms", NUM_ROUNDS, now.elapsed().as_millis());
}

#[test]
fn test_non_redisue() {
    let non_res = NON_RESIDUE;
    println!("Non-residue = {}", non_res);
    let non_res_sqrt = non_res.sqrt();
    assert!(non_res_sqrt.is_none());
}

#[test]
fn test_minus_one() {
    let non_res = MINUS_ONE;
    println!("Minus one = {}", non_res);
    let non_res_sqrt = non_res.sqrt();
    assert!(non_res_sqrt.is_some());
}

