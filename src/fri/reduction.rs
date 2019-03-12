/// Reduction domain is a FRI domain that with arbitrary power-of-two choice
/// for a choice of multivariate polynomial Q(X, Y) = P_0(Y) + X*P_1(Y)
/// 
/// FRI with two-fold reduction per step
/// Examples: f^{0}(x) = P^{0}(x) = P_0^{1}(x2) + x * P_1^{1}(x2) is of degree n
/// transforms into multivariate poly Q^{1}(X, Y) = P_0^{1}(Y) + X * P_1^{1}(Y), that is only 
/// valid when Y == X^2
/// deg_X(Q^{1}) < 2
/// deg_Y(Q^{1}) = n/2
/// defines a map q_0(X) = X^{2}
/// domain of Y becomes n/2 order multiplicative group
/// domain of X if F - the field itself
/// proving: sample x_0 at uniformly at random from F, ask a prover to commit to the values of 
/// Q^{1}(X, Y) at x0 as Q^{1}(x_0, Y)
/// verification: sample beta, such as alpha_0, alpha_1: alpha_0^2 = alpha_1^2 = beta,
/// and alpha_0, alpha_1 are in the multiplicative subgroup of a size 2^n of the field F,
/// so beta is in a multiplicative group or 2^(n-1)
/// sample answers from the f^{0} at points alpha_0, alpha_1 as p_0 and p_1 and from f^{1} at beta
/// interpolate the polynomial between p_0 and p_1 as p(X), and if p(x_0) == beta
/// 
/// Now for Y = X^4
/// Q^{1}(X, Y) = P_0^{1}(Y) + X * P_1^{1}(Y) + X^2 * P_2^{1}(Y) + X^3 * P_3^{1}(Y)

use ff::{
    Field, 
    PrimeField
};

use crate::utils::poly::*;

pub struct ReductionDomain<F: PrimeField> {
    coeffs: Vec<F>,
    exp: u32,
    domain_size: usize,
    omega: F,
    omegainv: F,
    geninv: F,
    minv: F,
    nonzero_coeffs: usize,
    reduction_degree: usize,
    rho_exponent: u32
}

impl<F: PrimeField> ReductionDomain<F> {
    pub fn as_ref(&self) -> &[F] {
        &self.coeffs
    }

    pub fn as_mut(&mut self) -> &mut [F] {
        &mut self.coeffs
    }

    pub fn into_coeffs(self) -> Vec<F> {
        self.coeffs
    }

    pub fn from_coeffs(mut coeffs: Vec<F>, reduction_degree: usize, rho_exponent: u32) -> Result<ReductionDomain<F>, ()>
    {
        let coeffs_len = coeffs.len();

        // we need an underlying evaluation domain to be
        // 2^(k + rho) , where k is the smallest integer for
        // 2^k > degree

        let mut m = 1;
        let mut exp = 0;
        let mut omega = F::root_of_unity();
        while m < coeffs_len {
            m *= 2;
            exp += 1;
            if exp > F::S {
                return Err(())
            }
        }

        // nonzero coeffs is power of 2;
        let nonzero_coeffs = m;

        // now exp is a smallest power of two to keep our original polynomial
        // and we can extend it by rho
        exp = exp + rho_exponent;

        let mut domain_size = m;

        for _ in 0..rho_exponent {
            domain_size *= 2;
        }

        let max_degree = (1 << F::S) - 1;

        if domain_size > max_degree {
            return Err(())
        }

        // reduce domain if it's excessive
        for _ in exp..F::S {
            omega.square();
        }

        // Extend the coeffs vector with zeroes if necessary
        coeffs.resize(domain_size, F::zero());

        // we perfectly know that only small portion of coefficients is zero
        // and keep it explicitly

        Ok(ReductionDomain {
            coeffs: coeffs,
            exp: exp,
            domain_size: domain_size,
            omega: omega,
            omegainv: omega.inverse().unwrap(),
            geninv: F::multiplicative_generator().inverse().unwrap(),
            minv: F::from_str(&format!("{}", m)).unwrap().inverse().unwrap(),
            nonzero_coeffs: nonzero_coeffs,
            reduction_degree: reduction_degree,
            rho_exponent: rho_exponent
        })
    }

    pub fn from_coeffs_into_sized(mut coeffs: Vec<F>, size: usize, reduction_degree: usize, rho_exponent: u32) -> Result<ReductionDomain<F>, ()>
    {
        coeffs.resize(size, F::zero());

        Self::from_coeffs(coeffs, reduction_degree, rho_exponent)
    }

    // return new 
    pub fn reduce(&self, challenge: F) -> Result<(ReductionDomain<F>, Vec<F>), ()> {
        // make components for multivariate polynomials
        assert!(self.coeffs.len() % self.reduction_degree == 0);

        let subsizes = self.nonzero_coeffs / self.reduction_degree;
        let mut subpoly = vec![vec![F::zero(); subsizes]; self.reduction_degree];
        let reduction_degree = self.reduction_degree;

        for i in 0..self.nonzero_coeffs {
            let element_index = i / reduction_degree;
            let subpoly_index = i % reduction_degree;
            subpoly[subpoly_index][element_index] = self.coeffs[i];
        }

        let expected_nonzero_elements = self.nonzero_coeffs / self.reduction_degree;

        // get a generator for a reduced domain
        let mut next_domain_generator = self.omega;
        for _ in 0..self.reduction_degree {
            next_domain_generator.mul_assign(&self.omega);
        } 

        // now make evaluation

        // distribute powers of challenge into each of the subpolys
        
        let mut q: Option<Vec<F>> = None;
        {
            let mut tmp = F::one();
            for subpoly in subpoly.into_iter() {
                if let Some(q_poly) = q.as_mut() {
                    mul_add_polynomials(&mut q_poly[..], &subpoly[..], tmp);
                } else {
                    q = Some(subpoly);
                }
                tmp.mul_assign(&challenge);
            }
        }

        // now perform an evaluation over the reduced domain
        let domain_size = self.domain_size as usize;
        let mut oracle_values = vec![F::zero(); domain_size];
        let mut tmp = next_domain_generator;
        let q = q.unwrap();

        for i in 0..domain_size {
            oracle_values[i] = evaluate_at_consequitive_powers(& q[..expected_nonzero_elements], F::one(), tmp);
            tmp.mul_assign(&next_domain_generator);
        }

        let next_domain = ReductionDomain::from_coeffs(q, self.reduction_degree, self.rho_exponent).unwrap();
        
        Ok((next_domain, oracle_values))
    }

    pub fn create_oracle(&self) -> Result<Vec<F>, ()> {
        let mut result = vec![F::zero(); self.domain_size];
        let mut tmp = self.omega;
        // this is a naive evaluation
        for i in 0..self.domain_size {
            result[i] = evaluate_at_consequitive_powers(& self.coeffs[..self.nonzero_coeffs], F::one(), tmp);
            tmp.mul_assign(&self.omega);
        }

        Ok(result)
    }

    pub fn query(&self, point: F) -> Result<F, ()> {
        // ensure that the point is in the domain
        let tmp = point.pow([self.domain_size as u64]);
        if tmp != F::one() {
            return Err(());
        }

        Ok(evaluate_at_consequitive_powers(& self.coeffs[..self.nonzero_coeffs], F::one(), point))       
    }

    pub fn domain_size(&self) -> usize {
        self.domain_size   
    }

    pub fn exp(&self) -> u32 {
        self.exp   
    }

    pub fn domain_generator(&self) -> F {
        self.omega
    }
}


#[test]
fn test_simple_reduction() {
    use crate::Fr;

    let poly = vec![
        Fr::from_str("1").unwrap(),
        Fr::from_str("2").unwrap(),
        Fr::from_str("3").unwrap(),
        Fr::from_str("4").unwrap(),
    ];

    let domain = ReductionDomain::from_coeffs(poly, 2, 4).unwrap();

    let gen = domain.domain_generator();

    let pow = 3u64;

    let z = gen.pow([pow]);

    let oracle = domain.create_oracle().unwrap();

    let domain_size = domain.domain_size();

    assert_eq!(gen.pow([domain_size as u64]), Fr::one(), "domain generator should indeed be for a right size group");
    let alpha_0_index = 5;
    let alpha_1_index = alpha_0_index + domain_size / 2;

    let alpha_0 = gen.pow([alpha_0_index as u64]);
    let alpha_1 = gen.pow([alpha_1_index as u64]);

    let mut beta = alpha_0;
    beta.square();

    let mut beta_1 = alpha_1;
    beta_1.square();

    assert_eq!(beta, beta_1, "alphas should correspond to beta");

    let s_0 = domain.query(alpha_0).unwrap();
    let s_1 = domain.query(alpha_1).unwrap();

    // index-1 is important, cause index is indeed a power (starts with 1), while vector is zero indexed
    let s_0_oracle = oracle[alpha_0_index-1];
    let s_1_oracle = oracle[alpha_1_index-1];

    assert_eq!(s_0, s_0_oracle, "query and oracle should match for first point");
    assert_eq!(s_1, s_1_oracle, "query and oracle should match for second point");

    let (next_domain, _next_oracle) = domain.reduce(z).unwrap();

    let y = next_domain.query(beta).unwrap();

    let interpolate = interpolate(& vec![(alpha_0, s_0), (alpha_1, s_1)]).unwrap();
    assert_eq!(interpolate.len(), 2, "interpolant should be of degree 1");
    let value = evaluate_at_consequitive_powers(& interpolate[..], Fr::one(), z);
    assert_eq!(value, y, "oracle must be consistent");
}

#[test]
fn test_simple_reduction_of_degree_4() {
    use crate::Fr;

    let poly = vec![
        Fr::from_str("1").unwrap(),
        Fr::from_str("2").unwrap(),
        Fr::from_str("3").unwrap(),
        Fr::from_str("4").unwrap(),
    ];

    let domain = ReductionDomain::from_coeffs(poly, 4, 4).unwrap();

    let gen = domain.domain_generator();

    let pow = 3u64;

    let z = gen.pow([pow]);

    let oracle = domain.create_oracle().unwrap();

    let domain_size = domain.domain_size();

    assert_eq!(gen.pow([domain_size as u64]), Fr::one(), "domain generator should indeed be for a right size group");
    let alpha_0_index = 5;
    let alpha_1_index = alpha_0_index + domain_size / 2;

    let alpha_0 = gen.pow([alpha_0_index as u64]);
    let alpha_1 = gen.pow([alpha_1_index as u64]);

    let mut beta = alpha_0;
    beta.square();

    let mut beta_1 = alpha_1;
    beta_1.square();

    assert_eq!(beta, beta_1, "alphas should correspond to beta");

    let s_0 = domain.query(alpha_0).unwrap();
    let s_1 = domain.query(alpha_1).unwrap();

    // index-1 is important, cause index is indeed a power (starts with 1), while vector is zero indexed
    let s_0_oracle = oracle[alpha_0_index-1];
    let s_1_oracle = oracle[alpha_1_index-1];

    assert_eq!(s_0, s_0_oracle, "query and oracle should match for first point");
    assert_eq!(s_1, s_1_oracle, "query and oracle should match for second point");

    let (next_domain, _next_oracle) = domain.reduce(z).unwrap();

    let y = next_domain.query(beta).unwrap();

    let interpolate = interpolate(& vec![(alpha_0, s_0), (alpha_1, s_1)]).unwrap();
    assert_eq!(interpolate.len(), 2, "interpolant should be of degree 1");
    let value = evaluate_at_consequitive_powers(& interpolate[..], Fr::one(), z);
    assert_eq!(value, y, "oracle must be consistent");
}