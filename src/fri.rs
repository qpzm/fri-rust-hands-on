use crate::merkle::MerkleTree;
use crate::utils::{
    evaluate, fold_polynomial, get_evaluation_points, get_omega, reconstruct_merkle_root,
};
use ark_ff::PrimeField;

pub struct FRI<E: PrimeField> {
    blowup_factor: u64,
    prover: Prover<E>,
    verifier: Verifier<E>,
}

// Items that need to be stored by the prover
pub struct Prover<E: PrimeField> {
    // Initial polynomial and all the intermediate folded polynomials (in coefficient form)
    polys: Vec<Vec<E>>,
    // Merkle trees constructed from initial polynomial and intermediate folder polynomials (in evaluation form)
    merkle_trees: Vec<MerkleTree<E>>,
}

impl<E: PrimeField> Default for Prover<E> {
    fn default() -> Self {
        Self {
            polys: Vec::new(),
            merkle_trees: Vec::new(),
        }
    }
}

pub struct Verifier<E: PrimeField> {
    // Random values that are sent to the prover during the commit phase (in reality this would be fiat shamir)
    random_values: Vec<E>,
    // Commitments that are provided by the prover
    commitments: Vec<E>,
}

impl<E: PrimeField> Default for Verifier<E> {
    fn default() -> Self {
        Self {
            random_values: Vec::new(),
            commitments: Vec::new(),
        }
    }
}

#[derive(Debug)]
pub struct RoundProof<E: PrimeField> {
    omega: E,
    f_g: E,
    f_negative_g: E,
    merkle_path: Vec<E>,
}

impl<E: PrimeField> FRI<E> {
    pub fn new(blowup_factor: u64) -> Self {
        Self {
            blowup_factor,
            prover: Prover::default(),
            verifier: Verifier::default(),
        }
    }

    pub fn commit(&mut self, mut poly: Vec<E>, random_values: Vec<E>) -> Vec<E> {
        // for a polynomial with length 2^n, it requires folding n times to get to a plaintext
        assert_eq!(2u32.pow(random_values.len() as u32), poly.len() as u32);

        // vector that stores the commitment at each round
        let mut commitments: Vec<E> = vec![];

        // NOTE: doesn't include the final plaintext
        let mut polys: Vec<Vec<E>> = vec![];

        // vector that stores the merkle tree for the initial polynomial and folded polynomials
        let mut merkle_trees: Vec<MerkleTree<E>> = vec![];

        // roots of unity point
        let mut omega = get_omega(&poly);

        let mut random_value_idx = 0;
        while poly.len() > 1 {
            let eval_points = get_evaluation_points(&poly, omega, self.blowup_factor);
            let merkle_tree = MerkleTree::new(eval_points);
            let merkle_root = merkle_tree.root;
            commitments.push(merkle_root);
            polys.push(poly.clone());
            merkle_trees.push(merkle_tree);

            poly = fold_polynomial(&poly, random_values[random_value_idx]);
            random_value_idx += 1;
            // evaluation points for the next round are used using the square of current roots of unity point
            omega = omega.pow([2]);
        }

        // the final commitment is a plaintext
        assert_eq!(poly.len(), 1);
        commitments.push(poly[0]);

        self.prover = Prover {
            polys,
            merkle_trees,
        };

        self.verifier = Verifier {
            random_values,
            commitments: commitments.clone(),
        };

        commitments
    }

    pub fn query(&self, mut index: u64) -> Vec<RoundProof<E>> {
        // when verifier passes the prover with an evaluation point (a point within the roots of unity)
        // the prover sends the proof for each round to the verifier
        let mut proofs = vec![];
        let mut omega = get_omega(&self.prover.polys[0]).pow([index]);
        for (poly, merkle_tree) in self
            .prover
            .polys
            .iter()
            .zip(self.prover.merkle_trees.iter())
        {
            if index >= poly.len() as u64 {
                index = index % poly.len() as u64;
            }

            dbg!(index);
            dbg!(poly.len());

            let proof = RoundProof {
                omega,
                f_g: evaluate(poly, omega),
                f_negative_g: evaluate(poly, -omega),
                merkle_path: merkle_tree.get_merkle_path(index),
            };

            proofs.push(proof);

            omega = omega.pow([2]);
        }

        proofs
    }

    pub fn verify(&self, proofs: Vec<RoundProof<E>>) -> bool {
        let mut i = 0;
        let mut previous_proof: Option<&RoundProof<E>> = None;
        for (proof, commitment) in proofs.iter().zip(self.verifier.commitments.iter()) {
            let merkle_root = reconstruct_merkle_root(proof.f_g, &proof.merkle_path);
            if &merkle_root != commitment {
                return false;
            }

            // check that the evaluation matches the evaluations of the previous round
            match previous_proof {
                Some(prev_round) => {
                    // f1(x^2) = (x+r1)(f0(x)) / 2x + (r1-x)(f0(-x))/2(-x)
                    if proof.f_g != (prev_round.omega + self.verifier.random_values[i - 1])
                        * prev_round.f_g
                        / (E::from(2) * prev_round.omega)
                        + (self.verifier.random_values[i - 1] - prev_round.omega)
                        * prev_round.f_negative_g
                        / (E::from(2) * -prev_round.omega)
                    {
                        return false
                    }
                }
                None => {}
            }

            previous_proof = Some(proof);
            i += 1;
        }

        let final_round = previous_proof.expect("final round proof must exist");
        let final_commitment = self
            .verifier
            .commitments
            .last()
            .expect("final commitment must exist")
            .to_owned();

        if final_commitment
            != (final_round.omega + self.verifier.random_values[i - 1]) * (final_round.f_g)
                / (E::from(2) * final_round.omega)
                + (self.verifier.random_values[i - 1] - final_round.omega)
                    * (final_round.f_negative_g)
                    / (E::from(2) * -final_round.omega)
        {
            return false;
        }

        true
    }
}
