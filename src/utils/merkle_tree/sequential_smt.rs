// Sparse Merkle tree with flexible hashing strategy

use std::collections::HashMap;
use super::TreeHasher;

// Tree of depth 0 should contain ONE element that is also a root
// Tree of depth 1 should contain TWO elements
// Tree of depth 20 should contain 2^20 elements

// [0, (2^TREE_DEPTH - 1)]
type ItemIndex = u32;

// [0, TREE_DEPTH]
type Depth = u32;

// Hash index determines on what level of the tree the hash is 
// and kept as level (where zero is a root) and item in a level indexed from 0
type HashIndex = (u32, u32);

type ItemIndexPacked = u64;

trait PackToIndex {
    fn pack(&self) -> ItemIndexPacked;
}

impl PackToIndex for HashIndex {
    fn pack(&self) -> ItemIndexPacked {
        let mut packed = 0u64;
        packed += u64::from(self.0);
        packed <<= 32u64;
        packed += u64::from(self.1);

        packed
    }
}

pub trait GetBytes {
    fn be_bytes(&self) -> Vec<u8>;
}

#[derive(Debug, Clone)]
pub struct SparseMerkleTree<T: GetBytes + Default, TreeHasher>
{
    tree_depth: Depth,
    pub prehashed: Vec<Vec<u8>>,
    pub items: HashMap<ItemIndex, T>,
    pub hashes: HashMap<ItemIndexPacked, Vec<u8>>,
    _marker: std::marker::PhantomData<TreeHasher>
}

impl<T: GetBytes + Default, H: TreeHasher> SparseMerkleTree<T, H>
{

    pub fn new(tree_depth: Depth) -> Self {
        let items = HashMap::new();
        let hashes = HashMap::new();
        // we need to make sparse hashes for tree depth levels
        let mut prehashed = Vec::with_capacity((tree_depth + 1) as usize);
        let mut cur = H::hash_leaf(&T::default().be_bytes());
        prehashed.push(cur.clone());

        for i in 0..tree_depth {
            cur = H::compress(&cur, &cur, i as usize);
            prehashed.push(cur.clone());
        }
        prehashed.reverse();

        assert_eq!(prehashed.len() - 1, tree_depth as usize);
        Self{
            tree_depth: tree_depth, 
            prehashed: prehashed, 
            items: items, 
            hashes: hashes, 
            _marker: std::marker::PhantomData
        }
    }

    // How many items can the tree hold
    pub fn capacity(&self) -> u32 {
        1 << self.tree_depth
    }

    pub fn insert(&mut self, index: ItemIndex, item: T) {
        assert!(index < self.capacity());
        let hash_index = (self.tree_depth, index);

        let item_bytes = item.be_bytes();

        let hash = H::hash_leaf(&item_bytes);

        self.hashes.insert(hash_index.pack(), hash);

        self.items.insert(index, item);

        let mut next_level = (hash_index.0, hash_index.1);

        for _ in 0..next_level.0 {
            next_level = (next_level.0 - 1, next_level.1 >> 1);
            self.update_hash(next_level);

        }

        assert_eq!(next_level.0, 0);
    }

    pub fn delete(&mut self, index: ItemIndex) {
        assert!(index < self.capacity());
        let hash_index = (self.tree_depth, index);

        let item = T::default();

        let item_bytes = item.be_bytes();

        let hash = H::hash_leaf(&item_bytes);

        self.hashes.insert(hash_index.pack(), hash);

        self.items.insert(index, item);

        let mut next_level = (hash_index.0, hash_index.1);

        for _ in 0..next_level.0 {
            next_level = (next_level.0 - 1, next_level.1 >> 1);
            self.update_hash(next_level);

        }
        
        assert_eq!(next_level.0, 0);
    }

    fn update_hash(&mut self, index: HashIndex) -> Vec<u8> {
        // should NEVER be used to update the leaf hash

        assert!(index.0 < self.tree_depth);
        assert!(index.1 < self.capacity());

        // indices for child nodes in the tree
        let lhs_index = (index.0 + 1, (index.1 << 1));
        let rhs_index = (index.0 + 1, (index.1 << 1) + 1);

        let lhs_hash = self.get_hash(lhs_index);
        let rhs_hash = self.get_hash(rhs_index);

        let hash = H::compress(&lhs_hash, &rhs_hash, (self.tree_depth - 1 - index.0) as usize);

        self.hashes.insert(index.pack(), hash.clone());
        hash

    }

    pub fn get_hash(&self, index: HashIndex) -> Vec<u8> {
        assert!(index.0 <= self.tree_depth);
        assert!(index.1 < self.capacity());

        if let Some(hash) = self.hashes.get(&index.pack()) {
            // if hash for this index exists, return it
            hash.clone()
        } else {
            // otherwise return pre-computed
            self.prehashed.get((index.0) as usize).unwrap().clone()
        }
    }

    pub fn merkle_path(&self, index: ItemIndex) -> Vec<(Vec<u8>, bool)> {
        assert!(index < self.capacity());
        let mut hash_index = (self.tree_depth, index);

        (0..self.tree_depth).rev().map(|_level| {
            let dir = (hash_index.1 & 1) > 0;
            let proof_index = (hash_index.0, hash_index.1 ^ 1);
            let hash = self.get_hash(proof_index);
            hash_index = (hash_index.0 - 1, hash_index.1 >> 1);
            (hash, dir)
        }).collect()
    }

    pub fn verify_proof(&self, index: ItemIndex, item: T, proof: Vec<(Vec<u8>, bool)>) -> bool {
        assert!(index < self.capacity());
        let item_bytes = item.be_bytes();
        let mut hash = H::hash_leaf(&item_bytes);
        let mut proof_index: ItemIndex = 0;

        for (i, e) in proof.clone().into_iter().enumerate() {
            if e.1 {
                // current is right
                proof_index |= 1 << i;
                hash = H::compress(&e.0, &hash, i);
            } else {
                // current is left
                hash = H::compress(&hash, &e.0, i);
            }
        }

        if proof_index != index {
            return false;
        }

        hash == self.root_hash()
    }

    pub fn root_hash(&self) -> Vec<u8> {
        self.get_hash((0, 0))
    }

}