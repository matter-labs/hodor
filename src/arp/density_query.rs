use crate::air::{DenseConstraint, SparseConstraint};

pub(crate) trait DensityQuery: Send + Sync + 'static {
    fn next_row(&mut self) -> Option<usize>;
}

impl Iterator for DensityQuery {
    type Item = usize;
    
    fn next(&mut self) -> Option<usize> {
        self.next_row()
    }
}

pub(crate) struct DenseConstraintQuery {
    current_row: usize,
    limit: usize
}

impl DensityQuery for DenseConstraintQuery {
    fn next_row(&mut self) -> Option<usize> {
        if self.current_row < self.limit {
            let t = self.current_row;
            self.current_row += 1;

            Some(t)
        } else {
            None
        }
    }
}

impl DenseConstraintQuery {
    pub(crate) fn new(dense: &DenseConstraint, num_rows: usize) -> Self {
        let start_at = dense.start_at;
        let span = dense.span;

        let limit = num_rows - span;

        Self {
            current_row: start_at,
            limit: limit
        }
    }
}