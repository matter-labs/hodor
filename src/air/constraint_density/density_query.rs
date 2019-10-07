pub trait DensityQuery: Send + Sync + 'static {
    fn next_row(&mut self) -> Option<usize>;
}

impl Iterator for dyn DensityQuery {
    type Item = usize;
    
    fn next(&mut self) -> Option<usize> {
        self.next_row()
    }
}