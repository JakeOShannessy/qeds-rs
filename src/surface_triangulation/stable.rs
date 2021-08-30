//! A surface triangulation is a triangulation with boundaries and no
//! retriangulation.
use super::SurfaceTriangulation;
use crate::point::*;
use crate::qeds::*;
#[cfg(serialize)]
use serde::Serialize;

#[derive(Debug)]
/// A Qeds data structure specialised to a 2d triangulation.
pub struct SurfaceTriangulationStable<'a, T> {
    pub st: &'a mut SurfaceTriangulation<T>,
}

impl<'a, T> Drop for SurfaceTriangulationStable<'a, T> {
    fn drop(&mut self) {
        self.st.retriangulate_all();
    }
}

impl<'a, T> SurfaceTriangulationStable<'a, T> {
    pub fn base_targets(&self) -> impl Iterator<Item = EdgeTarget> + '_ {
        self.st.base_targets()
    }
}

impl<'a, T: Clone + Default> SurfaceTriangulationStable<'a, T> {
    /// Add a point to a specified edge. If the point lies on one of the
    /// vertices just add it there.
    pub fn add_point_to_edge(
        &mut self,
        edge_target: EdgeTarget,
        point: Point,
        data: T,
    ) -> Vec<EdgeTarget> {
        self.st
            .add_point_to_edge_impl(edge_target, point, data, false)
    }
}
