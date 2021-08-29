//! A surface triangulation is a triangulation with boundaries and no
//! retriangulation.
use crate::point::*;
use crate::qeds::*;
use serde::{Deserialize, Serialize};

use super::SurfaceTriangulation;

#[derive(Clone, Debug, Serialize, Deserialize)]
/// A Qeds data structure specialised to a 2d triangulation.
pub struct SurfaceTriangulationStable<T> {
    pub st: SurfaceTriangulation<T>,
}

impl<T> SurfaceTriangulationStable<T> {
    pub fn unfreeze(mut self) -> SurfaceTriangulation<T> {
        self.st.retriangulate_all();
        self.st
    }
    pub fn base_targets(&self) -> impl Iterator<Item = EdgeTarget> + '_ {
        self.st.base_targets()
    }
}

impl<T: Clone + Default + Serialize> SurfaceTriangulationStable<T> {
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
