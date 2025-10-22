use crate::point::*;
#[cfg(feature = "serialize")]
use serde::{Deserialize, Serialize};
use slab::Slab;

pub struct BoundaryIter<'a, AData, BData> {
    qeds: &'a Qeds<AData, BData>,
    next_edge: Option<EdgeRefA<'a, AData, BData>>,
    first_edge_target: EdgeTarget,
}

impl<'a, AData, BData> BoundaryIter<'a, AData, BData> {
    pub(crate) fn new(qeds: &'a Qeds<AData, BData>, initial_boundary_edge: EdgeTarget) -> Self {
        Self {
            qeds,
            next_edge: None,
            first_edge_target: initial_boundary_edge,
        }
    }
}

impl<'a, AData, BData> Iterator for BoundaryIter<'a, AData, BData> {
    type Item = EdgeRefA<'a, AData, BData>;
    fn next(&mut self) -> Option<Self::Item> {
        let next_edge = self.next_edge;
        if let Some(next_edge) = next_edge {
            if next_edge.target() == self.first_edge_target {
                None
            } else {
                self.next_edge = Some(next_edge.r_prev());
                Some(next_edge)
            }
        } else {
            // This is the first time that next has been called.
            let next_edge = self.qeds.edge_a_ref(self.first_edge_target);
            self.next_edge = Some(next_edge.r_prev());
            Some(next_edge)
        }
    }
}

/// This data structure is a single instance of a quad-edge data structure.
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serialize", derive(Serialize, Deserialize))]
pub struct Qeds<AData, BData> {
    /// The vector of quads which we index into.
    pub quads: Slab<Quad<AData, BData>>,
}

impl<AData: Default, BData: Default> Qeds<AData, BData> {
    /// Create an edge in the [`Qeds`].
    pub fn make_edge(&mut self) -> EdgeRefA<'_, AData, BData> {
        let entry = self.quads.vacant_entry();
        let this_index = entry.key();
        let quad = Quad {
            edges_a: [
                // The base edge e.
                Edge {
                    next: EdgeTarget::new(this_index, 0, 0),
                    point: Default::default(),
                },
                // eSym
                Edge {
                    next: EdgeTarget::new(this_index, 2, 0),
                    point: Default::default(),
                },
            ],
            edges_b: [
                // eRot
                Edge {
                    next: EdgeTarget::new(this_index, 3, 0),
                    point: Default::default(),
                },
                // eSymRot
                Edge {
                    next: EdgeTarget::new(this_index, 1, 0),
                    point: Default::default(),
                },
            ],
        };
        entry.insert(quad);
        let t: &Qeds<AData, BData> = self;
        EdgeRefA {
            qeds: t,
            target: EdgeTarget::new(this_index, 0, 0),
        }
    }
}

impl<AData, BData: Default> Qeds<AData, BData> {
    /// Create an edge in the [`Qeds`].
    pub fn make_edge_with_a(&mut self, org: AData, dest: AData) -> EdgeRefA<'_, AData, BData> {
        let entry = self.quads.vacant_entry();
        let this_index = entry.key();
        let quad = Quad {
            edges_a: [
                // The base edge e.
                Edge {
                    next: EdgeTarget::new(this_index, 0, 0),
                    point: org,
                },
                // eSym
                Edge {
                    next: EdgeTarget::new(this_index, 2, 0),
                    point: dest,
                },
            ],
            edges_b: [
                // eRot
                Edge {
                    next: EdgeTarget::new(this_index, 3, 0),
                    point: Default::default(),
                },
                // eSymRot
                Edge {
                    next: EdgeTarget::new(this_index, 1, 0),
                    point: Default::default(),
                },
            ],
        };
        entry.insert(quad);
        let t: &Qeds<AData, BData> = self;
        EdgeRefA {
            qeds: t,
            target: EdgeTarget::new(this_index, 0, 0),
        }
    }
}

impl<AData, BData> Qeds<AData, BData> {
    /// Create an edge in the [`Qeds`].
    pub fn make_edge_with_ab(
        &mut self,
        a_org: AData,
        a_dest: AData,
        b_org: BData,
        b_dest: BData,
    ) -> EdgeRefA<'_, AData, BData> {
        let entry = self.quads.vacant_entry();
        let this_index = entry.key();
        let quad = Quad {
            edges_a: [
                // The base edge e.
                Edge {
                    next: EdgeTarget::new(this_index, 0, 0),
                    point: a_org,
                },
                // eSym
                Edge {
                    next: EdgeTarget::new(this_index, 2, 0),
                    point: a_dest,
                },
            ],
            edges_b: [
                // eRot
                Edge {
                    next: EdgeTarget::new(this_index, 3, 0),
                    point: b_org,
                },
                // eSymRot
                Edge {
                    next: EdgeTarget::new(this_index, 1, 0),
                    point: b_dest,
                },
            ],
        };
        entry.insert(quad);
        let t: &Qeds<AData, BData> = self;
        EdgeRefA {
            qeds: t,
            target: EdgeTarget::new(this_index, 0, 0),
        }
    }
}

impl<AData: Clone, BData: Default> Qeds<AData, BData> {
    /// Connect the Org of a with the Dest of b by creating a new edge. TODO: we
    /// need to special case infinite edges.
    pub fn connect(
        &mut self,
        edge_a: EdgeTarget,
        edge_b: EdgeTarget,
    ) -> EdgeRefA<'_, AData, BData> {
        // First, make the new edge.
        // Set the Org of e to the Dest of a
        let p1 = self.edge_a(edge_a.sym()).point.clone();
        // Set the Dest of e to the Org of b
        let p2 = self.edge_a(edge_b).point.clone();
        let q_target = self.make_edge_with_a(p1, p2).target;
        self.splice(q_target, self.edge_ref(edge_a).l_next().target);
        self.splice(q_target.sym(), edge_b);
        self.edge_a_ref(q_target)
    }
}

impl<AData, BData> Qeds<AData, BData> {
    /// Create the simplest [`Qeds`] with a single [`Edge`] and a single
    /// [`Face`]. Covers the whole sphere.
    pub fn new() -> Self {
        Self { quads: Slab::new() }
    }

    pub fn edge_a(&self, target: EdgeTarget) -> &Edge<AData> {
        if target.r == 0 || target.r == 2 {
            &self.quads[target.e].edges_a[(target.r / 2) as usize]
        } else {
            unreachable!()
        }
    }

    pub fn edge_b(&self, target: EdgeTarget) -> &Edge<BData> {
        if target.r == 0 || target.r == 2 {
            unreachable!()
        } else {
            &self.quads[target.e].edges_b[(target.r / 2) as usize]
        }
    }

    pub fn edge_ref(&self, target: EdgeTarget) -> EdgeRefAny<'_, AData, BData> {
        EdgeRefAny { qeds: self, target }
    }

    pub fn edge_a_ref(&self, target: EdgeTarget) -> EdgeRefA<'_, AData, BData> {
        EdgeRefA { qeds: self, target }
    }

    pub fn edge_b_ref(&self, target: EdgeTarget) -> EdgeRefB<'_, AData, BData> {
        EdgeRefB { qeds: self, target }
    }

    pub fn edge_mut(&mut self, target: EdgeTarget) -> EdgeABMut<'_, AData, BData> {
        if target.r == 0 || target.r == 2 {
            EdgeABMut::A(&mut self.quads[target.e].edges_a[(target.r / 2) as usize])
        } else {
            EdgeABMut::B(&mut self.quads[target.e].edges_b[(target.r / 2) as usize])
        }
    }

    pub fn edge_a_mut(&mut self, target: EdgeTarget) -> &mut Edge<AData> {
        if target.r == 0 || target.r == 2 {
            &mut self.quads[target.e].edges_a[(target.r / 2) as usize]
        } else {
            unreachable!()
        }
    }

    pub fn edge_b_mut(&mut self, target: EdgeTarget) -> &mut Edge<BData> {
        if target.r == 0 || target.r == 2 {
            unreachable!()
        } else {
            &mut self.quads[target.e].edges_b[(target.r / 2) as usize]
        }
    }

    pub fn splice(&mut self, edge_a: EdgeTarget, edge_b: EdgeTarget) {
        let alpha = self.edge_ref(edge_a).onext().rot().target;
        let beta = self.edge_ref(edge_b).onext().rot().target;

        // We want to swap aOnext with bOnext and αOnext with βONext
        let ta = self.edge_ref(edge_a).onext().target;
        let x = self.edge_ref(edge_b).onext().target;
        self.edge_mut(edge_a).set_next(x);
        self.edge_mut(edge_b).set_next(ta);

        let ta = self.edge_ref(alpha).onext().target;
        let y = self.edge_ref(beta).onext().target;
        self.edge_mut(alpha).set_next(y);
        self.edge_mut(beta).set_next(ta);
    }

    pub fn delete(&mut self, e: EdgeTarget) {
        let oprev = self.edge_ref(e).oprev().target();
        self.splice(e, oprev);
        let sym_oprev = self.edge_ref(e).sym().oprev().target();
        self.splice(e.sym(), sym_oprev);
        // TODO: verify that this is correct
        self.quads.remove(e.e);
    }

    pub fn base_edges(&self) -> BaseEdgeIter<'_, AData, BData> {
        BaseEdgeIter::new(self)
    }
}

impl<AData, BData> Default for Qeds<AData, BData> {
    fn default() -> Self {
        Self::new()
    }
}

pub struct BaseEdgeIter<'a, AData, BData> {
    qeds: &'a Qeds<AData, BData>,
    quad_iter: slab::Iter<'a, Quad<AData, BData>>,
}

impl<'a, AData, BData> BaseEdgeIter<'a, AData, BData> {
    fn new(qeds: &'a Qeds<AData, BData>) -> Self {
        Self {
            qeds,
            quad_iter: qeds.quads.iter(),
        }
    }
}

impl<'a, AData, BData> Iterator for BaseEdgeIter<'a, AData, BData> {
    type Item = EdgeRefA<'a, AData, BData>;
    fn next(&mut self) -> Option<Self::Item> {
        let quad = self.quad_iter.next()?;
        Some(self.qeds.edge_a_ref(EdgeTarget::new(quad.0, 0, 0)))
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd)]
#[cfg_attr(feature = "serialize", derive(Serialize, Deserialize))]
pub struct Quad<AData, BData> {
    pub edges_a: [Edge<AData>; 2],
    pub edges_b: [Edge<BData>; 2],
}

#[derive(Debug)]
pub struct EdgeRefAny<'a, AData, BData> {
    pub qeds: &'a Qeds<AData, BData>,
    pub target: EdgeTarget,
}

impl<'a, AData, BData> Clone for EdgeRefAny<'a, AData, BData> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<'a, AData, BData> Copy for EdgeRefAny<'a, AData, BData> {}

impl<'a, AData, BData> EdgeRefAny<'a, AData, BData> {
    pub fn qeds(&self) -> &Qeds<AData, BData> {
        self.qeds
    }
    pub fn target(&self) -> EdgeTarget {
        self.target
    }
    pub fn target_mut(&mut self) -> &mut EdgeTarget {
        &mut self.target
    }
    pub fn edge(&self) -> EdgeAB<'_, AData, BData> {
        // We know we can use this unsafe because lifetime gurantee that the
        // Qeds data structure has not been modified since this EdgeRef was
        // created.
        {
            if self.target().r % 2 == 0 {
                EdgeAB::A(self.qeds().edge_a(self.target()))
            } else {
                EdgeAB::B(self.qeds().edge_b(self.target()))
            }
        }
    }

    #[inline(always)]
    pub fn offset_r(&self, offset: u8) -> Self {
        let mut q = *self;
        let new_target = q.target().offset_r(offset);
        *q.target_mut() = new_target;
        q
    }

    #[inline(always)]
    pub fn onext(&self) -> Self {
        let edge = self.edge();
        let mut q = *self;
        *q.target_mut() = edge.next();
        q
    }

    /// Rot
    #[inline(always)]
    pub fn rot(&self) -> Self {
        self.offset_r(1)
    }

    #[inline(always)]
    pub fn sym(&self) -> Self {
        self.offset_r(2)
    }

    #[inline(always)]
    pub fn sym_rot(&self) -> Self {
        self.offset_r(3)
    }

    #[inline(always)]
    pub fn oprev(&self) -> Self {
        self.rot().onext().rot()
    }

    #[inline(always)]
    pub fn r_next(&self) -> Self {
        self.rot().onext().rot().rot().rot()
    }

    #[inline(always)]
    pub fn r_prev(&self) -> Self {
        self.sym().onext()
    }

    #[inline(always)]
    pub fn l_next(&self) -> Self {
        self.rot().rot().rot().onext().rot()
    }
}

#[derive(Debug)]
pub struct EdgeRefA<'a, AData, BData> {
    pub qeds: &'a Qeds<AData, BData>,
    pub target: EdgeTarget,
}

impl<'a, AData, BData> Clone for EdgeRefA<'a, AData, BData> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<'a, AData, BData> Copy for EdgeRefA<'a, AData, BData> {}

impl<'a, AData, BData> EdgeRefA<'a, AData, BData> {
    pub fn qeds(&self) -> &Qeds<AData, BData> {
        self.qeds
    }
    pub fn target(&self) -> EdgeTarget {
        self.target
    }
    pub fn target_mut(&mut self) -> &mut EdgeTarget {
        &mut self.target
    }
    pub fn edge(&self) -> &Edge<AData> {
        // We know we can use this unsafe because lifetime gurantee that the
        // Qeds data structure has not been modified since this EdgeRef was
        // created.
        self.qeds().edge_a(self.target())
    }
    // pub fn edge_mut(&mut self) -> &Edge<AData> {
    //     // We know we can use this unsafe because lifetime gurantee that the
    //     // Qeds data structure has not been modified since this EdgeRef was
    //     // created.
    //     unsafe { self.qeds().edge_a_mut(self.target()) }
    // }

    #[inline(always)]
    pub fn offset_r(&self, offset: u8) -> Self {
        let mut q = *self;
        let new_target = q.target().offset_r(offset);
        *q.target_mut() = new_target;
        q
    }

    #[inline(always)]
    pub fn onext(&self) -> Self {
        let edge = self.edge();
        let mut q = *self;
        *q.target_mut() = edge.next;
        q
    }

    #[inline(always)]
    pub fn d_prev(&self) -> Self {
        self.inv_rot().onext().inv_rot()
    }

    #[inline(always)]
    pub fn inv_rot(&self) -> EdgeRefB<'a, AData, BData> {
        self.rot().rot().rot()
    }

    /// Rot
    #[inline(always)]
    pub fn rot(&self) -> EdgeRefB<'a, AData, BData> {
        EdgeRefB {
            qeds: self.qeds,
            target: self.target().offset_r(1),
        }
    }

    #[inline(always)]
    pub fn sym(&self) -> Self {
        self.offset_r(2)
    }

    #[inline(always)]
    pub fn sym_rot(&self) -> EdgeRefB<'a, AData, BData> {
        EdgeRefB {
            qeds: self.qeds,
            target: self.target().offset_r(3),
        }
    }

    #[inline(always)]
    pub fn oprev(&self) -> Self {
        self.rot().onext().rot()
    }

    #[inline(always)]
    pub fn r_next(&self) -> Self {
        self.rot().onext().rot().rot().rot()
    }

    #[inline(always)]
    pub fn r_prev(&self) -> Self {
        self.sym().onext()
    }

    #[inline(always)]
    pub fn l_prev(&self) -> Self {
        self.onext().sym()
    }

    #[inline(always)]
    pub fn l_next(&self) -> Self {
        self.rot().rot().rot().onext().rot()
    }

    pub fn l_face(&self) -> Face<'_, AData, BData> {
        // Save the next edge of this face so that we know when to end the loop.
        let first_next_edge = self.l_next();
        let mut edges = Vec::new();
        // Loop around collecting the edges.
        let mut current_next_edge = self.l_next();
        loop {
            edges.push(current_next_edge);
            current_next_edge = current_next_edge.l_next();
            if current_next_edge.target() == first_next_edge.target() {
                break;
            }
        }
        Face { edges }
    }
    pub fn r_face(&self) -> Face<'_, AData, BData> {
        // Save the next edge of this face so that we know when to end the loop.
        let first_next_edge = self.r_next();
        let mut edges = Vec::new();
        // Loop around collecting the edges.
        let mut current_next_edge = self.r_next();
        loop {
            edges.push(current_next_edge);
            current_next_edge = current_next_edge.r_next();
            if current_next_edge.target() == first_next_edge.target() {
                break;
            }
        }
        Face { edges }
    }
}

#[derive(Debug)]
pub struct EdgeRefB<'a, AData, BData> {
    pub qeds: &'a Qeds<AData, BData>,
    pub target: EdgeTarget,
}

impl<'a, AData, BData> Clone for EdgeRefB<'a, AData, BData> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<'a, AData, BData> Copy for EdgeRefB<'a, AData, BData> {}

impl<'a, AData, BData> EdgeRefB<'a, AData, BData> {
    pub fn qeds(&self) -> &Qeds<AData, BData> {
        self.qeds
    }
    pub fn target(&self) -> EdgeTarget {
        self.target
    }
    pub fn target_mut(&mut self) -> &mut EdgeTarget {
        &mut self.target
    }

    pub fn edge(&self) -> &Edge<BData> {
        // We know we can use this unsafe because lifetime gurantee that the
        // Qeds data structure has not been modified since this EdgeRef was
        // created.
        self.qeds().edge_b(self.target())
    }

    #[inline(always)]
    pub fn offset_r(&self, offset: u8) -> Self {
        let mut q = *self;
        let new_target = q.target().offset_r(offset);
        *q.target_mut() = new_target;
        q
    }

    #[inline(always)]
    pub fn onext(&self) -> Self {
        let edge = self.edge();
        let mut q = *self;
        *q.target_mut() = edge.next;
        q
    }

    /// Rot
    #[inline(always)]
    pub fn rot(&self) -> EdgeRefA<'a, AData, BData> {
        EdgeRefA {
            qeds: self.qeds,
            target: self.target().offset_r(1),
        }
    }

    #[inline(always)]
    pub fn sym(&self) -> Self {
        self.offset_r(2)
    }

    #[inline(always)]
    pub fn sym_rot(&self) -> EdgeRefA<'a, AData, BData> {
        EdgeRefA {
            qeds: self.qeds,
            target: self.target().offset_r(3),
        }
    }

    #[inline(always)]
    pub fn oprev(&self) -> EdgeRefB<'a, AData, BData> {
        self.rot().onext().rot()
    }

    #[inline(always)]
    pub fn r_next(&self) -> Self {
        self.rot().onext().rot().rot().rot()
    }

    #[inline(always)]
    pub fn r_prev(&self) -> Self {
        self.sym().onext()
    }

    #[inline(always)]
    pub fn l_next(&self) -> Self {
        self.rot().rot().rot().onext().rot()
    }

    #[inline(always)]
    pub fn inv_rot(&self) -> EdgeRefA<'a, AData, BData> {
        self.rot().rot().rot()
    }

    pub fn l_face(&self) -> FaceB<'_, AData, BData> {
        // Save the next edge of this face so that we know when to end the loop.
        let first_next_edge = self.l_next();
        let mut edges = Vec::new();
        // Loop around collecting the edges.
        let mut current_next_edge = self.l_next();
        loop {
            edges.push(current_next_edge);
            current_next_edge = current_next_edge.l_next();
            if current_next_edge.target() == first_next_edge.target() {
                break;
            }
        }
        FaceB { edges }
    }
    pub fn r_face(&self) -> FaceB<'_, AData, BData> {
        // Save the next edge of this face so that we know when to end the loop.
        let first_next_edge = self.r_next();
        let mut edges = Vec::new();
        // Loop around collecting the edges.
        let mut current_next_edge = self.r_next();
        loop {
            edges.push(current_next_edge);
            current_next_edge = current_next_edge.r_next();
            if current_next_edge.target() == first_next_edge.target() {
                break;
            }
        }
        FaceB { edges }
    }
}

#[derive(Copy, Clone, Debug, Ord, PartialOrd, Eq, PartialEq)]
pub enum Direction {
    Left,
    Straight,
    Right,
}

/// Does p3 lie to the left or the right (or is collinear) of the line formed by
/// p1 and p2.
pub fn left_or_right(pa: Point, pb: Point, pc: Point) -> Direction {
    let r = robust::orient2d(pa.into(), pb.into(), pc.into());
    if !r.is_finite() {
        panic!("Non-finite determinant");
    } else if r > 0.0 {
        Direction::Left
    } else if r == 0.0 {
        Direction::Straight
    } else {
        Direction::Right
    }
}

impl<AData, BData> PartialEq for EdgeRefA<'_, AData, BData> {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.qeds, other.qeds) && self.target == other.target
    }
}

impl<AData, BData> Eq for EdgeRefA<'_, AData, BData> {}

impl<AData, BData> PartialEq for EdgeRefB<'_, AData, BData> {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.qeds, other.qeds) && self.target == other.target
    }
}

impl<AData, BData> Eq for EdgeRefB<'_, AData, BData> {}

// r and f can fit in one byte. If we limit the size of qeds we could fit it in
// the edge index.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serialize", derive(Serialize, Deserialize))]
pub struct EdgeTarget {
    // The index of the Quad in the Qeds.
    pub e: usize,
    // Can only be 0, 1, 2, or 3.
    pub r: u8,
    // Can only be 0 or 1.
    pub f: u8,
}

impl EdgeTarget {
    pub fn new(e: usize, r: u8, f: u8) -> Self {
        Self { e, r, f }
    }

    #[inline(always)]
    fn offset_r(&self, offset: u8) -> Self {
        let mut q = *self;
        q.r = (q.r + offset) % 4;
        q
    }

    /// Rot
    #[inline(always)]
    pub fn rot(&self) -> Self {
        self.offset_r(1)
    }

    #[inline(always)]
    pub fn sym(&self) -> Self {
        self.offset_r(2)
    }

    #[inline(always)]
    pub fn sym_rot(&self) -> Self {
        self.offset_r(3)
    }
}

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Vertex {
    /// The index of the edge which starts at this vertex.
    pub edge_index: usize,
    pub point: Point,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd)]
pub enum EdgeAB<'a, AData, BData> {
    A(&'a Edge<AData>),
    B(&'a Edge<BData>),
}

impl<'a, AData, BData> EdgeAB<'a, AData, BData> {
    fn next(&self) -> EdgeTarget {
        match *self {
            EdgeAB::A(edge) => edge.next,
            EdgeAB::B(edge) => edge.next,
        }
    }
}

#[derive(Debug, PartialEq, Eq, PartialOrd)]
pub enum EdgeABMut<'a, AData, BData> {
    A(&'a mut Edge<AData>),
    B(&'a mut Edge<BData>),
}

impl<'a, AData, BData> EdgeABMut<'a, AData, BData> {
    fn set_next(&mut self, target: EdgeTarget) {
        match self {
            EdgeABMut::A(edge) => edge.next = target,
            EdgeABMut::B(edge) => edge.next = target,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd)]
#[cfg_attr(feature = "serialize", derive(Serialize, Deserialize))]
pub struct Edge<Data> {
    pub next: EdgeTarget,
    pub point: Data,
}

#[derive(Clone, Debug)]
pub struct Face<'a, AData, BData> {
    /// An index to any of the CCW oriented edges surrounding this face.
    pub edges: Vec<EdgeRefA<'a, AData, BData>>,
}

pub type VertexPairIter<'a, AData, BData> = std::iter::Zip<
    FaceVerticesIter<'a, AData, BData>,
    std::iter::Skip<std::iter::Cycle<FaceVerticesIter<'a, AData, BData>>>,
>;

#[derive(Clone, Debug)]
pub struct FaceB<'a, AData, BData> {
    /// An index to any of the CCW oriented edges surrounding this face.
    pub edges: Vec<EdgeRefB<'a, AData, BData>>,
}

#[derive(Clone, Debug)]
pub struct FaceVerticesIter<'a, AData, BData> {
    pub(crate) face: &'a Face<'a, AData, BData>,
    pub(crate) next_index: usize,
}

impl<'a, AData, BData> FaceVerticesIter<'a, AData, BData> {
    pub(crate) fn new(face: &'a Face<'a, AData, BData>) -> Self {
        Self {
            face,
            next_index: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn targetting() {
        let t = EdgeTarget::new(0, 0, 0);
        assert_eq!(t.offset_r(2), EdgeTarget::new(0, 2, 0));
    }

    // First we test the Edge Function Properties as described by G&S 2.3
    #[test]
    fn e1() {
        // Step 1. Create a Qeds data structure.
        let mut qeds: Qeds<(), ()> = Qeds::new();
        // Step 2. Add some edges to it.
        let e = qeds.make_edge();
        // e with Rot applied 4 times is the same edge as e
        assert_eq!(e.rot().rot().rot().rot(), e);
    }

    /// Test the basic relationships of a quad that should always be true.
    #[test]
    fn basic_relationships() {
        // Step 1. Create a Qeds data structure.
        let mut qeds: Qeds<(), ()> = Qeds::new();
        // Step 2. Add some edges to it.
        let e = qeds.make_edge();
        // Step 3. Get the single edge in the data structure.
        // e with Rot applied 2 times is not the same edge as e
        assert_eq!(
            e.rot().rot().edge() as *const Edge<()>,
            e.sym().edge() as *const Edge<()>
        );
        assert_eq!(
            e.rot().sym().edge() as *const Edge<()>,
            e.sym().rot().edge() as *const Edge<()>
        );
        assert_eq!(
            e.onext().edge() as *const Edge<()>,
            e.edge() as *const Edge<()>
        );
    }

    /// Test the basic relationships of a quad that should always be true for an
    /// unconnected edge.
    #[test]
    fn floating_relationships() {
        // Step 1. Create a Qeds data structure.
        let mut qeds: Qeds<(), ()> = Qeds::new();
        // Step 2. Add some edges to it.
        let e = qeds.make_edge();
        // Step 3. Get the single edge in the data structure.
        // e with Rot applied 2 times is not the same edge as e
        assert_eq!(
            e.rot().rot().edge() as *const Edge<()>,
            e.sym().edge() as *const Edge<()>
        );
        assert_eq!(
            e.rot().sym().edge() as *const Edge<()>,
            e.sym().rot().edge() as *const Edge<()>
        );
        assert_eq!(
            e.onext().edge() as *const Edge<()>,
            e.edge() as *const Edge<()>
        );
        assert_eq!(
            e.sym().onext().edge() as *const Edge<()>,
            e.sym().edge() as *const Edge<()>
        );
        assert_eq!(
            e.rot().onext().edge() as *const Edge<()>,
            e.rot().sym().edge() as *const Edge<()>
        );
    }

    #[test]
    fn rot2_eq_sym() {
        // Step 1. Create a Qeds data structure.
        let mut qeds: Qeds<(), ()> = Qeds::new();
        // Step 2. Add some edges to it.
        let q1 = qeds.make_edge().target;
        // Step 3. Get the single edge in the data structure.
        // e with Rot applied 2 times is not the same edge as e
        assert_eq!(q1.rot().rot(), q1.sym());
    }

    #[test]
    fn d5() {
        let mut qeds = Qeds::new();
        let q1 = qeds.make_edge().target;
        let q1 = { qeds.edge_a_ref(q1) };
        assert_eq!(q1.rot().edge(), &qeds.quads[0].edges_b[0]);
        assert_eq!(q1.sym().edge(), &qeds.quads[0].edges_a[1]);
        assert_eq!(q1.sym().rot().edge(), &qeds.quads[0].edges_b[1]);
        assert_eq!(
            q1.rot().sym().edge() as *const Edge<()>,
            &qeds.quads[0].edges_b[1] as *const Edge<()>
        );
        assert_eq!(
            q1.rot().onext().edge() as *const Edge<()>,
            q1.rot().sym().edge() as *const Edge<()>
        );
        assert_eq!(
            q1.rot().sym().onext().edge() as *const Edge<()>,
            q1.rot().edge() as *const Edge<()>
        );
        assert_eq!(
            q1.onext().edge() as *const Edge<()>,
            &qeds.quads[0].edges_a[0] as *const Edge<()>
        );
    }

    #[test]
    fn splice_test() {
        // Step 1. Create a Qeds data structure.
        let mut qeds: Qeds<(), ()> = Qeds::new();
        // Step 2. Add the first two edges. Converting them into references as
        // it is the more useful format for testing.
        let q1 = qeds.make_edge().target;
        let q2 = qeds.make_edge().target;
        {
            qeds.splice(q1.sym(), q2);
            let q1 = qeds.edge_a_ref(q1);
            let q2 = qeds.edge_a_ref(q2);
            // 1. aSymOnext == b
            assert_eq!(q1.sym().onext(), q2);
            // 2. bOnext == aSym
            assert_eq!(q2.onext(), q1.sym());
            // 3. aRotOnext == bRot
            assert_ne!(q1.rot().onext(), q1.rot());
            // 4. bRotOnext == aRot
            assert_eq!(q2.rot().onext(), q1.rot());

            // bLnext = bSym
            assert_eq!(q2.l_next(), q2.sym());
        }
    }

    #[test]
    fn single_splice() {
        // Step 1. Create a Qeds data structure.
        let mut qeds: Qeds<(), ()> = Qeds::new();
        // Step 2. Add some edges to it.
        let a = qeds.make_edge().target;
        let b = qeds.make_edge().target;
        // Step 3. Splice those edges together so that we actually have
        // something of a network.
        {
            qeds.splice(a.sym(), b);
        }
    }
}
