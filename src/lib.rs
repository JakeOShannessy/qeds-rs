use core::ops::Add;
use core::ops::Sub;
use std::{ops::Mul};
use slab::Slab;

pub struct Triangulation {
    pub qeds: Qeds<Box<Option<Point>>,()>,
    pub free_point: Option<Point>,
    pub boundary_edge: Option<EdgeTarget>,
}

impl Triangulation where {
    pub fn new() -> Self {
        Self {
            qeds: Qeds::new(),
            free_point: None,
            boundary_edge: None,
        }
    }

    pub fn add_point(&mut self, point: Point) {
        // TODO: consider: points lying on other points, on lines, and within triangles.
        if self.qeds.quads.len() == 0 {
            // If there are no points, we can just add a free-floating point.
            if let Some(free_point) = self.free_point {
                let e = self.qeds.make_edge().target();
                // Assign the appropriate geometry to the edge.
                unsafe {
                    self.qeds.edge_a_mut(e).point = Box::new(Some(free_point));
                    self.qeds.edge_a_mut(e.sym()).point = Box::new(Some(point));
                }
                self.free_point = None;
                self.boundary_edge = Some(e);
            } else {
                self.free_point = Some(point);
            }
        } else {
            // TODO: do some proper locating. We iterate through each of the
            // boundary edges and connect the the point to all edges it lies to
            // the right of, this should maintain convecity.
            let mut relevant_edges = Vec::new();
            for edge in self.boundary().unwrap() {
                if edge.lies_right(point) {
                    relevant_edges.push(edge.target());
                }
            }
            // The edges need to be sorted such that their points are in CCW
            // order around the new point. We can do this by performing a CCW on
            // the base edge points around the new point.
            unsafe {
                self.add_external_point_unordered(relevant_edges, point);
            }
            // todo!()
        }
    }

    /// Assumes the EdgeTargets have already been put in the correct order.
    pub unsafe fn add_external_point_ordered(
        &mut self,
        boundary_edges: Vec<EdgeTarget>,
        point: Point,
    ) {
        let mut edges = boundary_edges.into_iter();
        let (g,e) = {
            let (a,b) = self
            .add_external_point(edges.next().unwrap(), point);
            (a.target(),b.target())
        };
        for f in edges {
            self.qeds.connect(e.sym(), f.sym());
        }
        self.boundary_edge = Some(g);
    }

    pub unsafe fn add_external_point_unordered(
        &mut self,
        mut boundary_edges: Vec<EdgeTarget>,
        point: Point,
    ) {
        // Reorder the edges so that their base points are CCW around the point.
        boundary_edges.sort_by(|a, b| {
            let pa = self.qeds.edge_a_ref(*a).edge().point.unwrap();
            let pb = self.qeds.edge_a_ref(*b).edge().point.unwrap();
            match left_or_right(point, pb, pa) {
                Direction::Left => std::cmp::Ordering::Less,
                Direction::Straight => std::cmp::Ordering::Equal,
                Direction::Right => std::cmp::Ordering::Greater,
            }
        });
        self.add_external_point_ordered(boundary_edges, point);
    }

    pub fn boundary(&self) -> Option<BoundaryIter<Box<Option<Point>>,()>> {
        if let Some(boundary_edge) = self.boundary_edge {
            Some(BoundaryIter::new(&self.qeds, boundary_edge))
        } else {
            None
        }
    }

    pub unsafe fn add_external_point(
        &mut self,
        boundary_edge: EdgeTarget,
        p: Point,
    ) -> (EdgeRefA<Box<Option<Point>>,()>, EdgeRefA<Box<Option<Point>>,()>) {
        let e = self.qeds.make_edge().target();
        self.qeds.edge_a_mut(e).point = Box::new(Some(p));
        self.qeds.edge_a_mut(e.sym()).point = self.qeds.edge_a_ref(boundary_edge).sym().edge().point.clone();
        self.qeds.splice(e.sym(), boundary_edge.sym());
        let f = self.qeds.connect(boundary_edge.sym(), e).target();
        (self.qeds.edge_a_ref(f), self.qeds.edge_a_ref(e))
    }
}

pub struct BoundaryIter<'a,AData,BData> {
    qeds: &'a Qeds<AData,BData>,
    next_edge: Option<EdgeRefA<'a,AData,BData>>,
    first_edge_target: EdgeTarget,
}

impl<'a,AData,BData> BoundaryIter<'a,AData,BData> {
    fn new(qeds: &'a Qeds<AData,BData>, initial_boundary_edge: EdgeTarget) -> Self {
        // let next_edge = unsafe {
        //     qeds.edge_a_ref(initial_boundary_edge)
        // };
        Self {
            qeds: qeds,
            next_edge: None,
            first_edge_target: initial_boundary_edge,
        }
    }
}

impl<'a,AData,BData> Iterator for BoundaryIter<'a,AData,BData> {
    type Item = EdgeRefA<'a,AData,BData>;
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
            let next_edge = unsafe { self.qeds.edge_a_ref(self.first_edge_target) };
            self.next_edge = Some(next_edge.r_prev());
            Some(next_edge)
        }
    }
}

/// This data structure is a single instance of a quad-edge data structure.
#[derive(Clone, Debug)]
pub struct Qeds<AData,BData> {
    /// The vector of quads which we index into.
    pub quads: Slab<Quad<AData,BData>>,
}

impl<AData:Default,BData:Default> Qeds<AData,BData> {
    /// Create an edge in the [`Qeds`].
    pub fn make_edge(&mut self) -> EdgeRefA<AData,BData> {
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
        let t: &Qeds<AData,BData> = self;
        EdgeRefA {
            qeds: t,
            target: EdgeTarget::new(this_index, 0, 0),
        }
    }

}

impl<AData:Default + Clone,BData:Default> Qeds<AData,BData> {
    /// Connect the Org of a with the Dest of b by creating a new edge.
    pub fn connect(&mut self, edge_a: EdgeTarget, edge_b: EdgeTarget) -> EdgeRefA<AData,BData> {
        // First, make the new edge.
        let q_target = self.make_edge().target;
        unsafe {
            // Set the Org of e to the Dest of a
            self.edge_a_mut(q_target).point = self.edge_a(edge_a.sym()).point.clone();
            // Set the Dest of e to the Org of b
            self.edge_a_mut(q_target.sym()).point = self.edge_a(edge_b).point.clone();

            self.splice(q_target, self.edge_ref(edge_a).l_next().target);
            self.splice(q_target.sym(), edge_b);
            self.edge_a_ref(q_target)
        }
    }
}

impl<AData,BData> Qeds<AData,BData> {
    /// Create the simplest [`Qeds`] with a single [`Edge`] and a single
    /// [`Face`]. Covers the whole sphere.
    pub fn new() -> Self {
        Self {
            quads: Slab::new(),
        }
    }

    pub unsafe fn edge_a(&self, target: EdgeTarget) -> &Edge<AData> {
        if target.r == 0 || target.r == 2 {
            &self.quads[target.e].edges_a[(target.r/2) as usize]
        } else {
            unreachable!()
        }
    }

    pub unsafe fn edge_b(&self, target: EdgeTarget) -> &Edge<BData> {
        if target.r == 0 || target.r == 2 {
            unreachable!()
        } else {
            &self.quads[target.e].edges_b[(target.r/2) as usize]
        }
    }

    pub unsafe fn edge_ref(&self, target: EdgeTarget) -> EdgeRefAny<AData,BData> {
        EdgeRefAny {
            qeds: &self,
            target,
        }
    }

    pub unsafe fn edge_a_ref(&self, target: EdgeTarget) -> EdgeRefA<AData,BData> {
        EdgeRefA {
            qeds: &self,
            target,
        }
    }

    pub unsafe fn edge_mut(&mut self, target: EdgeTarget) -> EdgeABMut<AData,BData> {
        if target.r == 0 || target.r == 2 {
            EdgeABMut::A(&mut self.quads[target.e].edges_a[(target.r/2) as usize])
        } else {
            EdgeABMut::B(&mut self.quads[target.e].edges_b[(target.r/2) as usize])
        }
    }

    pub unsafe fn edge_a_mut(&mut self, target: EdgeTarget) -> &mut Edge<AData> {
        if target.r == 0 || target.r == 2 {
            &mut self.quads[target.e].edges_a[(target.r/2) as usize]
        } else {
            unreachable!()
        }
    }

    pub unsafe fn edge_b_mut(&mut self, target: EdgeTarget) -> &mut Edge<BData> {
        if target.r == 0 || target.r == 2 {
            unreachable!()
        } else {
            &mut self.quads[target.e].edges_b[(target.r/2) as usize]
        }
    }

    pub unsafe fn splice(&mut self, edge_a: EdgeTarget, edge_b: EdgeTarget) {
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

    // pub unsafe fn delete(&mut self, e: &mut Edge) {
    //     // TODO: we don't actually free the memory here.
    //     self.splice(e, e.oprev_mut());
    //     self.splice(e.sym_mut(), e.sym_mut().oprev_mut());
    // }

    pub fn base_edges(&self) -> BaseEdgeIter<AData,BData> {
        BaseEdgeIter::new(self)
    }
}

pub struct BaseEdgeIter<'a,AData,BData> {
    qeds: &'a Qeds<AData,BData>,
    quad_iter: slab::Iter<'a,Quad<AData,BData>>,
}

impl<'a,AData,BData> BaseEdgeIter<'a,AData,BData> {
    fn new(qeds: &'a Qeds<AData,BData>) -> Self {
        Self {
            qeds: qeds,
            quad_iter: qeds.quads.iter(),
        }
    }
}

impl<'a,AData,BData> Iterator for BaseEdgeIter<'a,AData,BData> {
    type Item = EdgeRefA<'a,AData,BData>;
    fn next(&mut self) -> Option<Self::Item> {
        let quad = self.quad_iter.next()?;
        unsafe {
            Some(self.qeds.edge_a_ref(EdgeTarget::new(quad.0,0,0)))
        }
    }
}

#[derive(Clone, Debug, PartialEq, PartialOrd)]
pub struct Quad<AData,BData> {
    pub edges_a: [Edge<AData>; 2],
    pub edges_b: [Edge<BData>; 2],
}


#[derive(Debug)]
pub struct EdgeRefAny<'a,AData,BData> {
    pub qeds: &'a Qeds<AData,BData>,
    pub target: EdgeTarget,
}

impl<'a,AData,BData> Clone for EdgeRefAny<'a,AData,BData> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<'a,AData,BData> Copy for EdgeRefAny<'a,AData,BData> {}

impl<'a,AData,BData> EdgeRefAny<'a,AData,BData> {
    pub fn qeds(&self) -> &Qeds<AData,BData> {
        self.qeds
    }
    pub fn target(&self) -> EdgeTarget {
        self.target
    }
    pub fn target_mut(&mut self) -> &mut EdgeTarget {
        &mut self.target
    }
    pub fn edge(&self) -> EdgeAB<AData,BData> {
        // We know we can use this unsafe because lifetime gurantee that the
        // Qeds data structure has not been modified since this EdgeRef was
        // created.
        unsafe {
            if self.target().r % 2 == 0 {
                EdgeAB::A(self.qeds().edge_a(self.target()))
            } else {
                EdgeAB::B(self.qeds().edge_b(self.target()))
            }
        }

    }

    #[inline(always)]
    pub fn offset_r(&self, offset: u8) -> Self {
        let mut q = self.clone();
        let new_target = q.target().offset_r(offset);
        *q.target_mut() = new_target;
        q
    }

    #[inline(always)]
    pub fn onext(&self) -> Self {
        let edge = self.edge();
        let mut q = self.clone();
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
pub struct EdgeRefA<'a,AData,BData> {
    pub qeds: &'a Qeds<AData,BData>,
    pub target: EdgeTarget,
}

impl<'a,AData,BData> Clone for EdgeRefA<'a,AData,BData> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<'a,AData,BData> Copy for EdgeRefA<'a,AData,BData> {}

impl<'a,BData> EdgeRefA<'a,Box<Option<Point>>,BData> {
    pub fn midpoint(&self) -> Option<Point> {
        let this_point = (*self.edge().point)?;
        let next_point = (*self.sym().edge().point)?;
        Some(this_point.midpoint(next_point))
    }
    pub fn lies_right(&self, point: Point) -> bool {
        let pa = if let Some(p) = *self.edge().point {
            p
        } else {
            return false;
        };
        let pb = if let Some(p) = *self.sym().edge().point {
            p
        } else {
            return false;
        };
        left_or_right(pa, pb, point) == Direction::Right
    }
}

impl<'a,AData,BData> EdgeRefA<'a,AData,BData> {
    pub fn qeds(&self) -> &Qeds<AData,BData> {
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
        unsafe { self.qeds().edge_a(self.target()) }
    }

    #[inline(always)]
    pub fn offset_r(&self, offset: u8) -> Self {
        let mut q = self.clone();
        let new_target = q.target().offset_r(offset);
        *q.target_mut() = new_target;
        q
    }

    #[inline(always)]
    pub fn onext(&self) -> Self {
        let edge = self.edge();
        let mut q = self.clone();
        *q.target_mut() = edge.next;
        q
    }

    /// Rot
    #[inline(always)]
    pub fn rot(&self) -> EdgeRefB<'a,AData,BData> {
        EdgeRefB {
            qeds: self.qeds,
            target: self.target().offset_r(1)
        }
    }

    #[inline(always)]
    pub fn sym(&self) -> Self {
        self.offset_r(2)
    }

    #[inline(always)]
    pub fn sym_rot(&self) -> EdgeRefB<'a,AData,BData> {
        EdgeRefB {
            qeds: self.qeds,
            target: self.target().offset_r(3)
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
    pub fn l_next(&self) -> Self {
        self.rot().rot().rot().onext().rot()
    }

    pub fn l_face(&self) -> Face<AData,BData> {
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
    pub fn r_face(&self) -> Face<AData,BData> {
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
pub struct EdgeRefB<'a,AData,BData> {
    pub qeds: &'a Qeds<AData,BData>,
    pub target: EdgeTarget,
}

impl<'a,AData,BData> Clone for EdgeRefB<'a,AData,BData> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<'a,AData,BData> Copy for EdgeRefB<'a,AData,BData> {}

impl<'a,AData,BData> EdgeRefB<'a,AData,BData> {
    pub fn qeds(&self) -> &Qeds<AData,BData> {
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
        unsafe { self.qeds().edge_b(self.target()) }
    }

    #[inline(always)]
    pub fn offset_r(&self, offset: u8) -> Self {
        let mut q = self.clone();
        let new_target = q.target().offset_r(offset);
        *q.target_mut() = new_target;
        q
    }

    #[inline(always)]
    pub fn onext(&self) -> Self {
        let edge = self.edge();
        let mut q = self.clone();
        *q.target_mut() = edge.next;
        q
    }

    /// Rot
    #[inline(always)]
    pub fn rot(&self) -> EdgeRefA<'a,AData,BData> {
        EdgeRefA {
            qeds: self.qeds,
            target: self.target().offset_r(1)
        }
    }

    #[inline(always)]
    pub fn sym(&self) -> Self {
        self.offset_r(2)
    }

    #[inline(always)]
    pub fn sym_rot(&self) -> EdgeRefA<'a,AData,BData> {
        EdgeRefA {
            qeds: self.qeds,
            target: self.target().offset_r(3)
        }
    }

    #[inline(always)]
    pub fn oprev(&self) -> EdgeRefB<'a,AData,BData> {
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

    pub fn l_face(&self) -> FaceB<AData,BData> {
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
    pub fn r_face(&self) -> FaceB<AData,BData> {
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

pub fn left_or_right(
    Point { x: x1, y: y1 }: Point,
    Point { x: x2, y: y2 }: Point,
    Point { x: x3, y: y3 }: Point,
) -> Direction {
    use Direction::*;
    let determinant = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
    if determinant > 0.0 {
        Left
    } else if determinant == 0.0 {
        Straight
    } else {
        Right
    }
}

impl<AData,BData> PartialEq for EdgeRefA<'_,AData,BData> {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.qeds, other.qeds) && self.target == other.target
    }
}

impl<AData,BData> Eq for EdgeRefA<'_,AData,BData> {}

impl<AData,BData> PartialEq for EdgeRefB<'_,AData,BData> {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.qeds, other.qeds) && self.target == other.target
    }
}

impl<AData,BData> Eq for EdgeRefB<'_,AData,BData> {}

// r and f can fit in one byte. If we limit the size of qeds we could fit it in
// the edge index.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
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
        let mut q = self.clone();
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

#[derive(Clone, Debug, PartialEq, PartialOrd)]
pub enum EdgeAB<'a,AData,BData> {
    A(&'a Edge<AData>),
    B(&'a Edge<BData>),
}

impl<'a,AData,BData> EdgeAB<'a,AData,BData> {
    fn next(&self) -> EdgeTarget {
        match *self {
            EdgeAB::A(edge) => edge.next,
            EdgeAB::B(edge) => edge.next,
        }
    }
}


#[derive(Debug, PartialEq, PartialOrd)]
pub enum EdgeABMut<'a,AData,BData> {
    A(&'a mut Edge<AData>),
    B(&'a mut Edge<BData>),
}

impl<'a,AData,BData> EdgeABMut<'a,AData,BData> {
    // fn next(&self) -> EdgeTarget {
    //     match self {
    //         EdgeABMut::A(edge) => edge.next,
    //         EdgeABMut::B(edge) => edge.next,
    //     }
    // }
    fn set_next(&mut self, target: EdgeTarget) {
        match self {
            EdgeABMut::A(edge) => edge.next = target,
            EdgeABMut::B(edge) => edge.next = target,
        }
    }
}

#[derive(Clone, Debug, PartialEq, PartialOrd)]
pub struct Edge<Data> {
    pub next: EdgeTarget,
    pub point: Data,
}

#[derive(Clone, Debug)]
pub struct Face<'a,AData,BData> {
    /// An index to any of the CCW oriented edges surrounding this face.
    pub edges: Vec<EdgeRefA<'a,AData,BData>>,
}

impl<BData:Clone> Face<'_,Box<Option<Point>>,BData> {
    pub fn centroid(&self) -> Option<Point> {
        let signed_area = self.signed_area()?;
        if signed_area <= 0.0 {
            return None;
        }
        let mut xsum = 0.0;
        let mut ysum = 0.0;
        for (pa, pb) in self.vertices_pairs() {
            let pa = pa?;
            let pb = pb?;
            xsum += (pa.x + pb.x) * (pa.x * pb.y - pb.x * pa.y);
            ysum += (pa.y + pb.y) * (pa.x * pb.y - pb.x * pa.y);
        }
        Some(Point {
            x: xsum / 6.0 / signed_area,
            y: ysum / 6.0 / signed_area,
        })
    }

    pub fn signed_area(&self) -> Option<f64> {
        let mut sum = 0.0;
        for (pa, pb) in self.vertices_pairs() {
            let pa = pa?;
            let pb = pb?;
            sum += pa.x * pb.y - pb.x * pa.y;
        }
        Some(0.5 * sum)
    }
    pub fn vertices(&self) -> FaceVerticesIter<BData> {
        FaceVerticesIter::new(self)
    }
    pub fn vertices_pairs(
        &self,
    ) -> std::iter::Zip<FaceVerticesIter<'_,BData>, std::iter::Skip<std::iter::Cycle<FaceVerticesIter<'_,BData>>>>
    {
        self.vertices().zip(self.vertices().cycle().skip(1))
    }
}


#[derive(Clone, Debug)]
pub struct FaceB<'a,AData,BData> {
    /// An index to any of the CCW oriented edges surrounding this face.
    pub edges: Vec<EdgeRefB<'a,AData,BData>>,
}

#[derive(Clone, Debug)]
pub struct FaceVerticesIter<'a,BData> {
    face: &'a Face<'a,Box<Option<Point>>,BData>,
    next_index: usize,
}

impl<'a,BData> FaceVerticesIter<'a,BData> {
    fn new(face: &'a Face<'a,Box<Option<Point>>,BData>) -> Self {
        Self {
            face,
            next_index: 0,
        }
    }
}

impl<'a,BData> Iterator for FaceVerticesIter<'a,BData> {
    type Item = Option<Point>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.next_index >= self.face.edges.len() {
            None
        } else {
            let edge = self.face.edges.get(self.next_index);
            self.next_index += 1;
            Some(edge.and_then(|edge| *edge.edge().point))
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl Point {
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }
    pub fn midpoint(self, other: Self) -> Self {
        Self {
            x: (self.x + other.x) / 2.0,
            y: (self.y + other.y) / 2.0,
        }
    }

    pub fn is_finite(&self) -> bool {
        self.x.is_finite() && self.y.is_finite()
    }
}

impl Add for Point {
    type Output = Self;

    fn add(self, other: Point) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl Mul<f64> for Point {
    type Output = Self;

    fn mul(self, other: f64) -> Self {
        Self {
            x: self.x * other,
            y: self.y * other,
        }
    }
}

impl Sub for Point {
    type Output = Self;

    fn sub(self, other: Point) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
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
        let mut qeds: Qeds<(),()> = Qeds::new();
        // Step 2. Add some edges to it.
        let e = qeds.make_edge();
        // e with Rot applied 4 times is the same edge as e
        assert_eq!(e.rot().rot().rot().rot(), e);
    }

    /// Test the basic relationships of a quad that should always be true.
    #[test]
    fn basic_relationships() {
        // Step 1. Create a Qeds data structure.
        let mut qeds: Qeds<(),()> = Qeds::new();
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
        assert_eq!(e.onext().edge() as *const Edge<()>, e.edge() as *const Edge<()>);
    }

    /// Test the basic relationships of a quad that should always be true for an
    /// unconnected edge.
    #[test]
    fn floating_relationships() {
        // Step 1. Create a Qeds data structure.
        let mut qeds: Qeds<(),()> = Qeds::new();
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
        assert_eq!(e.onext().edge() as *const Edge<()>, e.edge() as *const Edge<()>);
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
    fn segment() {
        // Step 1. Create a Qeds data structure.
        let mut qeds: Qeds<Box<Option<Point>>,()> = Qeds::new();
        let point_a = Point::new(0.0, 0.0);
        let point_b = Point::new(5.0, 0.0);
        // Step 2. Add some edges to it.
        let q1 = qeds.make_edge().target;
        println!("q1 Target: {:?}", q1);
        unsafe {
            qeds.edge_a_mut(q1).point = Box::new(Some(point_a));
            qeds.edge_a_mut(q1.sym()).point = Box::new(Some(point_b));
            let e = &qeds.quads[0].edges_a[0];
            // e with sym applied 2 times is the same edge as e
            assert_eq!(qeds.edge_a_ref(q1).edge() as *const Edge<Box<Option<Point>>>, e as *const Edge<Box<Option<Point>>>);
            assert_eq!(qeds.edge_a_ref(q1).sym().sym().edge(), e);
            // The Org of the edge is a point and it is point_a
            assert_eq!(*qeds.edge_a_ref(q1).edge().point, Some(point_a));
            // The Dest of the edge (eSymOrg) is point_b
            println!("aa: {:?}", qeds.edge_a_ref(q1).target());
            println!("ab: {:?}", qeds.edge_a_ref(q1).sym().target());
            assert_eq!(*qeds.edge_a_ref(q1).sym().edge().point, Some(point_b));
        }
    }

    #[test]
    fn rot2_eq_sym() {
        // Step 1. Create a Qeds data structure.
        let mut qeds: Qeds<Box<Option<Point>>,()> = Qeds::new();
        // Step 2. Add some edges to it.
        let q1 = qeds.make_edge().target;
        let p1 = Point::new(0.0, 0.0);
        unsafe {
            qeds.edge_a_mut(q1).point = Box::new(Some(p1));
        }
        // Step 3. Get the single edge in the data structure.
        // e with Rot applied 2 times is not the same edge as e
        assert_eq!(q1.rot().rot(), q1.sym());
    }

    #[test]
    fn d5() {
        let mut qeds = Qeds::new();
        let q1 = qeds.make_edge().target;
        let q1 = unsafe { qeds.edge_a_ref(q1) };
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
        let mut qeds: Qeds<(),()> = Qeds::new();
        // Step 2. Add the first two edges. Converting them into references as
        // it is the more useful format for testing.
        let q1 = qeds.make_edge().target;
        let q2 = qeds.make_edge().target;
        unsafe {
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
    fn simple_triangle() {
        // Step 1. Create a Qeds data structure.
        let mut qeds: Qeds<Box<Option<Point>>,()> = Qeds::new();
        // Step 2. Add some edges to it.
        let a = qeds.make_edge().target;
        let b = qeds.make_edge().target;
        println!("Edges created");

        // Step 3. Splice those edges together so that we actually have
        // something of a network.
        unsafe {
            qeds.splice(a.sym(), b);
        }
    }

    #[test]
    fn triangle_face() {
        // Step 1. Create a Qeds data structure.
        let mut qeds: Qeds<Box<Option<Point>>,()> = Qeds::new();
        // Step 2. Add some edges to it.
        let q1 = qeds.make_edge().target;
        let q2 = qeds.make_edge().target;

        // Create some point to add geometry.
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(5.0, 0.0);
        let p3 = Point::new(2.5, 5.0);

        unsafe {
            qeds.edge_a_mut(q1).point = Box::new(Some(p1));
            qeds.edge_a_mut(q1.sym()).point = Box::new(Some(p2));
            qeds.edge_a_mut(q2).point = Box::new(Some(p2));
            qeds.edge_a_mut(q2.sym()).point = Box::new(Some(p3));

            // Step 3. Splice those edges together so that we actually have
            // something of a network.
            qeds.splice(q1.sym(), q2);
            assert_eq!(qeds.edge_a_ref(q1).l_next(), qeds.edge_a_ref(q2));
            assert_eq!(qeds.edge_a_ref(q2).onext(), qeds.edge_a_ref(q1).sym());
            assert_eq!(
                qeds.edge_a_ref(q1).rot().sym().onext().rot(),
                qeds.edge_a_ref(q2)
            );
            let q3 = qeds.connect(q2, q1).target;
            assert_eq!(qeds.edge_a_ref(q2).sym().onext(), qeds.edge_a_ref(q3));

            let edge = qeds.edge_a_ref(EdgeTarget::new(0, 0, 0));
            assert_eq!(edge.l_next(), qeds.edge_a_ref(q2));
            assert_eq!(edge.l_next().l_next(), edge.onext().sym());
            assert_eq!(edge.l_next().l_next().l_next(), edge);

            assert_eq!(*edge.edge().point, Some(p1));
            assert_eq!(*edge.sym().edge().point, Some(p2));
            assert_eq!(*edge.l_next().edge().point, Some(p2));
            assert_eq!(*edge.l_next().sym().edge().point, Some(p3));
            assert_eq!(*edge.l_next().l_next().edge().point, Some(p3));
            assert_eq!(*edge.l_next().l_next().sym().edge().point, Some(p1));
            // Get the face from it.
            let face = edge.l_face();
            for (edge, (p1, p2)) in face
                .edges
                .iter()
                .zip(vec![(p2, p3), (p3, p1), (p1, p2)].into_iter())
            {
                println!(
                    "Edge[]: {:?} -> {:?}",
                    edge.edge().point,
                    edge.sym().edge().point
                );
                assert_eq!(*edge.edge().point, Some(p1));
                assert_eq!(*edge.sym().edge().point, Some(p2));
            }
            // It should have 3 edges.
            assert_eq!(face.edges.len(), 3);
            assert_eq!(face.vertices().count(), 3);

            assert_eq!(
                face.centroid(),
                Some(Point {
                    x: 2.5,
                    y: 1.6666666666666665
                })
            );
        }
    }

    // #[ignore]
    #[test]
    fn quad() {
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(5.0, 0.0);
        let p3 = Point::new(2.5, 5.0);
        let p4 = Point::new(5.0, 5.0);
        // All the geometry data is set up here. This makes the whole thing
        // static of course.
        //
        // Step 1. Create a Qeds data structure.
        let mut qeds: Qeds<Box<Option<Point>>,()> = Qeds::new();
        // Step 2. Add the first two edges.
        let q1 = qeds.make_edge().target;
        let q2 = qeds.make_edge().target;
        let q4 = qeds.make_edge().target;
        // Step 3. Assign the appropriate geometry to the edges.
        unsafe {
            qeds.edge_a_mut(q1).point = Box::new(Some(p1));
            qeds.edge_a_mut(q1.sym()).point = Box::new(Some(p2));
            qeds.edge_a_mut(q2).point = Box::new(Some(p2));
            qeds.edge_a_mut(q2.sym()).point = Box::new(Some(p3));

            qeds.edge_a_mut(q4).point = Box::new(Some(p4));
            qeds.edge_a_mut(q4.sym()).point = Box::new(Some(p3));
        }

        // Step 4. Splice those edges together so that we actually have
        // something of a network. This adds the third edge.
        unsafe {
            qeds.splice(q1.sym(), q2);
            assert_eq!(qeds.edge_a_ref(q1).sym().onext().target, q2);
            let e = qeds.connect(q2, q1).target();
            // At this point the triangle has been created.
            {
                assert_eq!(qeds.edge_a_ref(q1).l_next().target, q2);
                assert_eq!(*qeds.edge_a(q1).point, Some(p1));
                assert_eq!(*qeds.edge_a(q1.sym()).point, Some(p2));
                assert_eq!(*qeds.edge_a_ref(q1).l_next().edge().point, Some(p2));
                assert_eq!(*qeds.edge_a_ref(q1).l_next().sym().edge().point, Some(p3));
                assert_eq!(*qeds.edge_a_ref(q1).l_next().l_next().edge().point, Some(p3));
                assert_eq!(
                    *qeds.edge_a_ref(EdgeTarget::new(3, 0, 0)).edge().point,
                    Some(p3)
                );
                assert_eq!(
                    *qeds.edge_a_ref(EdgeTarget::new(3, 0, 0)).sym().edge().point,
                    Some(p1)
                );
                assert_eq!(
                    *qeds.edge_a_ref(q1).l_next().l_next().sym().edge().point,
                    Some(p1)
                );
            }

            assert_eq!(qeds.edge_a_ref(e).onext().target, q2.sym());

            // Now we want to splice on the fourth edge.
            qeds.splice(q4.sym(), q2.sym());
            // 1. dSymOnext == c
            assert_eq!(qeds.edge_a_ref(q4).sym().onext().target(), e);
            // 2. bSymOnext == dSym
            assert_eq!(qeds.edge_a_ref(q2).sym().onext().target(), q4.sym());
            // 3. dSymRotOnext == bRot
            assert_eq!(qeds.edge_a_ref(q4).sym().rot().onext().target(), q2.rot());
            // 4. cRotOnext == dRot
            assert_eq!(qeds.edge_a_ref(e).rot().onext().target(), q4.rot());

            // Now we add in the fifth edge, closing the quad.
            qeds.connect(q2.sym(), q4);

            // Check the first triangle face
            assert_eq!(
                qeds.edge_a_ref(q1).rot().sym().onext().target(),
                q2.rot().sym()
            );
        }
        unsafe {
            // Get the first edge.
            let edge = qeds.edge_a_ref(EdgeTarget::new(0, 0, 0));

            println!(
                "Edge1[{:?}]: {:?} -> {:?}",
                edge.edge() as *const Edge<Box<Option<Point>>>,
                edge.edge().point,
                edge.sym().edge().point
            );
            println!(
                "Edge2[{:?}]: {:?} -> {:?}",
                edge.l_next().edge() as *const Edge<Box<Option<Point>>>,
                edge.l_next().edge().point,
                edge.l_next().sym().edge().point
            );
            println!(
                "Edge3[{:?}]: {:?} -> {:?}",
                edge.l_next().l_next().edge() as *const Edge<Box<Option<Point>>>,
                edge.l_next().l_next().edge().point,
                edge.l_next().l_next().sym().edge().point
            );
            println!(
                "Edge4[{:?}]: {:?} -> {:?}",
                edge.l_next().l_next().l_next().edge() as *const Edge<Box<Option<Point>>>,
                edge.l_next().l_next().l_next().edge().point,
                edge.l_next().l_next().l_next().sym().edge().point
            );

            assert_eq!(edge.l_next().target(), q2);
            assert_eq!(*edge.edge().point, Some(p1));
            assert_eq!(*edge.sym().edge().point, Some(p2));
            assert_eq!(*edge.l_next().edge().point, Some(p2));
            assert_eq!(*edge.l_next().sym().edge().point, Some(p3));
            assert_eq!(*edge.l_next().l_next().edge().point, Some(p3));
            assert_eq!(*edge.l_next().l_next().sym().edge().point, Some(p1));

            // Get the face from it.
            let face = edge.l_face();
            // It should have 3 edges.
            assert_eq!(face.edges.len(), 3);
            let edge1 = qeds.base_edges().next().unwrap();
            assert_eq!(edge1.l_face().edges.len(), 3);

            let mut midpoints = Vec::new();
            println!("looping through face");
            for (edge, (p1, p2)) in face
                .edges
                .iter()
                .zip(vec![(p2, p3), (p3, p1), (p1, p2)].into_iter())
            {
                println!(
                    "Edge[]: {:?} -> {:?}",
                    edge.edge().point,
                    edge.sym().edge().point
                );
                assert_eq!(*edge.edge().point, Some(p1));
                assert_eq!(*edge.sym().edge().point, Some(p2));
                midpoints.push(edge.midpoint());
            }
            let mut centre = Point::new(0.0, 0.0);
            let n = midpoints.len();
            println!("midpoints: {:?}", midpoints);
            for p in midpoints.into_iter() {
                centre.x += p.unwrap().x;
                centre.y += p.unwrap().y;
            }
            centre.x = centre.x / (n as f64);
            centre.y = centre.y / (n as f64);
            assert_eq!(centre, Point::new(2.5, 5.0 / 3.0))
        }
    }
}
