use core::ops::Add;
use core::ops::Sub;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::ops::Mul;

/// This data structure is a single instance of a quad-edge data structure.
#[derive(Clone, Debug)]
pub struct Qeds {
    /// The vector of quads which we index into.
    pub quads: Vec<Quad>,
    /// The priority queue of "dead" quads in the vector that we can freely
    /// replace.
    pub spaces: BinaryHeap<Reverse<usize>>,
}

impl Qeds {
    /// Create the simplest [`Qeds`] with a single [`Edge`] and a single
    /// [`Face`]. Covers the whole sphere.
    pub fn new() -> Self {
        Self {
            quads: Vec::new(),
            spaces: BinaryHeap::new(),
        }
    }

    fn take_next_index(&mut self) -> usize {
        // If there are any empty spaces in the list, use the smallest of those.
        // Otherwise it's simply the length of the vector.
        self.spaces.pop().map(|i| i.0).unwrap_or(self.quads.len())
    }

    fn insert(&mut self, index: usize, quad: Quad) {
        if let Some(elem) = self.quads.get_mut(index) {
            *elem = quad;
        } else {
            self.quads.push(quad);
        }
    }

    /// Create an edge in the [`Qeds`].
    pub fn make_edge(&mut self) -> EdgeRef {
        let this_index: usize = self.take_next_index();
        let quad = Quad {
            edges: [
                // The base edge e.
                Edge {
                    next: EdgeTarget::new(this_index, 0, 0),
                    point: Box::new(None),
                },
                // eRot
                Edge {
                    next: EdgeTarget::new(this_index, 3, 0),
                    point: Box::new(None),
                },
                // eSym
                Edge {
                    next: EdgeTarget::new(this_index, 2, 0),
                    point: Box::new(None),
                },
                // eSymRot
                Edge {
                    next: EdgeTarget::new(this_index, 1, 0),
                    point: Box::new(None),
                },
            ],
        };
        self.insert(this_index, quad);
        let t: &Qeds = self;
        EdgeRef {
            qeds: t,
            target: EdgeTarget::new(this_index, 0, 0),
        }
    }

    unsafe fn edge(&self, target: EdgeTarget) -> &Edge {
        &self.quads[target.e].edges[target.r as usize]
    }

    unsafe fn edge_ref(&self, target: EdgeTarget) -> EdgeRef {
        EdgeRef {
            qeds: &self,
            target,
        }
    }

    fn edge_mut(&mut self, target: EdgeTarget) -> &mut Edge {
        &mut self.quads[target.e].edges[target.r as usize]
    }

    pub unsafe fn splice(&mut self, edge_a: EdgeTarget, edge_b: EdgeTarget) {
        let alpha = self.edge_ref(edge_a).onext().rot().target;
        let beta = self.edge_ref(edge_b).onext().rot().target;

        // We want to swap aOnext with bOnext and αOnext with βONext
        let ta = self.edge_ref(edge_a).onext().target;
        self.edge_mut(edge_a).next = self.edge_ref(edge_b).onext().target;
        self.edge_mut(edge_b).next = ta;

        let ta = self.edge_ref(alpha).onext().target;
        self.edge_mut(alpha).next = self.edge_ref(beta).onext().target;
        self.edge_mut(beta).next = ta;
    }

    /// Connect the Org of a with the Dest of b by creating a new edge.
    pub fn connect(&mut self, edge_a: EdgeTarget, edge_b: EdgeTarget) -> EdgeRef {
        // First, make the new edge.
        let q_target = self.make_edge().target;
        unsafe {
            // Set the Org of e to the Dest of a
            self.edge_mut(q_target).point = self.edge(edge_a.sym()).point.clone();
            // Set the Dest of e to the Org of b
            self.edge_mut(q_target.sym()).point = self.edge(edge_b).point.clone();

            self.splice(q_target, self.edge_ref(edge_a).l_next().target);
            self.splice(q_target.sym(), edge_b);
            self.edge_ref(q_target)
        }
    }

    // pub unsafe fn delete(&mut self, e: &mut Edge) {
    //     // TODO: we don't actually free the memory here.
    //     self.splice(e, e.oprev_mut());
    //     self.splice(e.sym_mut(), e.sym_mut().oprev_mut());
    // }

    pub unsafe fn join(&mut self, edge_a: EdgeRef, edge_b: EdgeRef) {
        self.splice(edge_a.target, edge_b.target);
        // TODO: get rid of this clone.
        self.edge_mut(edge_a.target).point = self.edge(edge_b.target).point.clone();
    }
}

#[derive(Clone, Debug, PartialEq, PartialOrd)]
pub struct Quad {
    pub edges: [Edge; 4],
}

impl Quad {}

#[derive(Copy, Clone, Debug)]
pub struct EdgeRef<'a> {
    pub qeds: &'a Qeds,
    pub target: EdgeTarget,
}

impl EdgeRef<'_> {
    fn l_face(&self) -> Face {
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
}

impl AsEdgeRef for EdgeRef<'_> {
    fn qeds(&self) -> &Qeds {
        self.qeds
    }
    fn target(&self) -> EdgeTarget {
        self.target
    }
    fn target_mut(&mut self) -> &mut EdgeTarget {
        &mut self.target
    }
}

impl PartialEq for EdgeRef<'_> {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.qeds, other.qeds) && self.target == other.target
    }
}

impl Eq for EdgeRef<'_> {}

#[derive(Debug)]
pub struct EdgeRefMut<'a> {
    pub qeds: &'a mut Qeds,
    pub target: EdgeTarget,
}

// impl AsEdgeRef for EdgeRefMut<'_> {
//     fn qeds(&self) -> &Qeds {
//         self.qeds
//     }
//     fn target(&self) -> EdgeTarget {
//         self.target
//     }
//     fn target_mut(&mut self) -> &mut EdgeTarget {
//         &mut self.target
//     }
// }

impl PartialEq for EdgeRefMut<'_> {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.qeds, other.qeds) && self.target == other.target
    }
}

impl Eq for EdgeRefMut<'_> {}

pub trait AsEdgeRef: Sized + Copy {
    fn qeds(&self) -> &Qeds;
    fn target(&self) -> EdgeTarget;
    fn target_mut(&mut self) -> &mut EdgeTarget;
    fn edge(&self) -> &Edge {
        // We know we can use this unsafe because lifetime gurantee that the
        // Qeds data structure has not been modified since this EdgeRef was
        // created.
        unsafe { self.qeds().edge(self.target()) }
    }

    #[inline(always)]
    fn offset_r(&self, offset: u8) -> Self {
        let mut q = self.clone();
        *q.target_mut() = q.target().offset_r(offset);
        q
    }

    #[inline(always)]
    fn onext(&self) -> Self {
        let edge = self.edge();
        let mut q = self.clone();
        *q.target_mut() = edge.next;
        q
    }

    /// Rot
    #[inline(always)]
    fn rot(&self) -> Self {
        self.offset_r(1)
    }

    #[inline(always)]
    fn sym(&self) -> Self {
        self.offset_r(2)
    }

    #[inline(always)]
    fn sym_rot(&self) -> Self {
        self.offset_r(3)
    }

    #[inline(always)]
    fn oprev(&self) -> Self {
        self.rot().onext().rot()
    }

    #[inline(always)]
    fn r_next(&self) -> Self {
        self.rot().onext().rot()
    }

    #[inline(always)]
    fn l_next(&self) -> Self {
        self.rot().rot().rot().onext().rot()
    }

    fn midpoint(&self) -> Option<Point> {
        let this_point = (*self.edge().point)?;
        let next_point = (*self.sym().edge().point)?;
        Some(this_point.midpoint(next_point))
    }
}

// r and f can fit in one byte. If we limit the size of qeds we could fit it in
// the edge index.

impl<'a> EdgeRef<'a> {}

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
pub struct Edge {
    // / The index of the face to the left of this edge.
    // pub face_index: usize,
    // pub index: usize,
    pub next: EdgeTarget,
    pub point: Box<Option<Point>>,
}

impl Edge {
    // /// Rot
    // #[inline(always)]
    // pub fn rot(&self) -> &Edge {
    //     self.offset(1)
    // }

    // #[inline(always)]
    // pub fn rot_mut(&mut self) -> &mut Edge {
    //     self.offset_mut(1)
    // }

    // #[inline(always)]
    // pub fn sym(&self) -> &Edge {
    //     self.offset(2)
    // }

    // #[inline(always)]
    // pub fn sym_mut(&mut self) -> &mut Edge {
    //     self.offset_mut(2)
    // }

    // #[inline(always)]
    // pub fn sym_rot(&self) -> &Edge {
    //     self.offset(3)
    // }

    // #[inline(always)]
    // fn offset(&self, offset: isize) -> &Edge {
    //     let offset1 = (self as *const Edge).align_offset(std::mem::align_of::<Quad>());
    //     // i is the index in Quad, that is, 0, 1, 2, or 3 (effectively r)
    //     let i: isize = (4 - offset1 as isize) % 4;
    //     // d is the offset from the current edge to the edge we want.
    //     let d = (i + offset) % 4 - i;
    //     let ptr = (self as *const Edge).wrapping_offset(d);
    //     unsafe { &*ptr }
    // }

    // #[inline(always)]
    // fn offset_mut(&mut self, offset: isize) -> &mut Edge {
    //     let ptr: *const Edge = self.offset(offset);
    //     unsafe { &mut *(ptr as *mut Edge) }
    // }

    // #[inline(always)]
    // pub fn onext(&self) -> &Edge {
    //     unsafe { &*self.next }
    // }

    // #[inline(always)]
    // pub fn oprev(&self) -> &Edge {
    //     self.rot().onext().rot()
    // }

    // #[inline(always)]
    // pub fn oprev_mut(&mut self) -> &mut Edge {
    //     self.rot_mut().onext_mut().rot_mut()
    // }

    // #[inline(always)]
    // pub fn onext_mut(&mut self) -> &mut Edge {
    //     unsafe {
    //         let const_next = self.next as *mut Edge;
    //         &mut *const_next
    //     }
    // }

    // #[inline(always)]
    // pub fn r_next(&self) -> &Edge {
    //     self.rot().onext().rot()
    // }

    // #[inline(always)]
    // pub fn l_next(&self) -> &Edge {
    //     self.rot().rot().rot().onext().rot()
    // }

    // #[inline(always)]
    // pub fn l_next_mut(&mut self) -> &mut Edge {
    //     self.rot_mut().rot_mut().rot_mut().onext_mut().rot_mut()
    // }

    // pub fn l_face(&self) -> Face {
    //     // Save the next edge of this face so that we know when to end the loop.
    //     let first_next_edge = self.l_next();
    //     let mut edges = Vec::new();
    //     // Loop around collecting the edges.
    //     let mut current_next_edge = self.l_next();
    //     loop {
    //         edges.push(current_next_edge);
    //         current_next_edge = current_next_edge.l_next();
    //         // TODO: use
    //         if std::ptr::eq(current_next_edge, first_next_edge) {
    //             break;
    //         }
    //     }
    //     Face { edges }
    // }
}

#[derive(Clone, Debug)]
pub struct Face<'a> {
    /// An index to any of the CCW oriented edges surrounding this face.
    pub edges: Vec<EdgeRef<'a>>,
}

impl Face<'_> {
    pub fn midpoint(&self) -> Option<Point> {
        let mut midpoints = Vec::new();
        for edge in self.edges.iter() {
            midpoints.push(edge.midpoint());
        }
        let mut centre = Point::new(0.0, 0.0);
        let n = midpoints.len();
        for p in midpoints.into_iter() {
            centre.x += p?.x;
            centre.y += p?.y;
        }
        centre.x = centre.x / (n as f64);
        centre.y = centre.y / (n as f64);
        Some(centre)
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
        let t = EdgeTarget::new(0,0,0);
        assert_eq!(t.offset_r(2), EdgeTarget::new(0,2,0));
    }

    // First we test the Edge Function Properties as described by G&S 2.3
    #[test]
    fn e1() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        // Step 2. Add some edges to it.
        let e = qeds.make_edge();
        // e with Rot applied 4 times is the same edge as e
        assert_eq!(e.rot().rot().rot().rot(), e);
    }

    /// Test the basic relationships of a quad that should always be true.
    #[test]
    fn basic_relationships() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        // Step 2. Add some edges to it.
        let e = qeds.make_edge();
        // Step 3. Get the single edge in the data structure.
        // e with Rot applied 2 times is not the same edge as e
        assert_eq!(e.rot().rot().edge() as *const Edge, e.sym().edge() as *const Edge);
        assert_eq!(e.rot().sym().edge() as *const Edge, e.sym().rot().edge() as *const Edge);
        assert_eq!(e.onext().edge() as *const Edge, e.edge() as *const Edge);
    }

    /// Test the basic relationships of a quad that should always be true for an
    /// unconnected edge.
    #[test]
    fn floating_relationships() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        // Step 2. Add some edges to it.
        let e = qeds.make_edge();
        // Step 3. Get the single edge in the data structure.
        // e with Rot applied 2 times is not the same edge as e
        assert_eq!(e.rot().rot().edge() as *const Edge, e.sym().edge() as *const Edge);
        assert_eq!(e.rot().sym().edge() as *const Edge, e.sym().rot().edge() as *const Edge);
        assert_eq!(e.onext().edge() as *const Edge, e.edge() as *const Edge);
        assert_eq!(e.sym().onext().edge() as *const Edge, e.sym().edge() as *const Edge);
        assert_eq!(e.rot().onext().edge() as *const Edge, e.rot().sym().edge() as *const Edge);
    }


    #[test]
    fn segment() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        let point_a = Point::new(0.0, 0.0);
        let point_b = Point::new(5.0, 0.0);
        // Step 2. Add some edges to it.
        let q1 = qeds.make_edge().target;
        println!("q1 Target: {:?}", q1);
        qeds.edge_mut(q1).point = Box::new(Some(point_a));
        qeds.edge_mut(q1.sym()).point = Box::new(Some(point_b));
        let e = &qeds.quads[0].edges[0];
        // e with sym applied 2 times is the same edge as e
        unsafe {
            assert_eq!(qeds.edge_ref(q1).edge() as *const Edge, e as *const Edge);
            assert_eq!(qeds.edge_ref(q1).sym().sym().edge(), e);
            // The Org of the edge is a point and it is point_a
            assert_eq!(*qeds.edge_ref(q1).edge().point, Some(point_a));
            // The Dest of the edge (eSymOrg) is point_b
            println!("aa: {:?}", qeds.edge_ref(q1).target());
            println!("ab: {:?}", qeds.edge_ref(q1).sym().target());
            assert_eq!(*qeds.edge_ref(q1).sym().edge().point, Some(point_b));
        }
    }

    #[test]
    fn rot2_eq_sym() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        // Step 2. Add some edges to it.
        let q1 = qeds.make_edge().target;
        let p1 = Point::new(0.0, 0.0);
        qeds.edge_mut(q1).point = Box::new(Some(p1));
        // Step 3. Get the single edge in the data structure.
        // e with Rot applied 2 times is not the same edge as e
        assert_eq!(q1.rot().rot(), q1.sym());
    }

    #[test]
    fn d5() {
        let mut qeds = Qeds::new();
        let q1 = qeds.make_edge().target;
        let q1 = unsafe { qeds.edge_ref(q1) };
        assert_eq!(q1.rot().edge(), &qeds.quads[0].edges[1]);
        assert_eq!(q1.sym().edge(), &qeds.quads[0].edges[2]);
        assert_eq!(q1.sym().rot().edge(), &qeds.quads[0].edges[3]);
        assert_eq!(
            q1.rot().sym().edge() as *const Edge,
            &qeds.quads[0].edges[3] as *const Edge
        );
        assert_eq!(
            q1.rot().onext().edge() as *const Edge,
            q1.rot().sym().edge() as *const Edge
        );
        assert_eq!(
            q1.rot().sym().onext().edge() as *const Edge,
            q1.rot().edge() as *const Edge
        );
        assert_eq!(
            q1.onext().edge() as *const Edge,
            &qeds.quads[0].edges[0] as *const Edge
        );
    }

    #[test]
    fn splice_test() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        // Step 2. Add the first two edges. Converting them into references as
        // it is the more useful format for testing.
        let q1 = qeds.make_edge().target;
        let q2 = qeds.make_edge().target;
        unsafe {
            qeds.splice(q1.sym(), q2);
            let q1 = qeds.edge_ref(q1);
            let q2 = qeds.edge_ref(q2);
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
    fn triangle_face() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        // Step 2. Add some edges to it.
        let q1 = qeds.make_edge().target;
        let q2 = qeds.make_edge().target;

        // let q1 = unsafe {qeds.edge_ref(q1)};
        // let q2 = unsafe {qeds.edge_ref(q2)};
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(5.0, 0.0);
        let p3 = Point::new(2.5, 5.0);

        qeds.edge_mut(q1).point = Box::new(Some(p1));
        qeds.edge_mut(q1.sym()).point = Box::new(Some(p2));
        qeds.edge_mut(q2).point = Box::new(Some(p2));
        qeds.edge_mut(q2.sym()).point = Box::new(Some(p3));

        // Step 3. Splice those edges together so that we actually have
        // something of a network.
        unsafe {
            qeds.splice(q1.sym(), q2);
            assert_eq!(qeds.edge_ref(q1).l_next(), qeds.edge_ref(q2));
            assert_eq!(qeds.edge_ref(q2).onext(), qeds.edge_ref(q1).sym());
            assert_eq!(
                qeds.edge_ref(q1).rot().sym().onext().rot(),
                qeds.edge_ref(q2)
            );
            let q3 = qeds.connect(q2, q1).target;
            assert_eq!(qeds.edge_ref(q2).sym().onext(), qeds.edge_ref(q3));

            // Get the first edge.
            let quad = qeds.quads.iter().next().unwrap();
            // let quad: &Quad = unsafe { &**quad };
            let edge = qeds.edge_ref(EdgeTarget::new(0, 0, 0));
            assert_eq!(edge.l_next(), qeds.edge_ref(q2));
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
        let mut qeds = Qeds::new();
        // Step 2. Add the first two edges.
        let q1 = qeds.make_edge().target;
        let q2 = qeds.make_edge().target;
        let q4 = qeds.make_edge().target;
        // Step 3. Assign the appropriate geometry to the edges.
        unsafe {
            qeds.edge_mut(q1).point = Box::new(Some(p1));
            qeds.edge_mut(q1.sym()).point = Box::new(Some(p2));
            qeds.edge_mut(q2).point = Box::new(Some(p2));
            qeds.edge_mut(q2.sym()).point = Box::new(Some(p3));

            qeds.edge_mut(q4).point = Box::new(Some(p4));
            qeds.edge_mut(q4.sym()).point = Box::new(Some(p3));
        }

        // Step 4. Splice those edges together so that we actually have
        // something of a network. This adds the third edge.
        unsafe {
            qeds.splice(q1.sym(), q2);
            assert_eq!(qeds.edge_ref(q1).sym().onext().target, q2);
            let e = qeds.connect(q2, q1).target();
            // At this point the triangle has been created.
            {
                assert_eq!(qeds.edge_ref(q1).l_next().target, q2);
                assert_eq!(*qeds.edge(q1).point, Some(p1));
                assert_eq!(*qeds.edge(q1.sym()).point, Some(p2));
                assert_eq!(*qeds.edge_ref(q1).l_next().edge().point, Some(p2));
                assert_eq!(*qeds.edge_ref(q1).l_next().sym().edge().point, Some(p3));
                assert_eq!(*qeds.edge_ref(q1).l_next().l_next().edge().point, Some(p3));
                assert_eq!(*qeds.edge_ref(EdgeTarget::new(3,0,0)).edge().point, Some(p3));
                assert_eq!(*qeds.edge_ref(EdgeTarget::new(3,0,0)).sym().edge().point, Some(p1));
                assert_eq!(
                    *qeds.edge_ref(q1).l_next().l_next().sym().edge().point,
                    Some(p1)
                );
            }

            assert_eq!(qeds.edge_ref(e).onext().target, q2.sym());

            // Now we want to splice on the fourth edge.
            qeds.splice(q4.sym(), q2.sym());
            // 1. dSymOnext == c
            assert_eq!(qeds.edge_ref(q4).sym().onext().target(), e);
            // 2. bSymOnext == dSym
            assert_eq!(qeds.edge_ref(q2).sym().onext().target(), q4.sym());
            // 3. dSymRotOnext == bRot
            assert_eq!(qeds.edge_ref(q4).sym().rot().onext().target(), q2.rot());
            // 4. cRotOnext == dRot
            assert_eq!(qeds.edge_ref(e).rot().onext().target(), q4.rot());

            // Now we add in the fifth edge, closing the quad.
            qeds.connect(q2.sym(), q4);

            // Check the first triangle face
            assert_eq!(
                qeds.edge_ref(q1).rot().sym().onext().target(),
                q2.rot().sym()
            );
        }
        unsafe {
            // Get the first edge.
            let edge = qeds.edge_ref(EdgeTarget::new(0, 0, 0));

            println!(
                "Edge1[{:?}]: {:?} -> {:?}",
                edge.edge() as *const Edge,
                edge.edge().point,
                edge.sym().edge().point
            );
            println!(
                "Edge2[{:?}]: {:?} -> {:?}",
                edge.l_next().edge() as *const Edge,
                edge.l_next().edge().point,
                edge.l_next().sym().edge().point
            );
            println!(
                "Edge3[{:?}]: {:?} -> {:?}",
                edge.l_next().l_next().edge() as *const Edge,
                edge.l_next().l_next().edge().point,
                edge.l_next().l_next().sym().edge().point
            );
            println!(
                "Edge4[{:?}]: {:?} -> {:?}",
                edge.l_next().l_next().l_next().edge() as *const Edge,
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
