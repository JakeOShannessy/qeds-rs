#![feature(new_uninit)]
#![feature(ptr_wrapping_offset_from)]
use core::mem::MaybeUninit;
use core::ops::Add;
use core::ops::Sub;
use std::ops::Mul;

/// This data structure is a single instance of a quad-edge data structure.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Qeds {
    pub quads: Vec<*mut Quad>,
    // pub faces: Vec<Face>,
}

impl Qeds {
    /// Create the simplest [`Qeds`] with a single [`Edge`] and a single
    /// [`Face`]. Covers the whole sphere.
    pub fn new() -> Self {
        let quads = Vec::new();
        Self { quads }
    }

    /// Create an edge in the [`Qeds`].
    pub fn make_edge(&mut self) -> *mut Quad {
        // let quad_layout = std::alloc::Layout::from_size_align(std::mem::size_of::<Quad>(), 128).unwrap();
        // let quad_layout_standard = std::alloc::Layout::from_size_align(std::mem::size_of::<Quad>(), 64).unwrap();
        // println!("Quad Layout: {:?}", quad_layout);
        // println!("Quad Layout Std: {:?}",quad_layout_standard);
        // let quad = unsafe {std::alloc::alloc_zeroed(quad_layout)};
        // let mut quad: Box<Quad> = unsafe {Box::from_raw(quad as *mut Quad)};

        let quad: Box<MaybeUninit<Quad>> = Box::new_zeroed();
        let mut quad = unsafe { quad.assume_init() };

        // The base edge e.
        quad.edges[0].next = &quad.edges[0];
        quad.edges[0].point = Box::new(None);
        // eRot
        quad.edges[1].next = &quad.edges[3];
        quad.edges[1].point = Box::new(None);
        // eSym
        quad.edges[2].next = &quad.edges[2];
        quad.edges[2].point = Box::new(None);
        // eSymRot
        quad.edges[3].next = &quad.edges[1];
        quad.edges[3].point = Box::new(None);
        let p = Box::into_raw(quad);
        self.quads.push(p);
        p
    }

    pub unsafe fn splice(&mut self, edge_a: &mut Edge, edge_b: &mut Edge) {
        let alpha_p = edge_a.onext().rot() as *const Edge;
        let beta_p = edge_b.onext().rot() as *const Edge;
        let alpha = alpha_p as *mut Edge;
        let beta = beta_p as *mut Edge;

        // We want to swap aOnext with bOnext and αOnext with βONext
        let ta = edge_a.onext() as *const Edge;
        edge_a.next = edge_b.next;
        edge_b.next = ta;

        let ta = (*alpha).onext() as *const Edge;
        (*alpha).next = (*beta).next;
        (*beta).next = ta;
    }

    /// Connect the Org of a with the Dest of b by creating a new edge.
    pub unsafe fn connect(&mut self, edge_a: &mut Edge, edge_b: &mut Edge) -> *mut Quad {
        // First, make the new edge.
        let q = self.make_edge();
        // Set the Org of e to the Dest of a
        (*q).edges[0].point = edge_a.sym().point.clone();
        // Set the Dest of e to the Org of b
        ((*q).edges[2]).point = edge_b.point.clone();
        self.splice(&mut (*q).edges[0],edge_a.l_next_mut());
        self.splice(&mut (*q).edges[2], edge_b);
        q
    }

    // pub unsafe fn delete(&mut self, e: &mut Edge) {
    //     // TODO: we don't actually free the memory here.
    //     self.splice(e, e.oprev_mut());
    //     self.splice(e.sym_mut(), e.sym_mut().oprev_mut());
    // }

    pub unsafe fn join(&mut self, edge_a: &mut Edge, edge_b: &mut Edge) {
        self.splice(edge_a, edge_b);
        let edge_a_mut = edge_a as *mut Edge;
        // TODO: get rid of this clone.
        (*edge_a_mut).point = edge_b.point.clone();
    }
}

#[derive(Clone, Debug, PartialEq, PartialOrd)]
// TODO: adjust this based on pointer width.
#[repr(align(64))]
pub struct Quad {
    pub edges: [Edge; 4],
}

impl Quad {}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct QuadRef<'a> {
    pub qeds: &'a Qeds,
    pub target: QuadTarget,
}

// r and f can fit in one byte. If we limit the size of qeds we could fit it in
// the edge index.

impl<'a> QuadRef<'a> {}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct QuadTarget {
    // The index of the Quad in the Qeds.
    pub e: usize,
    // Can only be 0, 1, 2, or 3.
    pub r: u8,
    // Can only be 0 or 1.
    pub f: u8,
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
    pub next: *const Edge,
    pub point: Box<Option<Point>>,
}

impl Edge {
    /// Rot
    #[inline(always)]
    pub fn rot(&self) -> &Edge {
        self.offset(1)
    }

    #[inline(always)]
    pub fn rot_mut(&mut self) -> &mut Edge {
        self.offset_mut(1)
    }

    #[inline(always)]
    pub fn sym(&self) -> &Edge {
        self.offset(2)
    }

    #[inline(always)]
    pub fn sym_mut(&mut self) -> &mut Edge {
        self.offset_mut(2)
    }

    #[inline(always)]
    pub fn sym_rot(&self) -> &Edge {
        self.offset(3)
    }

    #[inline(always)]
    fn offset(&self, offset: isize) -> &Edge {
        let offset1 = (self as *const Edge).align_offset(std::mem::align_of::<Quad>());
        // i is the index in Quad, that is, 0, 1, 2, or 3 (effectively r)
        let i: isize = (4 - offset1 as isize) % 4;
        // d is the offset from the current edge to the edge we want.
        let d = (i + offset) % 4 - i;
        let ptr = (self as *const Edge).wrapping_offset(d);
        unsafe { &*ptr }
    }

    #[inline(always)]
    fn offset_mut(&mut self, offset: isize) -> &mut Edge {
        let ptr: *const Edge = self.offset(offset);
        unsafe { &mut *(ptr as *mut Edge) }
    }

    #[inline(always)]
    pub fn onext(&self) -> &Edge {
        unsafe { &*self.next }
    }

    #[inline(always)]
    pub fn oprev(&self) -> &Edge {
        self.rot().onext().rot()
    }

    #[inline(always)]
    pub fn oprev_mut(&mut self) -> &mut Edge {
        self.rot_mut().onext_mut().rot_mut()
    }

    #[inline(always)]
    pub fn onext_mut(&mut self) -> &mut Edge {
        unsafe {
            let const_next = self.next as *mut Edge;
            &mut *const_next
        }
    }

    #[inline(always)]
    pub fn r_next(&self) -> &Edge {
        self.rot().onext().rot()
    }

    #[inline(always)]
    pub fn l_next(&self) -> &Edge {
        self.rot().rot().rot().onext().rot()
    }

    #[inline(always)]
    pub fn l_next_mut(&mut self) -> &mut Edge {
        self.rot_mut().rot_mut().rot_mut().onext_mut().rot_mut()
    }

    pub fn l_face(&self) -> Face {
        // Save the next edge of this face so that we know when to end the loop.
        let first_next_edge = self.l_next();
        let mut edges = Vec::new();
        // Loop around collecting the edges.
        let mut current_next_edge = self.l_next();
        loop {
            edges.push(current_next_edge);
            current_next_edge = current_next_edge.l_next();
            // TODO: use
            if std::ptr::eq(current_next_edge, first_next_edge) {
                break;
            }
        }
        Face { edges }
    }

    pub fn midpoint(&self) -> Option<Point> {
        let this_point = (*self.point)?;
        let next_point = (*self.sym().point)?;
        Some(this_point.midpoint(next_point))
    }
}

#[derive(Clone, Debug, PartialEq, PartialOrd)]
pub struct Face<'a> {
    /// An index to any of the CCW oriented edges surrounding this face.
    pub edges: Vec<&'a Edge>,
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
            x: (self.x + other.x)/2.0,
            y: (self.y + other.y)/2.0,
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
            x: self.x*other,
            y: self.y*other,
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

    // First we test the Edge Function Properties as described by G&S 2.3
    #[test]
    fn e1() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        // Step 2. Add some edges to it.
        qeds.make_edge();
        // let p1 = Point::new(0.0, 0.0);
        // Step 3. Get the single edge in the data structure.
        let e: &Edge = unsafe { &(*qeds.quads[0]).edges[0] };
        // e with Rot applied 4 times is the same edge as e
        assert_eq!(e.rot().rot().rot().rot() as *const Edge, e as *const Edge);
        assert!(std::ptr::eq(e.rot().rot().rot().rot() as *const Edge, e as *const Edge));
    }

    #[test]
    fn e3() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        // Step 2. Add some edges to it.
        qeds.make_edge();
        // Step 3. Get the single edge in the data structure.
        let e: &Edge = unsafe { &(*qeds.quads[0]).edges[0] };
        // e with Rot applied 2 times is not the same edge as e
        assert_ne!(e.rot().rot() as *const Edge, e as *const Edge);
    }

    #[test]
    fn segment() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        let point_a = Point::new(0.0, 0.0);
        let point_b = Point::new(5.0, 0.0);
        // Step 2. Add some edges to it.
        let q1 = qeds.make_edge();
        // let q2 = qeds.make_edge();
        unsafe {
            (*q1).edges[0].point = Box::new(Some(point_a));
            (*q1).edges[2].point = Box::new(Some(point_b));
        }

        // Step 3. Get the single edge in the data structure.
        let e: &Edge = unsafe { &(*qeds.quads[0]).edges[0] };
        // e with sym applied 2 times is the same edge as e
        assert_eq!(e.sym().sym() as *const Edge, e as *const Edge);
        // The Org of the edge is a point and it is point_a
        assert_eq!((*e.point), Some(point_a));
        // The Dest of the edge (eSymOrg) is point_b
        assert_eq!((*e.sym().point), Some(point_b));
    }


    #[test]
    fn rot2_eq_sym() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        // Step 2. Add some edges to it.
        let q1 = qeds.make_edge();
        let p1 = Point::new(0.0, 0.0);
        unsafe {
            (*q1).edges[0].point = Box::new(Some(p1));
        }
        // Step 3. Get the single edge in the data structure.
        let e: &Edge = unsafe { &(*qeds.quads[0]).edges[0] };
        // e with Rot applied 2 times is not the same edge as e
        assert_eq!(e.rot().rot() as *const Edge, e.sym() as *const Edge);
    }



    #[test]
    fn d5() {
        let mut qeds = Qeds::new();
        let q1 = qeds.make_edge();
        unsafe {
            assert_eq!((*q1).edges[0].rot(), &(*q1).edges[1]);
            assert_eq!((*q1).edges[0].sym(), &(*q1).edges[2]);
            assert_eq!((*q1).edges[0].sym().rot(), &(*q1).edges[3]);
            assert_eq!((*q1).edges[0].rot().sym(), &(*q1).edges[3]);
            assert_eq!((*q1).edges[0].rot().onext() as *const Edge, (*q1).edges[0].rot().sym() as *const Edge);
            assert_eq!((*q1).edges[0].rot().sym().onext() as *const Edge, (*q1).edges[0].rot() as *const Edge);
            assert_eq!((*q1).edges[1].onext() as *const Edge, &(*q1).edges[3] as *const Edge);
        }
    }

    #[test]
    fn splice_test() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        // Step 2. Add the first two edges.
        let q1 = qeds.make_edge();
        let q2 = qeds.make_edge();
        unsafe {
            println!("q1: {:?}", &(*q1).edges[0] as *const Edge);
            println!("q1Sym: {:?}", &(*q1).edges[2] as *const Edge);
            println!("q1Onext: {:?}", (*q1).edges[0].onext() as *const Edge);
            println!("q1Rot: {:?}", (*q1).edges[0].rot() as *const Edge);
            println!("q1RotSym: {:?}", (*q1).edges[0].rot().sym() as *const Edge);
            println!("q1RotOnext: {:?}", (*q1).edges[0].rot().onext() as *const Edge);
            println!("q2: {:?}", &(*q2).edges[0] as *const Edge);
            println!("q2Onext: {:?}", (*q2).edges[0].onext() as *const Edge);
            println!("q2Rot: {:?}", (*q2).edges[0].rot() as *const Edge);
            println!("q2RotSym: {:?}", (*q2).edges[0].rot().sym() as *const Edge);
            println!("q2RotOnext: {:?}", (*q2).edges[0].rot().onext() as *const Edge);
            assert_eq!((*q1).edges[1].next, &(*q1).edges[3] as *const Edge);
            qeds.splice(&mut (*q1).edges[2], &mut (*q2).edges[0]);
            assert_eq!((*q1).edges[2].onext(), &(*q2).edges[0]);
            println!("q1: {:?}", &(*q1).edges[0] as *const Edge);
            println!("q1Sym: {:?}", &(*q1).edges[2] as *const Edge);
            println!("q1Onext: {:?}", (*q1).edges[0].onext() as *const Edge);
            println!("q1Rot: {:?}", (*q1).edges[0].rot() as *const Edge);
            println!("q1RotSym: {:?}", (*q1).edges[0].rot().sym() as *const Edge);
            println!("q1RotOnext: {:?}", (*q1).edges[0].rot().onext() as *const Edge);
            println!("q2: {:?}", &(*q2).edges[0] as *const Edge);
            println!("q2Onext: {:?}", (*q2).edges[0].onext() as *const Edge);
            println!("q2Rot: {:?}", (*q2).edges[0].rot() as *const Edge);
            println!("q2RotSym: {:?}", (*q2).edges[0].rot().sym() as *const Edge);
            println!("q2RotOnext: {:?}", (*q2).edges[0].rot().onext() as *const Edge);

            // 1. aSymOnext == b
            assert_eq!((*q1).edges[0].sym().onext(), &(*q2).edges[0]);
            // 2. bOnext == aSym
            assert_eq!((*q2).edges[0].onext(),(*q1).edges[0].sym());
            // 3. aRotOnext == bRot
            assert_ne!((*q1).edges[0].rot().onext(),(*q1).edges[0].rot());
            // 4. bRotOnext == aRot
            assert_eq!((*q2).edges[0].rot().onext(),(*q1).edges[0].rot());

            // bLenxt = bSym
            assert_eq!((*q2).edges[0].l_next(),&(*q2).edges[2]);


        }
    }

    #[test]
    fn triangle_face() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        // Step 2. Add some edges to it.
        let q1 = qeds.make_edge();
        let q2 = qeds.make_edge();

        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(5.0, 0.0);
        let p3 = Point::new(2.5, 5.0);

        unsafe {
            (*q1).edges[0].point = Box::new(Some(p1));
            (*q1).edges[2].point = Box::new(Some(p2));
            (*q2).edges[0].point = Box::new(Some(p2));
            (*q2).edges[2].point = Box::new(Some(p3));
        }

        // Step 3. Splice those edges together so that we actually have
        // something of a network.
        unsafe {
            qeds.splice(&mut (*q1).edges[2], &mut (*q2).edges[0]);
            assert_eq!((*q1).edges[0].l_next(), &(*q2).edges[0]);
            assert_eq!((*q2).edges[0].onext(), &(*q1).edges[2]);
            assert_eq!((*q1).edges[0].rot().sym().onext().rot(), &(*q2).edges[0]);
            let q3 = qeds.connect(&mut (*q2).edges[0], &mut (*q1).edges[0]);
            assert_eq!((*q2).edges[0].sym().onext(), &(*q3).edges[0]);
        }

        // Get the first edge.
        let quad = qeds.quads.iter().next().unwrap();
        let quad: &Quad = unsafe { &**quad };
        let edge = &quad.edges[0];
        unsafe{
            println!("Edge1: {:?} -> {:?}", edge, edge.sym().point);
            println!("Edge2: {:?} -> {:?}", edge.l_next(), edge.l_next().sym().point);
            println!("Edge3: {:?} -> {:?}", edge.l_next().l_next(), edge.l_next().l_next().sym().point);
            println!("Edge4: {:?} -> {:?}", edge.l_next().l_next().l_next(), edge.l_next().l_next().l_next().sym().point);

            assert_eq!(edge.l_next() as *const Edge, &(*q2).edges[0] as *const Edge);
            assert_eq!(edge.l_next().l_next() as *const Edge, &(*qeds.quads[2]).edges[0] as *const Edge);
            assert_eq!(edge.l_next().l_next().l_next() as *const Edge, edge);

            assert_eq!(*edge.point, Some(p1));
            assert_eq!(*edge.sym().point, Some(p2));
            assert_eq!(*edge.l_next().point, Some(p2));
            assert_eq!(*edge.l_next().sym().point, Some(p3));
            assert_eq!(*edge.l_next().l_next().point, Some(p3));
            assert_eq!(*edge.l_next().l_next().sym().point, Some(p1));
        }
        // Get the face from it.
        let face = edge.l_face();
        for (edge, (p1,p2)) in face.edges.iter().zip(vec![(p2,p3),(p3,p1),(p1,p2)].into_iter()) {
            println!("Edge[]: {:?} -> {:?}", edge.point, edge.sym().point);
            assert_eq!(*edge.point, Some(p1));
            assert_eq!(*edge.sym().point, Some(p2));
        }
        // It should have 3 edges.
        assert_eq!(face.edges.len(), 3);
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
        let q1 = qeds.make_edge();
        let q2 = qeds.make_edge();
        let q4 = qeds.make_edge();
        // Step 3. Assign the appropriate geometry to the edges.
        unsafe {
            (*q1).edges[0].point = Box::new(Some(p1));
            (*q1).edges[2].point = Box::new(Some(p2));
            (*q2).edges[0].point = Box::new(Some(p2));
            (*q2).edges[2].point = Box::new(Some(p3));

            (*q4).edges[0].point = Box::new(Some(p4));
            (*q4).edges[2].point = Box::new(Some(p3));
        }

        // Step 4. Splice those edges together so that we actually have
        // something of a network. This adds the third edge.
        unsafe {
            qeds.splice(&mut (*q1).edges[2], &mut (*q2).edges[0]);
            assert_eq!((*q1).edges[2].onext(), &(*q2).edges[0]);
            println!("q1: {:?}", &(*q1).edges[0] as *const Edge);
            println!("q1Sym: {:?}", &(*q1).edges[2] as *const Edge);
            println!("q1Onext: {:?}", (*q1).edges[0].onext() as *const Edge);
            println!("q1Rot: {:?}", (*q1).edges[0].rot() as *const Edge);
            println!("q1RotSym: {:?}", (*q1).edges[0].rot().sym() as *const Edge);
            println!("q1RotOnext: {:?}", (*q1).edges[0].rot().onext() as *const Edge);
            println!("q2: {:?}", &(*q2).edges[0] as *const Edge);
            println!("q2Onext: {:?}", (*q2).edges[0].onext() as *const Edge);
            println!("q2Rot: {:?}", (*q2).edges[0].rot() as *const Edge);
            println!("q2RotSym: {:?}", (*q2).edges[0].rot().sym() as *const Edge);
            println!("q2RotOnext: {:?}", (*q2).edges[0].rot().onext() as *const Edge);

            qeds.connect(&mut (*q2).edges[0], &mut (*q1).edges[0]);
            // At this point the triangle has been created.

            {
                assert_eq!((*q1).edges[0].l_next() as *const Edge, &(*q2).edges[0] as *const Edge);
                assert_eq!(*(*q1).edges[0].point, Some(p1));
                assert_eq!(*(*q1).edges[0].sym().point, Some(p2));
                assert_eq!(*(*q1).edges[0].l_next().point, Some(p2));
                assert_eq!(*(*q1).edges[0].l_next().sym().point, Some(p3));
                assert_eq!(*(*q1).edges[0].l_next().l_next().point, Some(p3));
                assert_eq!(*(*q1).edges[0].l_next().l_next().sym().point, Some(p1));
            }

            // This is the third edge from the triangle.
            let e_p = (*q2).edges[2].next as *const Edge;
            let e = e_p as *mut Edge;

            assert_eq!((&*e).onext(), &(*q2).edges[2]);

            // Now we want to splice on the fourth edge.
            qeds.splice(&mut (*q4).edges[2], &mut (*q2).edges[2]);
            // 1. dSymOnext == c
            assert_eq!((*q4).edges[2].onext(), &*e);
            // 2. bSymOnext == dSym
            assert_eq!((*q2).edges[2].onext(), &(*q4).edges[2]);
            // 3. dSymRotOnext == bRot
            assert_eq!((*q4).edges[2].rot().onext(), (*q2).edges[0].rot());
            // 4. cRotOnext == dRot
            assert_eq!((&*e).rot().onext(), (*q4).edges[0].rot());

            // Now we add in the fifth edge, closing the quad.
            qeds.connect(&mut (*q2).edges[2], &mut (*q4).edges[0]);

            // Check the first triangle face
            assert_eq!((*q1).edges[0].rot().sym().onext(), (*q2).edges[0].rot().sym());

        }
        // Get the first edge.
        let quad = qeds.quads.iter().next().unwrap();
        let quad: &Quad = unsafe { &**quad };
        let edge = &quad.edges[0];

        unsafe{
            println!("Edge1[{:?}]: {:?} -> {:?}", edge as *const Edge, edge.point, edge.sym().point);
            println!("Edge2[{:?}]: {:?} -> {:?}", edge.l_next() as *const Edge, edge.l_next().point, edge.l_next().sym().point);
            println!("Edge3[{:?}]: {:?} -> {:?}", edge.l_next().l_next() as *const Edge, edge.l_next().l_next().point, edge.l_next().l_next().sym().point);
            println!("Edge4[{:?}]: {:?} -> {:?}", edge.l_next().l_next().l_next() as *const Edge, edge.l_next().l_next().l_next().point, edge.l_next().l_next().l_next().sym().point);

            assert_eq!(edge.l_next() as *const Edge, &(*q2).edges[0] as *const Edge);
            assert_eq!(*edge.point, Some(p1));
            assert_eq!(*edge.sym().point, Some(p2));
            assert_eq!(*edge.l_next().point, Some(p2));
            assert_eq!(*edge.l_next().sym().point, Some(p3));
            assert_eq!(*edge.l_next().l_next().point, Some(p3));
            assert_eq!(*edge.l_next().l_next().sym().point, Some(p1));
        }

        // Get the face from it.
        let face = edge.l_face();
        // It should have 3 edges.
        assert_eq!(face.edges.len(), 3);
        let mut midpoints = Vec::new();
        println!("looping through face");
        for (edge, (p1,p2)) in face.edges.iter().zip(vec![(p2,p3),(p3,p1),(p1,p2)].into_iter()) {
            println!("Edge[]: {:?} -> {:?}", edge.point, edge.sym().point);
            assert_eq!(*edge.point, Some(p1));
            assert_eq!(*edge.sym().point, Some(p2));
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
        assert_eq!(centre, Point::new(2.5,5.0/3.0))
    }
}
