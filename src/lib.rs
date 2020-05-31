#![feature(new_uninit)]
#![feature(ptr_wrapping_offset_from)]
use core::mem::MaybeUninit;
use core::ops::Add;
use core::ops::Sub;

// fn main() {
//     println!("Edge Size: {:?}(0x{:x?})", std::mem::size_of::<Edge>(), std::mem::size_of::<Edge>());
//     println!("Quad Size: {:?}(0x{:x?})", std::mem::size_of::<Quad>(), std::mem::size_of::<Quad>());
//     println!("Quad Alignment: {:?}(0x{:x?})", std::mem::align_of::<Quad>(), std::mem::align_of::<Quad>());
//     // let quad_layout = std::alloc::Layout::from_size_align(std::mem::size_of::<Quad>(), 64).unwrap();
//     // let quad_layout_standard = std::alloc::Layout::from_size_align(std::mem::size_of::<Quad>(), 64).unwrap();
//     // println!("Quad Layout: {:?}", quad_layout);
//     // println!("Quad Layout Std: {:?}",quad_layout_standard);

//     // Step 1. Create a Qeds data structure.
//     let mut qeds = Qeds::new();
//     // Step 2. Add some edges to it.
//     let q1 = qeds.make_edge();
//     let q2 = qeds.make_edge();
//     println!("QEDS1: {:?}", qeds);
//     for quad in qeds.quads.iter(){
//         let quad: &Quad = unsafe {&**quad};
//         println!("Quad: {:?}", quad);
//         for edge in quad.edges.iter() {
//             let e = edge;
//             println!("Edge: {:?}", e);
//             let e_rot = e.rot();
//             println!("EdgeRot: {:?}", e_rot);
//             let e_rot2 = e_rot.rot();
//             println!("EdgeRotRot: {:?}", e_rot2);
//             let e_rot3 = e_rot2.rot();
//             println!("EdgeRotRotRot: {:?}", e_rot3);
//             let e_rot4 = e_rot3.rot();
//             println!("EdgeRotRotRotRot: {:?}", e_rot4);
//             break;
//         }
//         break;
//     }
//     unsafe {
//         let e1 = &mut (*q1).edges[0];
//         let e2= &mut (*q2).edges[0];
//         qeds.splice(e1,e2)
//     }
//     println!("QEDS2: {:?}", qeds);
// }

/// This data structure is a single instance of a quad-edge data structure.
/// Rather than using pointers it uses offsets into a vectors.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Qeds {
    pub quads: Vec<*mut Quad>,
    // pub faces: Vec<Face>,
}

impl Qeds {
    /// Create the simplest [`Qeds`] with a single [`Edge`] and a single [`Face`].
    /// Covers the whole sphere.
    pub fn new() -> Self {
        let quads = Vec::new();
        Self { quads }
    }

    /// Create an edge in the [`Qeds`].
    pub fn make_edge(&mut self, point: Point) -> *mut Quad {
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
        quad.edges[0].point = Box::new(Some(point));
        // eRot
        quad.edges[1].next = &quad.edges[3];
        quad.edges[1].point = Box::new(None);
        // eSym
        quad.edges[2].next = &quad.edges[2];
        quad.edges[2].point = Box::new(None);
        // eSymRot
        quad.edges[3].next = &quad.edges[2];
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
        // println!("offset(d): {:?}", d);
        let ptr = (self as *const Edge).wrapping_offset(d);
        unsafe { &*ptr }
    }
    #[inline(always)]
    fn offset_mut(&mut self, offset: isize) -> &mut Edge {
        let ptr: *const Edge = self.offset(offset);
        unsafe { &mut *(ptr as *mut Edge) }
    }

    pub fn onext(&self) -> &Edge {
        unsafe { &*self.next }
    }

    pub fn r_next(&self) -> &Edge {
        self.rot().onext().rot()
    }

    pub fn l_face(&self) -> Face {
        // Save the next edge of this face so that we know when to end the loop.
        let first_next_edge = self.sym().r_next();
        let mut edges = Vec::new();
        // Loop around collecting the edges.
        let mut current_next_edge = self.sym().r_next();
        loop {
            edges.push(current_next_edge);
            current_next_edge = current_next_edge.sym().r_next();
            if current_next_edge == first_next_edge {
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
        qeds.make_edge(Point::new(0.0, 0.0));
        // Step 3. Get the single edge in the data structure.
        let e: &Edge = unsafe { &(*qeds.quads[0]).edges[0] };
        // e with Rot applied 4 times is the same edge as e
        assert_eq!(e.rot().rot().rot().rot() as *const Edge, e as *const Edge);
    }

    #[test]
    fn e3() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        // Step 2. Add some edges to it.
        qeds.make_edge(Point::new(0.0, 0.0));
        // Step 3. Get the single edge in the data structure.
        let e: &Edge = unsafe { &(*qeds.quads[0]).edges[0] };
        // e with Rot applied 2 times is not the same edge as e
        assert_ne!(e.rot().rot() as *const Edge, e as *const Edge);
    }

    #[test]
    fn e4() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        let point_a = Point::new(0.0, 0.0);
        let point_b = Point::new(5.0, 0.0);
        // Step 2. Add some edges to it.
        let q1 = qeds.make_edge(point_a);
        let q2 = qeds.make_edge(point_b);
        unsafe {
            qeds.join(&mut (*q1).edges[2], &mut (*q2).edges[0]);
        }

        // Step 3. Get the single edge in the data structure.
        let e: &Edge = unsafe { &(*qeds.quads[0]).edges[0] };
        let e2: &Edge = unsafe { &(*qeds.quads[1]).edges[0] };
        // e with sym applied 2 times is not the same edge as e
        assert_eq!(e.sym().sym() as *const Edge, e as *const Edge);
        assert_eq!((*e.point).unwrap(), point_a);
        assert_eq!((*e2.point).unwrap(), point_b);
        assert_eq!((*e.sym().r_next().point).unwrap(), point_b);
        assert_eq!((*e.sym().point).unwrap(), point_b);
    }


    #[test]
    fn rot2_eq_sym() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        // Step 2. Add some edges to it.
        qeds.make_edge(Point::new(0.0, 0.0));
        // Step 3. Get the single edge in the data structure.
        let e: &Edge = unsafe { &(*qeds.quads[0]).edges[0] };
        // e with Rot applied 2 times is not the same edge as e
        assert_eq!(e.rot().rot() as *const Edge, e.sym() as *const Edge);
    }

    #[test]
    fn triangle_face() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        // Step 2. Add some edges to it.
        let q1 = qeds.make_edge(Point::new(0.0, 0.0));
        let q2 = qeds.make_edge(Point::new(5.0, 0.0));
        let q3 = qeds.make_edge(Point::new(2.5, 5.0));

        // Step 3. Splice those edges together so that we actually have
        // something of a network.
        unsafe {
            qeds.splice(&mut (*q1).edges[2], &mut (*q2).edges[0]);
            qeds.splice(&mut (*q2).edges[2], &mut (*q3).edges[0]);
            qeds.splice(&mut (*q3).edges[2], &mut (*q1).edges[0]);
        }
        // Get the first edge.
        let quad = qeds.quads.iter().next().unwrap();
        let quad: &Quad = unsafe { &**quad };
        let edge = &quad.edges[0];
        // Get the face from it.
        let face = edge.l_face();
        // It should have 3 edges.
        assert_eq!(face.edges.len(), 3);
    }

    #[test]
    fn triangle_face_centre() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        // Step 2. Add some edges to it.
        let q1 = qeds.make_edge(Point::new(0.0, 0.0));
        let q2 = qeds.make_edge(Point::new(5.0, 0.0));
        let q3 = qeds.make_edge(Point::new(2.5, 5.0));

        // Step 3. Splice those edges together so that we actually have
        // something of a network.
        unsafe {
            qeds.join(&mut (*q1).edges[2], &mut (*q2).edges[0]);
            qeds.join(&mut (*q2).edges[2], &mut (*q3).edges[0]);
            qeds.join(&mut (*q3).edges[2], &mut (*q1).edges[0]);
        }
        // Get the first edge.
        let quad = qeds.quads.iter().next().unwrap();
        let quad: &Quad = unsafe { &**quad };
        let edge = &quad.edges[0];
        // Get the face from it.
        let face = edge.l_face();
        // It should have 3 edges.
        assert_eq!(face.edges.len(), 3);
        let mut midpoints = Vec::new();
        for edge in face.edges.iter() {
            midpoints.push(edge.midpoint());
            // println!("edge: {:?}, {:?}, {:?}", edge, edge.sym(), edge.midpoint());
        }
        let mut centre = Point::new(0.0, 0.0);
        let n = midpoints.len();
        for p in midpoints.into_iter() {
            centre.x += p.unwrap().x;
            centre.y += p.unwrap().y;
        }
        centre.x = centre.x / (n as f64);
        centre.y = centre.y / (n as f64);
        assert_eq!(centre, Point::new(2.5,5.0/3.0))
    }
}
