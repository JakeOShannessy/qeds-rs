#![feature(new_uninit)]
#![feature(ptr_wrapping_offset_from)]
use core::mem::MaybeUninit;

fn main() {
    println!("Edge Size: {:?}(0x{:x?})", std::mem::size_of::<Edge>(), std::mem::size_of::<Edge>());
    println!("Quad Size: {:?}(0x{:x?})", std::mem::size_of::<Quad>(), std::mem::size_of::<Quad>());
    println!("Quad Alignment: {:?}(0x{:x?})", std::mem::align_of::<Quad>(), std::mem::align_of::<Quad>());
    // Step 1. Create a Qeds data structure.
    let mut qeds = Qeds::new();
    // Step 2. Add some edges to it.
    let q1 = qeds.make_edge();
    let q2 = qeds.make_edge();
    println!("QEDS1: {:?}", qeds);
    for quad in qeds.quads.iter(){
        let quad: &Quad = unsafe {&**quad};
        println!("Quad: {:?}", quad);
        for edge in quad.edges.iter() {
            let e = edge;
            println!("Edge: {:?}", e);
            let e_rot = e.rot();
            println!("EdgeRot: {:?}", e_rot);
            let e_rot2 = e_rot.rot();
            println!("EdgeRotRot: {:?}", e_rot2);
            let e_rot3 = e_rot2.rot();
            println!("EdgeRotRotRot: {:?}", e_rot3);
            let e_rot4 = e_rot3.rot();
            println!("EdgeRotRotRotRot: {:?}", e_rot4);
            break;
        }
        break;
    }
    unsafe {
        let e1 = &mut (*q1).edges[0];
        let e2= &mut (*q2).edges[0];
        qeds.splice(e1,e2)
    }
    println!("QEDS2: {:?}", qeds);
}

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
    pub fn make_edge(&mut self) -> *mut Quad {
        let quad: Box<MaybeUninit<Quad>> = Box::new_zeroed();
        let mut quad = unsafe {quad.assume_init()};
        // The base edge e.
        quad.edges[0].next = &quad.edges[0];
        // eRot
        quad.edges[1].next = &quad.edges[1];
        // eSym
        quad.edges[2].next = &quad.edges[2];
        // eSymRot
        quad.edges[3].next = &quad.edges[3];
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
}

#[derive(Clone, Debug, PartialEq, PartialOrd)]
// TODO: adjust this based on pointer width.
#[repr(align(64))]
pub struct Quad {
    pub edges: [Edge; 4],
    // pub data: [Box<Point>; 4],
}

impl Quad {

}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct QuadRef<'a> {
    pub qeds: &'a Qeds,
    pub target: QuadTarget
}

// r and f can fit in one byte. If we limit the size of qeds we could fit it in
// the edge index.

impl<'a> QuadRef<'a> {

}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct QuadTarget {
    // The index of the Quad in the Qeds.
    pub e: usize,
    // Can only be 0, 1, 2, or 3.
    pub r: u8,
    // Can only be 0 or 1.
    pub f: u8,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Vertex {
    /// The index of the edge which starts at this vertex.
    pub edge_index: usize,
    // pub pos: Point,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Edge {
    // / The index of the face to the left of this edge.
    // pub face_index: usize,
    // pub index: usize,
    pub next: *const Edge,
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
        let ptr = if (offset1 == 0) || (offset1 > (4 + offset as usize)) {
            (self as *const Edge).wrapping_offset(offset)
        } else {
            let offset = offset - (offset1 as isize - 1);
            (self as *const Edge).wrapping_offset(offset)
        };
        unsafe {
            &*ptr
        }
    }
    #[inline(always)]
    fn offset_mut(&mut self, offset: isize) -> &mut Edge {
        let offset1 = (self as *const Edge).align_offset(std::mem::align_of::<Quad>());
        println!("offset1: {:?}", offset1);
        let ptr = if (offset1 == 0) || (offset1 > (4 + offset as usize)) {
            (self as *mut Edge).wrapping_offset(offset)
        } else {
            let offset = offset - (offset1 as isize - 1);
            (self as *mut Edge).wrapping_offset(offset)
        };
        unsafe {
            &mut *ptr
        }
    }


    pub fn onext(&self) -> &Edge {
        unsafe { &*self.next }
    }

    pub fn r_next(&self) -> &Edge {
        self.rot().onext().rot()
    }
}




#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Face {
    /// An index to any of the CCW oriented edges surrounding this face.
    pub edge_index: usize,
}

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Point {
    x: f64,
    y: f64,
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
        qeds.make_edge();
        // Step 3. Get the single edge in the data structure.
        let e: &Edge = unsafe { &(*qeds.quads[0]).edges[0] };
        // e with Rot applied 2 times is not the same edge as e
        assert_ne!(e.rot().rot() as *const Edge, e as *const Edge);
    }

    #[test]
    fn rot2_eq_sym() {
        // Step 1. Create a Qeds data structure.
        let mut qeds = Qeds::new();
        // Step 2. Add some edges to it.
        qeds.make_edge();
        // Step 3. Get the single edge in the data structure.
        let e: &Edge = unsafe { &(*qeds.quads[0]).edges[0] };
        // e with Rot applied 2 times is not the same edge as e
        assert_eq!(e.rot().rot() as *const Edge, e.sym() as *const Edge);
    }
}
