use core::ops::Add;
use core::ops::Sub;
use nalgebra::Matrix4;
use slab::Slab;
use std::{
    fmt::{Display, Formatter},
    ops::Mul,
};

#[derive(Clone, Debug)]
/// A Qeds data structure specialised to a 2d triangulation.
pub struct Triangulation {
    /// The quad-edge data structure we use as the basis for the triangulation.
    pub qeds: Qeds<Point, ()>,
}

impl Triangulation {
    pub fn new() -> Self {
        Self { qeds: Qeds::new() }
    }

    pub fn some_edge_a(&self) -> Option<EdgeRefA<Point, ()>> {
        let (i, _) = self.qeds.quads.iter().next()?;
        unsafe { Some(self.qeds.edge_a_ref(EdgeTarget::new(i, 0, 0))) }
    }

    // TODO: we need to be able to locate on an edge etc.
    // TODO: This has not yet been proved to be stable. It may also loop inifintely
    pub fn locate(&self, point: Point) -> EdgeRefA<Point, ()> {
        use rand::Rng;
        let mut e = self.some_edge_a().unwrap();
        let mut rng = rand::thread_rng();
        let mut current_iterations = 0;
        loop {
            current_iterations += 1;
            if current_iterations > 100 {
                panic!("locating failed");
            }
            // A random variable to determine if onext is tested first. If not
            // it is tested second.
            let onext_first: bool = rng.gen();
            if point == e.edge().point || point == e.sym().edge().point {
                return e;
            } else if e.lies_right(point) == Lies::Yes {
                e = e.sym();
            } else {
                if onext_first {
                    if e.d_prev().lies_right(point) != Lies::Yes {
                        e = e.d_prev();
                    } else if e.onext().lies_right(point) != Lies::Yes {
                        e = e.onext();
                    } else {
                        return e;
                    }
                } else {
                    if e.d_prev().lies_right(point) != Lies::Yes {
                        e = e.d_prev();
                    } else if e.onext().lies_right(point) != Lies::Yes {
                        e = e.onext();
                    } else {
                        return e;
                    }
                }
            }
        }
    }

    // pub fn locate(&self, point: Point) -> EdgeRefA<Point, ()> {
    //     let mut matches = Vec::new();
    //     for e1 in self.qeds.base_edges() {
    //         let e2 = e1.sym();
    //         if e1.in_left_face(point) {
    //             matches.push(e1.target());
    //         }
    //         if e2.in_left_face(point) {
    //             matches.push(e2.target());
    //         }
    //     }
    //     // Just return the first one for now. TODO: don't do this.
    //     unsafe {
    //         let e = self.qeds.edge_a_ref(*matches.get(0).expect("No triangle matches"));
    //         let face_points: Vec<(Point,Point)> = e.l_face().vertices_pairs().collect();
    //         assert_eq!(face_points.len(), 3);
    //         println!("{} lies in the face described by {:?}, as {}", point, face_points, e.in_left_face(point));
    //         for edge in e.l_face().edges.iter() {
    //             println!("PEdge: {} {}, lies left: {}", edge.edge().point, edge.sym().edge().point, edge.lies_left(point));
    //         }
    //         for edge in e.l_face().edges.iter() {
    //             println!("FEdge: {} {}, lies left: {}", edge.edge().point, edge.sym().edge().point, !edge.lies_left(point));
    //             // if !edge.lies_left(point) {
    //             //     return false;
    //             // }
    //         }
    //         e
    //     }
    // }

    pub fn add_point(&mut self, mut point: Point) {
        point.x *= 1.0e6;
        point.x = point.x.round();
        point.x *= 1.0e-6;

        point.y *= 1.0e6;
        point.y = point.y.round();
        point.y *= 1.0e-6;

        // TODO: consider: points lying on other points, on lines, and within triangles.
        if self.qeds.quads.len() == 0 {
            // If there are no edges, we can just add a point with cardinal
            // edges to infinity.
            let inf_north = Point::new(point.x, std::f64::INFINITY);
            let inf_south = Point::new(point.x, std::f64::NEG_INFINITY);
            let inf_east = Point::new(std::f64::INFINITY, point.y);
            let inf_west = Point::new(std::f64::NEG_INFINITY, point.y);

            let e_north = self.qeds.make_edge_with_a(point, inf_north).target();
            let e_east = self.qeds.make_edge_with_a(point, inf_east).target();
            let e_west = self.qeds.make_edge_with_a(point, inf_west).target();
            let e_south = self.qeds.make_edge_with_a(point, inf_south).target();

            unsafe {
                // Splice the finite ends together.
                self.qeds.splice(e_north, e_south);
                self.qeds.splice(e_north.sym(), e_south.sym());
                self.qeds.splice(e_west, e_north);
                self.qeds.splice(e_west.sym(), e_north.sym());
                self.qeds.splice(e_east, e_south);
                self.qeds.splice(e_east.sym(), e_south.sym());
                // Connect the infinite edges.
                self.qeds.connect(e_east, e_north.sym());
                self.qeds.connect(e_north, e_west.sym());
                self.qeds.connect(e_west, e_south.sym());
                self.qeds.connect(e_south, e_east.sym());
            }
        } else {
            // There are already edges in the triangulation, we know there must
            // be at least the original 4.
            let edge_target = self.locate(point).target();
            self.add_to_l_face(edge_target, point);
        }
    }

    fn add_to_l_face(&mut self, mut edge_target: EdgeTarget, point: Point) {
        unsafe {
            if self.qeds.edge_a_ref(edge_target).lies_right(point) == Lies::On {
                // todo!("We haven't yet considered on-line");
                let oprev = self.qeds.edge_a_ref(edge_target).oprev().target();
                let to_delete = edge_target;
                self.qeds.delete(to_delete);
                edge_target = oprev;
            }
            let first = self.qeds.edge_a_ref(edge_target).edge().point;
            let mut base = self.qeds.make_edge_with_a(first, point).target();
            self.qeds.splice(base, edge_target);
            loop {
                let base_ref = self.qeds.connect(edge_target, base.sym());
                edge_target = base_ref.oprev().target();
                base = base_ref.target();
                if self.qeds.edge_a(edge_target.sym()).point == first {
                    break;
                }
            }
            edge_target = self.qeds.edge_a_ref(base).oprev().target();
            // TODO: perform swapping
        }
    }

    // /// Assumes the EdgeTargets have already been put in the correct order.
    // pub unsafe fn add_external_point_ordered(
    //     &mut self,
    //     boundary_edges: Vec<EdgeTarget>,
    //     point: Point,
    // ) {
    //     let mut edges = boundary_edges.into_iter();
    //     let (g, e) = {
    //         let (a, b) = self.add_external_point(edges.next().unwrap(), point);
    //         (a.target(), b.target())
    //     };
    //     for f in edges {
    //         self.qeds.connect(e.sym(), f.sym());
    //     }
    //     self.boundary_edge = Some(g);
    // }

    // pub unsafe fn add_external_point_unordered(
    //     &mut self,
    //     mut boundary_edges: Vec<EdgeTarget>,
    //     point: Point,
    // ) {
    //     // Reorder the edges so that their base points are CCW around the point.
    //     boundary_edges.sort_by(|a, b| {
    //         let pa = self.qeds.edge_a_ref(*a).edge().point;
    //         let pb = self.qeds.edge_a_ref(*b).edge().point;
    //         match left_or_right(point, pb, pa) {
    //             Direction::Left => std::cmp::Ordering::Less,
    //             Direction::Straight => std::cmp::Ordering::Equal,
    //             Direction::Right => std::cmp::Ordering::Greater,
    //         }
    //     });
    //     self.add_external_point_ordered(boundary_edges, point);
    // }

    // pub fn boundary(&self) -> Option<BoundaryIter<Point, ()>> {
    //     if let Some(boundary_edge) = self.boundary_edge {
    //         Some(BoundaryIter::new(&self.qeds, boundary_edge))
    //     } else {
    //         None
    //     }
    // }

    // pub unsafe fn add_external_point(
    //     &mut self,
    //     boundary_edge: EdgeTarget,
    //     p: Point,
    // ) -> (EdgeRefA<Point, ()>, EdgeRefA<Point, ()>) {
    //     let dest = self
    //         .qeds
    //         .edge_a_ref(boundary_edge)
    //         .sym()
    //         .edge()
    //         .point
    //         .clone();
    //     let e = self.qeds.make_edge_with_a(p, dest).target();
    //     self.qeds.splice(e.sym(), boundary_edge.sym());
    //     let f = self.qeds.connect(boundary_edge.sym(), e).target();
    //     (self.qeds.edge_a_ref(f), self.qeds.edge_a_ref(e))
    // }

    /// Determine whether an Edge (i.e. two NavTris) satisfies the Delaunay
    /// criterion. An edge with only a single NavTri will return True.
    unsafe fn del_test(&self, e: EdgeTarget) -> bool {
        // Get the edge.
        let edge = self.qeds.edge_a_ref(e);
        // Get all of the vertices around this egdge in a CCW order.
        let a = edge.oprev().edge().point;
        let b = edge.r_prev().sym().edge().point;
        let c = edge.l_next().edge().point;
        let d = edge.onext().sym().edge().point;
        if !a.is_finite() && d.is_finite() && b.is_finite() {
            cross_product_2d(d - c, b - c) > 0.0
        } else if !c.is_finite() && d.is_finite() && b.is_finite() {
            cross_product_2d(b - a, d - a) > 0.0
        } else if !a.is_finite() || !b.is_finite() || !c.is_finite() || !d.is_finite() {
            false
        } else {
            del_test_ccw(a, b, c, d)
        }
    }

    // TODO: remove this
    pub fn retriangulate_all(&mut self) {
        let edge_targets: Vec<EdgeTarget> =
            self.qeds.base_edges().map(|edge| edge.target()).collect();
        for e in edge_targets.into_iter() {
            unsafe {
                if self.del_test(e) {
                    self.swap(e);
                }
            }
        }
    }

    pub unsafe fn swap(&mut self, e: EdgeTarget) {
        let a = self.qeds.edge_a_ref(e).oprev().target();
        let b = self.qeds.edge_a_ref(e).sym().oprev().target();

        self.qeds.splice(e, a);
        self.qeds.splice(e.sym(), b);

        let a_lnext = self.qeds.edge_a_ref(a).l_next().target();
        self.qeds.splice(e, a_lnext);

        let b_lnext = self.qeds.edge_a_ref(b).l_next().target();
        self.qeds.splice(e.sym(), b_lnext);

        let a_dest = self.qeds.edge_a_ref(a).sym().edge().point;
        let b_dest = self.qeds.edge_a_ref(b).sym().edge().point;
        self.qeds.edge_a_mut(e).point = a_dest;
        self.qeds.edge_a_mut(e.sym()).point = b_dest;
    }

    // /// If at any point the boundary loops back on itself, we know it to be
    // /// linear. Technically it only shows that one edge is, but one of our
    // /// invariant is that only one linear edge is not possible.
    // pub fn is_linear(&self) -> bool {
    //     for edge in self.boundary().unwrap() {
    //         if edge.sym().onext() == edge.sym() {
    //             return true;
    //         }
    //     }
    //     false
    // }
}

fn cross_product_2d(Point { x: x1, y: y1 }: Point, Point { x: x2, y: y2 }: Point) -> f64 {
    x1 * y2 - x2 * y1
}
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub struct Triangle {
    /// The minimum bounding point.
    pub min: Point,
    /// The maximum bounding point.
    pub max: Point,
    /// An edge target which lies on the triangle. It can be any one of the
    /// three that make up the triangle.
    pub target: EdgeTarget,
}

/// Determine whether a set of 4 points satisfies the Delaunay criterion. This
/// assumes that the pointes are sorted in a CCW fashion.
fn del_test_ccw(a: Point, b: Point, c: Point, d: Point) -> bool {
    let matrix: Matrix4<f64> = Matrix4::new(
        a.x,
        a.y,
        a.x.powi(2) + a.y.powi(2),
        1.0,
        b.x,
        b.y,
        b.x.powi(2) + b.y.powi(2),
        1.0,
        c.x,
        c.y,
        c.x.powi(2) + c.y.powi(2),
        1.0,
        d.x,
        d.y,
        d.x.powi(2) + d.y.powi(2),
        1.0,
    );
    matrix.determinant() > 0.0
}

pub struct BoundaryIter<'a, AData, BData> {
    qeds: &'a Qeds<AData, BData>,
    next_edge: Option<EdgeRefA<'a, AData, BData>>,
    first_edge_target: EdgeTarget,
}

impl<'a, AData, BData> BoundaryIter<'a, AData, BData> {
    fn new(qeds: &'a Qeds<AData, BData>, initial_boundary_edge: EdgeTarget) -> Self {
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
            let next_edge = unsafe { self.qeds.edge_a_ref(self.first_edge_target) };
            self.next_edge = Some(next_edge.r_prev());
            Some(next_edge)
        }
    }
}

/// This data structure is a single instance of a quad-edge data structure.
#[derive(Clone, Debug)]
pub struct Qeds<AData, BData> {
    /// The vector of quads which we index into.
    pub quads: Slab<Quad<AData, BData>>,
}

impl<AData: Default, BData: Default> Qeds<AData, BData> {
    /// Create an edge in the [`Qeds`].
    pub fn make_edge(&mut self) -> EdgeRefA<AData, BData> {
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
    pub fn make_edge_with_a(&mut self, org: AData, dest: AData) -> EdgeRefA<AData, BData> {
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
    ) -> EdgeRefA<AData, BData> {
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
    pub fn connect(&mut self, edge_a: EdgeTarget, edge_b: EdgeTarget) -> EdgeRefA<AData, BData> {
        unsafe {
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
}

impl<AData, BData> Qeds<AData, BData> {
    /// Create the simplest [`Qeds`] with a single [`Edge`] and a single
    /// [`Face`]. Covers the whole sphere.
    pub fn new() -> Self {
        Self { quads: Slab::new() }
    }

    pub unsafe fn edge_a(&self, target: EdgeTarget) -> &Edge<AData> {
        if target.r == 0 || target.r == 2 {
            &self.quads[target.e].edges_a[(target.r / 2) as usize]
        } else {
            unreachable!()
        }
    }

    pub unsafe fn edge_b(&self, target: EdgeTarget) -> &Edge<BData> {
        if target.r == 0 || target.r == 2 {
            unreachable!()
        } else {
            &self.quads[target.e].edges_b[(target.r / 2) as usize]
        }
    }

    pub unsafe fn edge_ref(&self, target: EdgeTarget) -> EdgeRefAny<AData, BData> {
        EdgeRefAny {
            qeds: &self,
            target,
        }
    }

    pub unsafe fn edge_a_ref(&self, target: EdgeTarget) -> EdgeRefA<AData, BData> {
        EdgeRefA {
            qeds: &self,
            target,
        }
    }

    pub unsafe fn edge_mut(&mut self, target: EdgeTarget) -> EdgeABMut<AData, BData> {
        if target.r == 0 || target.r == 2 {
            EdgeABMut::A(&mut self.quads[target.e].edges_a[(target.r / 2) as usize])
        } else {
            EdgeABMut::B(&mut self.quads[target.e].edges_b[(target.r / 2) as usize])
        }
    }

    pub unsafe fn edge_a_mut(&mut self, target: EdgeTarget) -> &mut Edge<AData> {
        if target.r == 0 || target.r == 2 {
            &mut self.quads[target.e].edges_a[(target.r / 2) as usize]
        } else {
            unreachable!()
        }
    }

    pub unsafe fn edge_b_mut(&mut self, target: EdgeTarget) -> &mut Edge<BData> {
        if target.r == 0 || target.r == 2 {
            unreachable!()
        } else {
            &mut self.quads[target.e].edges_b[(target.r / 2) as usize]
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

    pub unsafe fn delete(&mut self, e: EdgeTarget) {
        let oprev = self.edge_ref(e).oprev().target();
        self.splice(e, oprev);
        let sym_oprev = self.edge_ref(e).sym().oprev().target();
        self.splice(e.sym(), sym_oprev);
        // TODO: verify that this is correct
        self.quads.remove(e.e);
    }

    pub fn base_edges(&self) -> BaseEdgeIter<AData, BData> {
        BaseEdgeIter::new(self)
    }
}

pub struct BaseEdgeIter<'a, AData, BData> {
    qeds: &'a Qeds<AData, BData>,
    quad_iter: slab::Iter<'a, Quad<AData, BData>>,
}

impl<'a, AData, BData> BaseEdgeIter<'a, AData, BData> {
    fn new(qeds: &'a Qeds<AData, BData>) -> Self {
        Self {
            qeds: qeds,
            quad_iter: qeds.quads.iter(),
        }
    }
}

impl<'a, AData, BData> Iterator for BaseEdgeIter<'a, AData, BData> {
    type Item = EdgeRefA<'a, AData, BData>;
    fn next(&mut self) -> Option<Self::Item> {
        let quad = self.quad_iter.next()?;
        unsafe { Some(self.qeds.edge_a_ref(EdgeTarget::new(quad.0, 0, 0))) }
    }
}

#[derive(Clone, Debug, PartialEq, PartialOrd)]
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
    pub fn edge(&self) -> EdgeAB<AData, BData> {
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

impl<'a, BData> EdgeRefA<'a, Point, BData> {
    pub fn midpoint(&self) -> Point {
        let this_point = self.edge().point;
        let next_point = self.sym().edge().point;
        this_point.midpoint(next_point)
    }

    pub fn in_left_face(&self, point: Point) -> bool {
        for edge in self.l_face().edges.iter() {
            let is_left = edge.lies_left(point);
            if !is_left {
                return false;
            }
        }
        true
    }

    pub fn lies_right(&self, point: Point) -> Lies {
        use Lies::*;
        let pa = self.edge().point;
        let pb = self.sym().edge().point;
        // We have special treatment of infinite edges.
        if !pa.is_finite() && !pb.is_finite() {
            if pa.x == std::f64::INFINITY && pb.y == std::f64::INFINITY {
                // CCW
                No
            } else if pa.y == std::f64::INFINITY && pb.x == std::f64::NEG_INFINITY {
                // CCW
                No
            } else if pa.x == std::f64::NEG_INFINITY && pb.y == std::f64::NEG_INFINITY {
                // CCW
                No
            } else if pa.y == std::f64::NEG_INFINITY && pb.x == std::f64::INFINITY {
                // CCW
                No
            } else {
                // asssuming CW
                Yes
            }
        } else if pa.x == std::f64::INFINITY {
            if point.y == pb.y {
                On
            } else if point.y > pb.y {
                Yes
            } else {
                No
            }
        } else if pa.x == std::f64::NEG_INFINITY {
            if point.y == pb.y {
                On
            } else if point.y < pb.y {
                Yes
            } else {
                No
            }
        } else if pa.y == std::f64::INFINITY {
            if point.x == pb.x {
                On
            } else if point.x < pb.x {
                Yes
            } else {
                No
            }
        } else if pa.y == std::f64::NEG_INFINITY {
            if point.x == pb.x {
                On
            } else if point.x > pb.x {
                Yes
            } else {
                No
            }
        } else if pb.x == std::f64::INFINITY {
            if point.y == pa.y {
                On
            } else if point.y < pa.y {
                Yes
            } else {
                No
            }
        } else if pb.x == std::f64::NEG_INFINITY {
            if point.y == pa.y {
                On
            } else if point.y > pa.y {
                Yes
            } else {
                No
            }
        } else if pb.y == std::f64::INFINITY {
            if point.x == pb.x {
                On
            } else if point.x > pa.x {
                Yes
            } else {
                No
            }
        } else if pb.y == std::f64::NEG_INFINITY {
            if point.x == pb.x {
                On
            } else if point.x < pa.x {
                Yes
            } else {
                No
            }
        } else {
            match left_or_right(pa, pb, point) {
                Direction::Right => Yes,
                Direction::Straight => On,
                Direction::Left => No,
            }
        }
    }

    pub fn lies_left(&self, point: Point) -> bool {
        let pa = self.edge().point;
        let pb = self.sym().edge().point;
        lies_left(point, pa, pb)
    }
}

#[derive(Copy, Clone, Debug, Ord, PartialOrd, Eq, PartialEq, Hash)]
pub enum Lies {
    Yes,
    No,
    On,
}

fn lies_left(point: Point, pa: Point, pb: Point) -> bool {
    print!("Does {} lie to the left of: {} - {}", point, pa, pb);
    // We have special treatment of infinite edges.
    if !pa.is_finite() && !pb.is_finite() {
        print!(" infinite");
        if pa.x == std::f64::INFINITY && pb.y == std::f64::INFINITY {
            // CCW
            println!(" yes");
            true
        } else if pa.y == std::f64::INFINITY && pb.x == std::f64::NEG_INFINITY {
            // CCW
            println!(" yes");
            true
        } else if pa.x == std::f64::NEG_INFINITY && pb.y == std::f64::NEG_INFINITY {
            // CCW
            println!(" yes");
            true
        } else if pa.y == std::f64::NEG_INFINITY && pb.x == std::f64::INFINITY {
            // CCW
            println!(" yes");
            true
        } else {
            // asssuming CW
            println!(" no");
            false
        }
    } else if pa.x == std::f64::INFINITY {
        println!("a");
        point.y < pb.y
    } else if pa.x == std::f64::NEG_INFINITY {
        println!("b");
        point.y > pb.y
    } else if pa.y == std::f64::INFINITY {
        println!("c({} < {}): {}", point.x, pb.x, point.x < pb.x);
        point.x > pb.x
    } else if pa.y == std::f64::NEG_INFINITY {
        println!("d");
        point.x < pb.x
    } else if pb.x == std::f64::INFINITY {
        println!("e");
        point.y > pa.y
    } else if pb.x == std::f64::NEG_INFINITY {
        println!("f");
        point.y < pa.y
    } else if pb.y == std::f64::INFINITY {
        println!("g");
        point.x < pb.x
    } else if pb.y == std::f64::NEG_INFINITY {
        println!("h");
        point.x > pb.x
    } else {
        println!("..Direction: {:?}", left_or_right(pa, pb, point));
        left_or_right(pa, pb, point) == Direction::Left
    }
}

pub enum LineClassification {
    CCWInfinite,
    CWInfinite,
    Finite,
}

pub fn classify_segment(pa: Point, pb: Point) -> LineClassification {
    use LineClassification::*;
    const INF: f64 = std::f64::INFINITY;
    const NEG_INF: f64 = std::f64::NEG_INFINITY;
    if pa.is_finite() && pb.is_finite() {
        Finite
    } else {
        todo!()
    }
    // match (pa.x,pa.y,pb.x,pb.y) {
    //     (inf,_,_,inf) => CCWInfinite,
    //     (_,inf,neg_inf,_) => CCWInfinite,
    //     (neg_inf,_,_,neg_inf,) => CCWInfinite,
    //     (_,neg_inf,inf,_) => CCWInfinite,

    //     _ => Finite
    // }
}

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
    pub fn l_next(&self) -> Self {
        self.rot().rot().rot().onext().rot()
    }

    pub fn l_face(&self) -> Face<AData, BData> {
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
    pub fn r_face(&self) -> Face<AData, BData> {
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

    pub fn l_face(&self) -> FaceB<AData, BData> {
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
    pub fn r_face(&self) -> FaceB<AData, BData> {
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

#[derive(Debug, PartialEq, PartialOrd)]
pub enum EdgeABMut<'a, AData, BData> {
    A(&'a mut Edge<AData>),
    B(&'a mut Edge<BData>),
}

impl<'a, AData, BData> EdgeABMut<'a, AData, BData> {
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
pub struct Face<'a, AData, BData> {
    /// An index to any of the CCW oriented edges surrounding this face.
    pub edges: Vec<EdgeRefA<'a, AData, BData>>,
}

impl<BData: Clone> Face<'_, Point, BData> {
    pub fn centroid(&self) -> Option<Point> {
        let signed_area = self.signed_area();
        if signed_area <= 0.0 {
            return None;
        }
        let mut xsum = 0.0;
        let mut ysum = 0.0;
        for (pa, pb) in self.vertices_pairs() {
            let pa = pa;
            let pb = pb;
            xsum += (pa.x + pb.x) * (pa.x * pb.y - pb.x * pa.y);
            ysum += (pa.y + pb.y) * (pa.x * pb.y - pb.x * pa.y);
        }
        Some(Point {
            x: xsum / 6.0 / signed_area,
            y: ysum / 6.0 / signed_area,
        })
    }

    pub fn signed_area(&self) -> f64 {
        let mut sum = 0.0;
        for (pa, pb) in self.vertices_pairs() {
            let pa = pa;
            let pb = pb;
            sum += pa.x * pb.y - pb.x * pa.y;
        }
        0.5 * sum
    }
    pub fn vertices(&self) -> FaceVerticesIter<BData> {
        FaceVerticesIter::new(self)
    }
    pub fn vertices_pairs(
        &self,
    ) -> std::iter::Zip<
        FaceVerticesIter<'_, BData>,
        std::iter::Skip<std::iter::Cycle<FaceVerticesIter<'_, BData>>>,
    > {
        self.vertices().zip(self.vertices().cycle().skip(1))
    }
}

#[derive(Clone, Debug)]
pub struct FaceB<'a, AData, BData> {
    /// An index to any of the CCW oriented edges surrounding this face.
    pub edges: Vec<EdgeRefB<'a, AData, BData>>,
}

#[derive(Clone, Debug)]
pub struct FaceVerticesIter<'a, BData> {
    face: &'a Face<'a, Point, BData>,
    next_index: usize,
}

impl<'a, BData> FaceVerticesIter<'a, BData> {
    fn new(face: &'a Face<'a, Point, BData>) -> Self {
        Self {
            face,
            next_index: 0,
        }
    }
}

impl<'a, BData> Iterator for FaceVerticesIter<'a, BData> {
    type Item = Point;
    fn next(&mut self) -> Option<Self::Item> {
        if self.next_index >= self.face.edges.len() {
            None
        } else {
            let edge = self.face.edges.get(self.next_index)?;
            self.next_index += 1;
            Some(edge.edge().point)
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

    pub fn distance(&self, other: Self) -> f64 {
        ((other.x - self.x).powi(2) + (other.y - self.y).powi(2)).sqrt()
    }
}

impl Display for Point {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "({},{})", self.x, self.y)
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
    use std::collections::HashSet;

    #[test]
    fn bad_collinear() {
        let b = Point::new(-10000.0, -1.0);
        let c = Point::new(-9000.0, -1.0);
        let d = Point::new(std::f64::MIN, -1.0 - std::f64::EPSILON);
        let p = Point::new(10000.0, -1.0);
        // This is a problem as it means we can have a triangle with both sides
        // collinear with another point. This indicates that our collinearity
        // test is insufficient.
        assert_ne!(b, d);
        assert_eq!(left_or_right(b, c, p), Direction::Straight);
        assert_eq!(left_or_right(d, c, p), Direction::Straight);
    }

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
    fn segment() {
        // Step 1. Create a Qeds data structure.
        let mut qeds: Qeds<Point, ()> = Qeds::new();
        let point_a = Point::new(0.0, 0.0);
        let point_b = Point::new(5.0, 0.0);
        // Step 2. Add some edges to it.
        let q1 = qeds.make_edge_with_a(point_a, point_b).target;
        println!("q1 Target: {:?}", q1);
        unsafe {
            let e = &qeds.quads[0].edges_a[0];
            // e with sym applied 2 times is the same edge as e
            assert_eq!(
                qeds.edge_a_ref(q1).edge() as *const Edge<Point>,
                e as *const Edge<Point>
            );
            assert_eq!(qeds.edge_a_ref(q1).sym().sym().edge(), e);
            // The Org of the edge is a point and it is point_a
            assert_eq!(qeds.edge_a_ref(q1).edge().point, point_a);
            // The Dest of the edge (eSymOrg) is point_b
            println!("aa: {:?}", qeds.edge_a_ref(q1).target());
            println!("ab: {:?}", qeds.edge_a_ref(q1).sym().target());
            assert_eq!(qeds.edge_a_ref(q1).sym().edge().point, point_b);
        }
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
        let mut qeds: Qeds<(), ()> = Qeds::new();
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
    fn single_splice() {
        // Step 1. Create a Qeds data structure.
        let mut qeds: Qeds<(), ()> = Qeds::new();
        // Step 2. Add some edges to it.
        let a = qeds.make_edge().target;
        let b = qeds.make_edge().target;

        // Step 3. Splice those edges together so that we actually have
        // something of a network.
        unsafe {
            qeds.splice(a.sym(), b);
        }
    }

    #[test]
    fn triangle_face() {
        // Step 1. Create a Qeds data structure.
        let mut qeds: Qeds<Point, ()> = Qeds::new();
        // Create some point to add geometry.
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(5.0, 0.0);
        let p3 = Point::new(2.5, 5.0);
        // Step 2. Add some edges to it.
        let q1 = qeds.make_edge_with_a(p1, p2).target;
        let q2 = qeds.make_edge_with_a(p2, p3).target;

        unsafe {
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

            assert_eq!(edge.edge().point, p1);
            assert_eq!(edge.sym().edge().point, p2);
            assert_eq!(edge.l_next().edge().point, p2);
            assert_eq!(edge.l_next().sym().edge().point, p3);
            assert_eq!(edge.l_next().l_next().edge().point, p3);
            assert_eq!(edge.l_next().l_next().sym().edge().point, p1);
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
                assert_eq!(edge.edge().point, p1);
                assert_eq!(edge.sym().edge().point, p2);
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

    fn has_edge(triangulation: &Triangulation, pa: Point, pb: Point) -> bool {
        get_edge(triangulation, pa, pb).is_some()
    }

    fn get_edge(
        triangulation: &Triangulation,
        pa: Point,
        pb: Point,
    ) -> Option<EdgeRefA<'_, Point, ()>> {
        for (i, quad) in triangulation.qeds.quads.iter() {
            let edge1 = quad.edges_a[0].point;
            let edge2 = quad.edges_a[1].point;
            let edge_target = if edge1 == pa && edge2 == pb {
                EdgeTarget::new(i, 0, 0)
            } else if edge1 == pb && edge2 == pa {
                EdgeTarget::new(i, 2, 0)
            } else {
                continue;
            };
            return Some(unsafe { triangulation.qeds.edge_a_ref(edge_target) });
        }
        None
    }

    #[test]
    fn one_point_triangulation_only() {
        let mut triangulation = Triangulation::new();
        let p1 = Point::new(0.0, 0.0);
        triangulation.add_point(p1);
        assert_eq!(triangulation.qeds.quads.len(), 8, "wrong number of quads");
        // Test that it has all four of the correct edges
        assert!(
            has_edge(
                &triangulation,
                Point::new(0.0, 0.0),
                Point::new(0.0, std::f64::INFINITY)
            ),
            "missing north edge"
        );
        assert!(
            has_edge(
                &triangulation,
                Point::new(0.0, 0.0),
                Point::new(0.0, std::f64::NEG_INFINITY)
            ),
            "missing south edge"
        );
        assert!(
            has_edge(
                &triangulation,
                Point::new(0.0, 0.0),
                Point::new(std::f64::INFINITY, 0.0)
            ),
            "missing east edge"
        );
        assert!(
            has_edge(
                &triangulation,
                Point::new(0.0, 0.0),
                Point::new(std::f64::NEG_INFINITY, 0.0)
            ),
            "missing west edge"
        );
        let north_edge = get_edge(
            &triangulation,
            Point::new(0.0, 0.0),
            Point::new(0.0, std::f64::INFINITY),
        )
        .unwrap();

        assert_eq!(north_edge.edge().point, Point::new(0.0, 0.0));
        assert_eq!(
            north_edge.sym().edge().point,
            Point::new(0.0, std::f64::INFINITY)
        );

        assert_eq!(
            north_edge.sym().onext().edge().point,
            Point::new(0.0, std::f64::INFINITY)
        );
        assert_eq!(
            north_edge.sym().onext().sym().edge().point,
            Point::new(std::f64::INFINITY, 0.0)
        );
    }

    #[test]
    fn lies_left_1() {
        assert!(!lies_left(
            Point::new(1.0, -1.0),
            Point { x: 0.0, y: 0.0 },
            Point {
                x: std::f64::INFINITY,
                y: 0.0
            }
        ));
    }

    #[test]
    fn one_point_triangulation_location() {
        let mut triangulation = Triangulation::new();
        let p1 = Point::new(0.0, 0.0);
        triangulation.add_point(p1);
        {
            // North-East
            let edge = triangulation.locate(Point::new(1.0, 1.0));
            let face_verts: Vec<Point> = edge.l_face().vertices().collect();
            assert!(face_verts
                .iter()
                .find(|x| **x == Point::new(std::f64::INFINITY, p1.y))
                .is_some());
            assert!(face_verts
                .iter()
                .find(|x| **x == Point::new(p1.x, std::f64::INFINITY))
                .is_some());
            assert!(face_verts.iter().find(|x| **x == p1).is_some());
        }
        {
            // North-West
            let edge = triangulation.locate(Point::new(-1.0, 1.0));
            let face_verts: Vec<Point> = edge.l_face().vertices().collect();
            assert!(face_verts
                .iter()
                .find(|x| **x == Point::new(std::f64::NEG_INFINITY, p1.y))
                .is_some());
            assert!(face_verts
                .iter()
                .find(|x| **x == Point::new(p1.x, std::f64::INFINITY))
                .is_some());
            assert!(face_verts.iter().find(|x| **x == p1).is_some());
        }
        {
            // South-East
            let point = Point::new(1.0, -1.0);
            let edge = triangulation.locate(point);
            let face_verts: Vec<Point> = edge.l_face().vertices().collect();
            assert!(!lies_left(
                point,
                Point { x: 0.0, y: 0.0 },
                Point {
                    x: std::f64::INFINITY,
                    y: 0.0
                }
            ));

            for edge in edge.l_face().edges.iter() {
                println!(
                    "Edge: {} {}, lies left: {}",
                    edge.edge().point,
                    edge.sym().edge().point,
                    edge.lies_left(point)
                );
                // if !self.lies_left(point) {
                //     return false;
                // }
            }

            assert!(face_verts
                .iter()
                .find(|x| **x == Point::new(std::f64::INFINITY, p1.y))
                .is_some());
            assert!(face_verts
                .iter()
                .find(|x| **x == Point::new(p1.x, std::f64::NEG_INFINITY))
                .is_some());
            assert!(face_verts.iter().find(|x| **x == p1).is_some());
        }
        {
            // South-West
            let edge = triangulation.locate(Point::new(-1.0, -1.0));
            let face_verts: Vec<Point> = edge.l_face().vertices().collect();
            assert!(face_verts
                .iter()
                .find(|x| **x == Point::new(std::f64::NEG_INFINITY, p1.y))
                .is_some());
            assert!(face_verts
                .iter()
                .find(|x| **x == Point::new(p1.x, std::f64::NEG_INFINITY))
                .is_some());
            assert!(face_verts.iter().find(|x| **x == p1).is_some());
        }
    }

    #[test]
    fn two_point_triangulation() {
        let mut triangulation = Triangulation::new();
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(1.0, 1.0);
        triangulation.add_point(p1);
        for (i, edge) in triangulation.qeds.base_edges().enumerate() {
            println!(
                "Edge[{}]: {}-{}",
                i,
                edge.edge().point,
                edge.sym().edge().point
            );
        }
        triangulation.add_point(p2);
        println!("Triangulation: {:?}", triangulation);
        for (i, edge) in triangulation.qeds.base_edges().enumerate() {
            println!(
                "Edge[{}]: {}-{}",
                i,
                edge.edge().point,
                edge.sym().edge().point
            );
        }
        assert_eq!(triangulation.qeds.quads.len(), 11);
        assert!(has_edge(&triangulation, p1, p2), "missing edge");
        assert!(
            has_edge(&triangulation, p2, Point::new(0.0, std::f64::INFINITY)),
            "missing edge"
        );
        assert!(
            has_edge(&triangulation, p2, Point::new(std::f64::INFINITY, 0.0)),
            "missing edge"
        );
    }

    #[test]
    fn triangle_triangulation() {
        let mut triangulation = Triangulation::new();
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(5.0, 0.10);
        let p3 = Point::new(2.5, 5.1);
        triangulation.add_point(p1);
        triangulation.add_point(p2);
        triangulation.add_point(p3);
    }

    fn valid_triangulation(triangulation: &Triangulation) {
        for (i, quad) in triangulation.qeds.quads.iter() {
            let target1 = EdgeTarget::new(i, 0, 0);
            let target2 = EdgeTarget::new(i, 2, 0);
            unsafe {
                assert_ne!(
                    triangulation.qeds.edge_a_ref(target1).onext().target(),
                    target1,
                    "Edge {:?} with points {} - {} has a floating end",
                    target1,
                    triangulation.qeds.edge_a_ref(target1).edge().point,
                    triangulation.qeds.edge_a_ref(target1).sym().edge().point
                );
                assert_ne!(
                    triangulation.qeds.edge_a_ref(target2).onext().target(),
                    target2,
                    "Edge {:?} with points {} - {} has a floating end",
                    target2,
                    triangulation.qeds.edge_a_ref(target2).edge().point,
                    triangulation.qeds.edge_a_ref(target2).sym().edge().point
                );
            }
        }
    }

    #[test]
    fn detailed_triangulation() {
        let mut triangulation = Triangulation::new();
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(10.0, 1.0);
        let p3 = Point::new(5.0, 10.0);
        let p4 = Point::new(11.0, 11.0);
        let p5 = Point::new(15.0, 2.5);
        let p6 = Point::new(12.0, 12.5);
        let p7 = Point::new(15.0, 12.5);
        triangulation.add_point(p1);
        triangulation.add_point(p2);
        triangulation.add_point(p3);
        valid_triangulation(&triangulation);
        // triangulation.add_point(p4);
    }

    #[test]
    fn spiral_triangulation() {
        let mut triangulation = Triangulation::new();
        // triangulation.add_point(Point::new(0.0, 0.0));
        let a = 1.0;
        let n_max = 100;
        for n in 0..n_max {
            let t = (n as f64 / 6.0) * std::f64::consts::PI;
            let x = a * t * t.cos();
            let y = a * t * t.sin();
            triangulation.add_point(Point::new(x, y));
        }
        valid_triangulation(&triangulation);
    }

    #[test]
    fn line_triangulation() {
        let mut triangulation = Triangulation::new();
        let p1 = Point::new(0.0, 0.0);
        triangulation.add_point(p1);
        {
            let point = Point::new(1.0, 0.0);
            println!("locating");
            let edge = triangulation.locate(point);
            println!("located");
            println!(
                "OnEdge[{:?}]: {} - {}",
                edge.lies_right(point),
                edge.edge().point,
                edge.sym().edge().point
            );
        }
        // assert_eq!(triangulation.qeds.quads.len(), 0, "wrong number of quads");
        let p2 = Point::new(5.0, 0.0);
        let p3 = Point::new(7.0, 0.0);
        triangulation.add_point(p2);
        triangulation.add_point(p3);
        for edge in triangulation.qeds.base_edges() {
            println!(
                "Base Edge: {} - {}",
                edge.edge().point,
                edge.sym().edge().point
            );
        }
    }

    #[test]
    fn problem_triangulation_swap() {
        let mut triangulation = Triangulation::new();
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(10.0, 0.0);
        let p3 = Point::new(5.0, 10.0);
        triangulation.add_point(p1);
        triangulation.add_point(p2);
        triangulation.add_point(p3);
        let e = get_edge(&triangulation, p3, Point::new(0.0, std::f64::INFINITY))
            .expect("no expected edge");
        println!("e: {:?} - del_test: {}", e.target(), unsafe {
            triangulation.del_test(e.target())
        });
        triangulation.retriangulate_all();
    }

    #[test]
    fn simple_triangulation_swap() {
        let mut triangulation = Triangulation::new();
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(-1.0, 1.0);
        let p3 = Point::new(1.0, 1.0);
        triangulation.add_point(p1);
        triangulation.add_point(p2);
        triangulation.add_point(p3);
        let e = get_edge(&triangulation, p1, Point::new(0.0, std::f64::INFINITY))
            .expect("no expected edge");
        println!("Cross Product: {}", cross_product_2d(p3 - p1, p2 - p1));
        println!("e: {:?} - del_test: {}", e.target(), unsafe {
            triangulation.del_test(e.target())
        });
        triangulation.retriangulate_all();
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
        let mut qeds: Qeds<Point, ()> = Qeds::new();
        // Step 2. Add the first two edges.
        let q1 = qeds.make_edge_with_a(p1, p2).target;
        let q2 = qeds.make_edge_with_a(p2, p3).target;
        let q4 = qeds.make_edge_with_a(p4, p3).target;

        // Step 4. Splice those edges together so that we actually have
        // something of a network. This adds the third edge.
        unsafe {
            qeds.splice(q1.sym(), q2);
            assert_eq!(qeds.edge_a_ref(q1).sym().onext().target, q2);
            let e = qeds.connect(q2, q1).target();
            // At this point the triangle has been created.
            {
                assert_eq!(qeds.edge_a_ref(q1).l_next().target, q2);
                assert_eq!(qeds.edge_a(q1).point, p1);
                assert_eq!(qeds.edge_a(q1.sym()).point, p2);
                assert_eq!(qeds.edge_a_ref(q1).l_next().edge().point, p2);
                assert_eq!(qeds.edge_a_ref(q1).l_next().sym().edge().point, p3);
                assert_eq!(qeds.edge_a_ref(q1).l_next().l_next().edge().point, p3);
                assert_eq!(qeds.edge_a_ref(EdgeTarget::new(3, 0, 0)).edge().point, p3);
                assert_eq!(
                    qeds.edge_a_ref(EdgeTarget::new(3, 0, 0)).sym().edge().point,
                    p1
                );
                assert_eq!(qeds.edge_a_ref(q1).l_next().l_next().sym().edge().point, p1);
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
                edge.edge() as *const Edge<Point>,
                edge.edge().point,
                edge.sym().edge().point
            );
            println!(
                "Edge2[{:?}]: {:?} -> {:?}",
                edge.l_next().edge() as *const Edge<Point>,
                edge.l_next().edge().point,
                edge.l_next().sym().edge().point
            );
            println!(
                "Edge3[{:?}]: {:?} -> {:?}",
                edge.l_next().l_next().edge() as *const Edge<Point>,
                edge.l_next().l_next().edge().point,
                edge.l_next().l_next().sym().edge().point
            );
            println!(
                "Edge4[{:?}]: {:?} -> {:?}",
                edge.l_next().l_next().l_next().edge() as *const Edge<Point>,
                edge.l_next().l_next().l_next().edge().point,
                edge.l_next().l_next().l_next().sym().edge().point
            );

            assert_eq!(edge.l_next().target(), q2);
            assert_eq!(edge.edge().point, p1);
            assert_eq!(edge.sym().edge().point, p2);
            assert_eq!(edge.l_next().edge().point, p2);
            assert_eq!(edge.l_next().sym().edge().point, p3);
            assert_eq!(edge.l_next().l_next().edge().point, p3);
            assert_eq!(edge.l_next().l_next().sym().edge().point, p1);

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
                assert_eq!(edge.edge().point, p1);
                assert_eq!(edge.sym().edge().point, p2);
                midpoints.push(edge.midpoint());
            }
            let mut centre = Point::new(0.0, 0.0);
            let n = midpoints.len();
            println!("midpoints: {:?}", midpoints);
            for p in midpoints.into_iter() {
                centre.x += p.x;
                centre.y += p.y;
            }
            centre.x = centre.x / (n as f64);
            centre.y = centre.y / (n as f64);
            assert_eq!(centre, Point::new(2.5, 5.0 / 3.0))
        }
    }
}
