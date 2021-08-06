use crate::point::*;
use crate::qeds::*;
use nalgebra::Matrix4;

#[derive(Copy, Clone, Debug, Ord, PartialOrd, Eq, PartialEq, Hash)]
pub enum Intersection {
    /// Contains the edge that is intersectect.
    Edge(EdgeTarget),
    /// Contains the edge that has it's origin at the point that is intersected.
    Point(EdgeTarget),
}

impl Intersection {
    fn target(self) -> EdgeTarget {
        match self {
            Intersection::Edge(e) => e,
            Intersection::Point(e) => e,
        }
    }
}

pub fn edge_from_target<T: Clone>(
    target: EdgeTarget,
    triangulation: &SurfaceTriangulation<T>,
) -> Edge<Segment<T>> {
    let edge_ref = unsafe { triangulation.qeds.edge_a_ref(target) };
    edge_ref.edge().clone()
}

pub fn triangle_across_from_target<T>(
    target: EdgeTarget,
    triangulation: &SurfaceTriangulation<T>,
) -> NodeTarget {
    let edge_ref = unsafe { triangulation.qeds.edge_a_ref(target) };
    let triangle_across = edge_ref.triangle_across();
    NodeTarget(triangle_across.target())
}

/// A NodeTarget is an EdgeTarget that is guaranteed to be a canonical triangle
/// EdgeTarget.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeTarget(EdgeTarget);

impl NodeTarget {
    pub fn new_unchecked(target: EdgeTarget) -> Self {
        NodeTarget(target)
    }

    pub fn as_edge(self) -> EdgeTarget {
        self.0
    }

    pub fn get_edge<T>(self, triangulation: &SurfaceTriangulation<T>, i: u8) -> EdgeTarget {
        let edges = self.edges(triangulation);
        edges[i as usize]
    }

    pub fn edges<T>(self, triangulation: &SurfaceTriangulation<T>) -> Vec<EdgeTarget> {
        let edge_ref = unsafe { triangulation.qeds.edge_a_ref(self.into()) };
        edge_ref
            .l_face()
            .edges
            .into_iter()
            .map(|e| e.target())
            .collect()
    }

    pub fn adjacent_tris<T>(self, triangulation: &SurfaceTriangulation<T>) -> Vec<NodeTarget> {
        let edges = self.edges(triangulation);
        let mut adjacents = Vec::new();
        for edge in edges {
            let tri = triangle_across_from_target(edge, triangulation);
            adjacents.push(tri);
        }
        adjacents
    }
}

impl From<NodeTarget> for EdgeTarget {
    fn from(nt: NodeTarget) -> Self {
        nt.0
    }
}

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Segment<T> {
    /// The point at the origin of this edge.
    pub point: Point,
    pub data: T,
}

impl<T> Segment<T> {
    pub fn new(point: Point, data: T) -> Self {
        Self { point, data }
    }
}

impl<T> HasPoint for Segment<T> {
    fn point(&self) -> Point {
        self.point
    }
}

#[derive(Clone, Debug)]
/// A Qeds data structure specialised to a 2d triangulation.
pub struct SurfaceTriangulation<T> {
    /// The quad-edge data structure we use as the basis for the triangulation.
    pub qeds: Qeds<Segment<T>, Space>,
    pub boundary_edge: EdgeTarget,
    pub bounds: Option<(Point, Point)>,
}

impl<T: Default + Clone> SurfaceTriangulation<T> {
    fn add_to_quad_unchecked(
        &mut self,
        mut edge_target: EdgeTarget,
        point: Point,
        data: T,
        retriangulate: bool
    ) -> EdgeTarget {
        unsafe {
            let first = self.qeds.edge_a_ref(edge_target).edge().point.point;
            let mut base = self
                .qeds
                .make_edge_with_a(Segment::new(first, data.clone()), Segment::new(point, data))
                .target();
            let return_value = base.sym();
            self.qeds.splice(base, edge_target);

            // TODO: not 100% sure this works
            let left_b = self.qeds.edge_a_ref(base).onext().rot().edge().point;
            let right_b = self.qeds.edge_a_ref(base).oprev().sym().rot().edge().point;

            self.qeds.edge_b_mut(base.sym().rot()).point = left_b;
            self.qeds.edge_b_mut(base.rot()).point = right_b;
            // {
            //     let aval = unsafe { self.qeds.edge_a_ref(base).onext().sym().rot().edge().point };
            //     unsafe {
            //         let ma = self.qeds.edge_b_mut(base.rot());
            //         ma.point = aval;
            //     }
            //     let bval = unsafe { self.qeds.edge_a_ref(base).oprev().sym().rot().edge().point };
            //     unsafe {
            //         let target = self.qeds.edge_a_ref(base).sym().rot().target();
            //         let mb = self.qeds.edge_b_mut(target);
            //         mb.point = bval;
            //     }
            // }
            let mut t = 0;
            loop {
                let base_ref = self.connect(edge_target, base.sym());
                edge_target = base_ref.oprev().target();
                base = base_ref.target();
                // TODO: the below BData settings need checking
                let left_b = base_ref.onext().rot().edge().point;
                let right_b = base_ref.oprev().sym().rot().edge().point;
                self.qeds.edge_b_mut(base.sym().rot()).point = left_b;
                self.qeds.edge_b_mut(base.rot()).point = right_b;
                if self.qeds.edge_a(edge_target.sym()).point.point == first {
                    break;
                }
            }
            edge_target = self.qeds.edge_a_ref(base).oprev().target();
            let e = edge_target;
            // We have to turn this off when doing bulk changes.
            if retriangulate {
                self.retriangulate_suspect_edges(e, point, first);
            }
            debug_assert_eq!(
                self.qeds.edge_a_ref(return_value).edge().point.point(),
                point
            );
            return_value
        }
    }
    fn connect<'a>(&'a mut self, a: EdgeTarget, b: EdgeTarget) -> EdgeRefA<'a, Segment<T>, Space> {
        let e = self.qeds.connect(a, b).target();
        // let aval = unsafe { self.qeds.edge_a_ref(e).sym().oprev().rot().edge().point };
        // unsafe {
        //     let ma = self.qeds.edge_b_mut(e.rot());
        //     ma.point = aval;
        // }
        // let bval = unsafe { self.qeds.edge_a_ref(e).oprev().rot().edge().point };
        // unsafe {
        //     let target = self.qeds.edge_a_ref(e).sym().rot().target();
        //     let mb = self.qeds.edge_b_mut(target);
        //     mb.point = bval;
        // }
        unsafe { self.qeds.edge_a_ref(e) }
    }
    /// This function accounts for the point lying on an existing edge or point.
    fn add_to_l_face(&mut self, edge_target: EdgeTarget, point: Point, data: T,retriangulate:bool) -> EdgeTarget {
        unsafe {
            let point_a = self.qeds.edge_a_ref(edge_target).edge().point.point;
            let point_b = self.qeds.edge_a_ref(edge_target).sym().edge().point.point;
            if point_a == point {
                // The point already lies on an existing point so just return
                // that edge.
                edge_target
            } else if point_b == point {
                // The point already lies on an existing point so just return
                // that edge.
                edge_target.sym()
            } else if self.qeds.edge_a_ref(edge_target).lies_right(point) == Lies::On {
                // The point lies on an edge, so we must invoke a procedure that
                // deletes and replaces that edge.
                self.add_point_to_edge_unchecked(edge_target, point, data,retriangulate)
            } else {
                self.add_to_quad_unchecked(edge_target, point, data, retriangulate)
            }
        }
    }
    /// The edge this returns should always have the added point at its origin.
    pub fn add_point(&mut self, mut point: Point, data: T,retriangulate:bool) -> Option<EdgeTarget> {
        point.snap();
        self.update_bounds(point);
        if let Some(edge_ref) = self.locate(point) {
            let edge_target = edge_ref.target();
            Some(self.add_to_l_face(edge_target, point, data,retriangulate))
        } else {
            // Point was out of bounds (probably)
            None
        }
    }

    /// Same as [`add_point_to_edge`] but does not check if the point is on one
    /// of the vertices of the edge. Will also restore constraints where
    /// necessary.
    fn add_point_to_edge_unchecked(
        &mut self,
        mut edge_target: EdgeTarget,
        point: Point,
        data: T,
        retriangulate:bool
    ) -> EdgeTarget {
        unsafe {
            {
                let oprev_ref = self.qeds.edge_a_ref(edge_target).oprev();

                let oprev = oprev_ref.target();
                drop(oprev_ref);

                self.qeds.delete(edge_target);
                edge_target = oprev;
            }
            let return_value = self.add_to_quad_unchecked(edge_target, point, data,retriangulate);
            return_value
        }
    }

    /// Add a point to a specified edge. If the point lies on one of the
    /// vertices just add it there.
    pub fn add_point_to_edge(
        &mut self,
        edge_target: EdgeTarget,
        point: Point,
        data: T,
        retriangulate:bool
    ) -> EdgeTarget {
        unsafe {
            let point_a = self.qeds.edge_a_ref(edge_target).edge().point.point;
            let point_b = self.qeds.edge_a_ref(edge_target).sym().edge().point.point;

            if point_a == point {
                edge_target
            } else if point_b == point {
                edge_target.sym()
            } else {
                self.add_point_to_edge_unchecked(edge_target, point, data,retriangulate)
            }
        }
    }
}

impl<T: Clone> SurfaceTriangulation<T> {
    /// Create a new surface from a series of points. The points are in 2d, so
    /// we don't need to worry about them all being co-planar at the moment.
    pub fn from_points(points: &[Point]) -> Self {
        todo!()
    }
    pub fn new(a: (Point, T), b: (Point, T), c: (Point, T)) -> Self {
        assert!(is_ccw(a.0, b.0, c.0));
        let mut qeds = Qeds::new();
        let edge_a = qeds
            .make_edge_with_ab(
                Segment::new(a.0, a.1.clone()),
                Segment::new(b.0, b.1.clone()),
                Space::Out,
                Space::In,
            )
            // .sym()
            .target();
        let edge_b = qeds
            .make_edge_with_ab(
                Segment::new(b.0, b.1),
                Segment::new(c.0, c.1.clone()),
                Space::Out,
                Space::In,
            )
            .target();
        // let edge_c = qeds
        //     .make_edge_with_a(Segment::new(c.0, c.1), Segment::new(a.0, a.1))
        //     .target();
        let result = unsafe {
            qeds.splice(edge_a.sym(), edge_b);
            // for (i, quad) in qeds.quads.iter() {
            //     println!("AQuad[{}]: {:?}", i, quad.edges_b);
            // }
            qeds.connect(edge_b, edge_a);
            qeds.quads[2].edges_b[0].point = Space::Out;
            // qeds.splice(edge_c.sym(), edge_a);
            SurfaceTriangulation {
                qeds,
                boundary_edge: edge_a,
                bounds: None,
            }
        };
        result
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

        let a_dest = self.qeds.edge_a_ref(a).sym().edge().point.clone();
        let b_dest = self.qeds.edge_a_ref(b).sym().edge().point.clone();
        self.qeds.edge_a_mut(e).point = a_dest;
        self.qeds.edge_a_mut(e.sym()).point = b_dest;
    }
    unsafe fn retriangulate_suspect_edges(
        &mut self,
        mut e: EdgeTarget,
        point: Point,
        first: Point,
    ) {
        // The suspect edges are e(.Onext.Lprev)^k for k=0,1,2,3...
        loop {
            let t = self.qeds.edge_a_ref(e).oprev().target();
            let t_dest = self.qeds.edge_a_ref(t).sym().edge().point.point;
            let e_dest = self.qeds.edge_a_ref(e).sym().edge().point.point;
            let e_org = self.qeds.edge_a_ref(e).edge().point.point;

            // TODO: we need to be cautious of infinite loops now that we're constrained.
            if self.qeds.edge_a_ref(e).lies_right_strict(t_dest)
                && del_test_ccw(e_org, t_dest, e_dest, point)
                && self.concave_test(e)
            {
                self.swap(e);
                // This is different from the algorithm in the paper
                e = self.qeds.edge_a_ref(e).oprev().target();
            } else if e_org == first {
                break;
            } else {
                e = self.qeds.edge_a_ref(e).onext().l_prev().target();
            }
        }
    }

    /// Warning: this is very inefficient and just for testing.
    fn retriangulate_all_single_pass(&mut self) -> usize {
        let mut swaps = 0;
        let edge_targets: Vec<EdgeTarget> =
            self.qeds.base_edges().map(|edge| edge.target()).collect();
        for e in edge_targets.into_iter() {
            unsafe {
                if self.del_test(e) && self.concave_test(e) {
                    swaps += 1;
                    self.swap(e);
                }
            }
        }
        swaps
    }

    /// Perform Delaunay swapping on the entire triangulation until complete.
    /// Should not be necessary, mainly included for testing.
    pub fn retriangulate_all(&mut self) -> usize {
        let mut iterations = 0;
        let mut total_swaps = 0;
        let total_swaps = loop {
            if iterations > 100 {
                panic!("too many triangulation iterations");
            }
            let swaps = self.retriangulate_all_single_pass();
            total_swaps += swaps;
            if swaps == 0 {
                break total_swaps;
            }
            iterations += 1;
        };
        // // TODO: Don't do this, this is a hack
        // unsafe {
        //     if let Some(out_target) = self.get_outside() {
        //         let start = out_target;
        //         let mut current = start;
        //         loop {
        //             let current_edge = self.qeds.edge_b_mut(current);
        //             println!("setting {:?} to {:?}", current_edge.point, Space::Out);
        //             current_edge.point = Space::Out;
        //             current = self.qeds.edge_b_ref(current).oprev().target;
        //             if current == start {
        //                 break;
        //             }
        //         }
        //     } else {
        //         panic!("no out")
        //     }
        // }

        total_swaps
    }

    fn get_outside(&self) -> Option<EdgeTarget> {
        unsafe {
            // Find at least one edge with an Out property.
            for (i, quad) in self.qeds.quads.iter() {
                if quad.edges_b[0].point == Space::Out {
                    return Some(EdgeTarget::new(i, 1, 0));
                }
                if quad.edges_b[1].point == Space::Out {
                    return Some(EdgeTarget::new(i, 3, 0));
                }
            }
        }
        None
    }

    // TODO: This has not yet been proved to be stable. It may also loop
    // inifintely, particularly with constrianed triangulations.
    pub fn locate(&self, point: Point) -> Option<EdgeRefA<Segment<T>, Space>> {
        use rand::Rng;
        let mut e = self.some_edge_a().unwrap();
        let mut rng = rand::thread_rng();
        let mut current_iterations = 0;
        let edge = loop {
            current_iterations += 1;
            if current_iterations > 2000 {
                for (a, b, c) in self.triangles() {
                    println!("Triangle: {:?}", (a.point, b.point, c.point));
                }
                panic!("locating failed for: {}", point);
            }
            // A random variable to determine if onext is tested first. If not
            // it is tested second.
            let onext_first: bool = rng.gen();
            if point == e.edge().point.point {
                break Some(e);
            } else if point == e.sym().edge().point.point {
                // If we lie on the eDest, return eSym(). This is a departure
                // from the paper which just returns e.
                break Some(e.sym());
            } else if e.lies_right_strict(point) {
                e = e.sym();
            } else {
                if onext_first {
                    if !e.onext().lies_right_strict(point) {
                        e = e.onext();
                    } else if !e.d_prev().lies_right_strict(point) {
                        e = e.d_prev();
                    } else {
                        break Some(e);
                    }
                } else {
                    if !e.d_prev().lies_right_strict(point) {
                        e = e.d_prev();
                    } else if !e.onext().lies_right_strict(point) {
                        e = e.onext();
                    } else {
                        break Some(e);
                    }
                }
            }
        };
        edge
    }

    /// Return the canonical tri within which this point is located.
    pub fn locate_tri(&self, point: Point) -> Option<EdgeRefA<Segment<T>, Space>> {
        self.locate(point).map(|edge| edge.get_tri_canonical())
    }
}

pub struct SegmentsIter<'a, T> {
    slab_iter: slab::Iter<'a, Quad<Segment<T>, Space>>,
    current_quad: Option<Quad<Segment<T>, Space>>,
}

impl<'a, T> SegmentsIter<'a, T> {
    pub fn new(triangulation: &'a SurfaceTriangulation<T>) -> Self {
        Self {
            slab_iter: triangulation.qeds.quads.iter(),
            current_quad: None,
        }
    }
}

impl<'a, T: Clone> Iterator for SegmentsIter<'a, T> {
    type Item = Edge<Segment<T>>;
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(current_quad) = &self.current_quad {
            let edge = current_quad.edges_a[1].clone();
            self.current_quad = None;
            Some(edge)
        } else {
            if let Some((_, quad)) = self.slab_iter.next() {
                let edge = quad.edges_a[0].clone();
                self.current_quad = Some(quad.clone());
                Some(edge)
            } else {
                None
            }
        }
    }
}

pub struct SegmentsIterMut<'a, T> {
    slab_iter: slab::IterMut<'a, Quad<Segment<T>, Space>>,
    next_segment: Option<&'a mut Segment<T>>,
}

impl<'a, T> SegmentsIterMut<'a, T> {
    pub fn new(triangulation: &'a mut SurfaceTriangulation<T>) -> Self {
        Self {
            slab_iter: triangulation.qeds.quads.iter_mut(),
            next_segment: None,
        }
    }
}

impl<'a, T: Clone> Iterator for SegmentsIterMut<'a, T> {
    type Item = &'a mut Segment<T>;
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(next_segment) = self.next_segment.take() {
            Some(next_segment)
        } else {
            if let Some((_, quad)) = self.slab_iter.next() {
                let (this, next) = quad.edges_a.split_at_mut(1);
                self.next_segment = Some(&mut next[0].point);
                Some(&mut this[0].point)
            } else {
                None
            }
        }
    }
}

impl<T> SurfaceTriangulation<T> {
    pub fn segments(&self) -> SegmentsIter<T> {
        SegmentsIter::new(self)
    }
    pub fn base_targets(&self) -> impl Iterator<Item=EdgeTarget> + '_ {
        self.qeds.quads.iter().map(|(i,_)|EdgeTarget::new(i,0,0))
    }
    pub fn segments_mut(&mut self) -> SegmentsIterMut<T> {
        SegmentsIterMut::new(self)
    }
    pub fn raw_triangles(&self) -> RawTriangleIter<T> {
        RawTriangleIter::new(self)
    }
    pub fn triangles(&self) -> TriangleIter<T> {
        TriangleIter::new(self)
    }
    pub fn triangle_edges(&self) -> TriangleEdgeIter<T> {
        TriangleEdgeIter::new(self)
    }
    // pub fn triangles_mut(&mut self) -> TriangleIterMut<T> {
    //     TriangleIterMut::new(self)
    // }
    pub fn nodes(&self) -> NodeIter<T> {
        NodeIter::new(self)
    }

    pub fn qeds(&self) -> Option<&Qeds<Segment<T>, Space>> {
        Some(&self.qeds)
    }

    pub fn some_edge_a(&self) -> Option<EdgeRefA<Segment<T>, Space>> {
        let (i, _) = self.qeds.quads.iter().next()?;
        unsafe { Some(self.qeds.edge_a_ref(EdgeTarget::new(i, 0, 0))) }
    }

    fn update_bounds(&mut self, point: Point) {
        match self.bounds {
            None => self.bounds = Some((point, point)),
            Some(ref mut bounds) => {
                if point.x < bounds.0.x {
                    bounds.0.x = point.x
                }
                if point.y < bounds.0.y {
                    bounds.0.y = point.y
                }
                if point.x > bounds.1.x {
                    bounds.1.x = point.x
                }
                if point.y > bounds.1.y {
                    bounds.1.y = point.y
                }
            }
        }
    }

    /// TODO: remove the recursion from this function as we don't have tail
    /// recursions. Edges are returned with each edge pointing from right to
    /// left across the intersection line. Second return value is an edge from
    /// the final point. It is necessary that the two points are already in the
    /// trianglulation as it is required for the end condition. If we pass
    /// through a point on a constraint segment, we also report that as an
    /// intersection. TODO: this needs extensive revision.
    fn find_intersections_between_points(&self, a: Point, b: Point) -> IntersectionIter<T> {
        IntersectionIter::new(self, a, b)
    }

    pub fn boundary(&self) -> BoundaryIter<Segment<T>, Space> {
        BoundaryIter::new(&self.qeds, self.boundary_edge)
    }

    /// Determine whether an Edge (i.e. two NavTris) satisfies the Delaunay
    /// criterion. An edge with only a single NavTri will return True.
    unsafe fn del_test(&self, e: EdgeTarget) -> bool {
        // Get the edge.
        let edge = self.qeds.edge_a_ref(e);
        // Get all of the vertices around this egdge in a CCW order.
        let a = edge.oprev().edge().point.point;
        let b = edge.r_prev().sym().edge().point.point;
        let c = edge.l_next().edge().point.point;
        let d = edge.onext().sym().edge().point.point;
        del_test_ccw(a, b, c, d)
    }

    /// Test whether an edge is the diagonal of a concave quad.
    unsafe fn concave_test(&self, e: EdgeTarget) -> bool {
        // Get the edge.
        let edge = self.qeds.edge_a_ref(e);
        // Get all of the vertices around this egdge in a CCW order.
        let a = edge.oprev().edge().point.point;
        let b = edge.r_prev().sym().edge().point.point;
        let c = edge.l_next().edge().point.point;
        let d = edge.onext().sym().edge().point.point;
        if left_or_right(a, b, c) == Direction::Right {
            return false;
        }
        if left_or_right(b, c, d) == Direction::Right {
            return false;
        }
        if left_or_right(c, d, a) == Direction::Right {
            return false;
        }
        if left_or_right(d, a, b) == Direction::Right {
            return false;
        }
        true
    }

    /// If at any point the boundary loops back on itself, we know it to be
    /// linear. Technically it only shows that one edge is, but one of our
    /// invariant is that only one linear edge is not possible.
    pub fn is_linear(&self) -> bool {
        for edge in self.boundary() {
            if edge.sym().onext() == edge.sym() {
                return true;
            }
        }
        false
    }
}

#[derive(Clone, Debug)]
pub struct IntersectionIter<'a, T> {
    triangulation: &'a SurfaceTriangulation<T>,
    current_intersection: Option<Intersection>,
    a: Point,
    b: Point,
}

impl<'a, T> IntersectionIter<'a, T> {
    pub fn new(triangulation: &'a SurfaceTriangulation<T>, a: Point, b: Point) -> Self {
        Self {
            triangulation,
            current_intersection: None,
            a,
            b,
        }
    }
}

impl<'a, T: Clone> Iterator for IntersectionIter<'a, T> {
    type Item = Intersection;
    fn next(&mut self) -> Option<Self::Item> {
        use Direction::*;
        // TODO: this belongs in the base trianglulation.

        if let Some(current_intersection) = self.current_intersection {
            // We have now taken care of the beginning and can just iterate
            // through triangles.

            // We can unwrap here as we know that the vector as we have pushed a
            // value but don't pop.
            match current_intersection {
                Intersection::Point(intersecting_edge) => {
                    // println!(
                    //     "Current Intersection: {}",
                    //     intersecting_edge.edge(self.triangulation).point.point()
                    // );
                    if edge_from_target(intersecting_edge, self.triangulation)
                        .point
                        .point()
                        == self.b
                    {
                        return None;
                    }
                    // This is special case, where the point [`a`] lies on the
                    // origin point of [`start`]. First we need to place
                    // ourselves on the correct triangle. To do that we iterate
                    // around the point until we are sure that the line passes
                    // through a given triangle. This is indicated by subsequent
                    // "spokes" switching from left to right. If the line passes
                    // through a vertex we can deal with it here.
                    let initial_start =
                        unsafe { self.triangulation.qeds.edge_a_ref(intersecting_edge) };
                    let mut spoke = initial_start;
                    // First we loop around until we are on a spoke where b is
                    // to the left.
                    loop {
                        let start_dir = left_or_right(
                            spoke.edge().point.point(),
                            spoke.sym().edge().point.point(),
                            self.b,
                        );
                        match start_dir {
                            Direction::Left => {
                                // We have found a spoke where b is to the left.
                                // Now we need to progress around the origin
                                // until we have a left/right pair.
                                loop {
                                    let next_edge = spoke.onext();
                                    let next_dir = left_or_right(
                                        next_edge.edge().point.point(),
                                        next_edge.sym().edge().point.point(),
                                        self.b,
                                    );
                                    match next_dir {
                                        Straight => {
                                            // We have a straight line, so we add that point
                                            // (on the other end) as an intersection point.
                                            let intersection =
                                                Intersection::Point(next_edge.sym().target());
                                            self.current_intersection = Some(intersection);
                                            return Some(intersection);
                                        }
                                        Left => {
                                            spoke = next_edge;
                                        }
                                        Right => break,
                                    }
                                }
                                break;
                            }
                            // TODO: this does not account for the 180 case
                            Direction::Straight => {
                                if self.b.distance(spoke.sym().edge().point.point())
                                    <= self.b.distance(spoke.edge().point.point())
                                {
                                    // We have a straight line, so we add that point
                                    // (on the other end) as an intersection point.
                                    let intersection = Intersection::Point(spoke.sym().target());
                                    self.current_intersection = Some(intersection);
                                    return Some(intersection);
                                } else {
                                    spoke = spoke.onext();
                                }
                            }
                            Direction::Right => {
                                spoke = spoke.onext();
                                if spoke == initial_start {
                                    panic!("looped around ring");
                                }
                            }
                        }
                    }
                    // println!(
                    //     "RightSpoke: {}-{}",
                    //     spoke.edge().point.point(),
                    //     spoke.sym().edge().point.point()
                    // );
                    // Point a is on startOrg. This changes our initial question
                    // somewhat. We know that the line either intersects the edge
                    // opposite a, or also passes through one of the other vertices.
                    let opposite_edge = spoke.l_next();
                    // println!(
                    //     "OppositeEdge: {}-{}",
                    //     opposite_edge.edge().point.point(),
                    //     opposite_edge.sym().edge().point.point()
                    // );
                    // Does the line pass through the vertex to the right?
                    let right = opposite_edge;
                    let right_point = right.edge().point.point();
                    let passes_through_right =
                        left_or_right(self.a, right_point, self.b) == Direction::Straight;
                    if passes_through_right {
                        // If so we need to start the calculation from that vertex.
                        let intersection = Intersection::Point(right.target());
                        self.current_intersection = Some(intersection);
                        // println!("Passes through right");
                        return Some(intersection);
                    }
                    // Does the line pass through the vertex to the left?
                    let left = opposite_edge.sym();
                    let left_point = left.edge().point.point();
                    let passes_through_left =
                        left_or_right(self.a, left_point, self.b) == Direction::Straight;
                    if passes_through_left {
                        // If so we need to start the calculation from that vertex.
                        let intersection = Intersection::Point(left.target());
                        self.current_intersection = Some(intersection);
                        // println!("Passes through left");
                        return Some(intersection);
                    }
                    // If it passes through neither, then the intersecting edges
                    // is simply the opposite edge, BUT, if pb is to the left of
                    // the opposite edge, there is no intersection.
                    if opposite_edge.lies_left(self.b) {
                        // println!("Past End");
                        return None;
                    }
                    let intersection = Intersection::Edge(opposite_edge.target());
                    self.current_intersection = Some(intersection);
                    // println!(
                    //     "NewCurrentIntersection: {}-{}",
                    //     opposite_edge.edge().point.point(),
                    //     opposite_edge.sym().edge().point.point()
                    // );
                    return Some(intersection);
                }
                Intersection::Edge(intersecting_edge) => {
                    // println!(
                    //     "Current Intersection: {}-{}",
                    //     intersecting_edge.edge(self.triangulation).point.point(),
                    //     intersecting_edge
                    //         .sym()
                    //         .edge(self.triangulation)
                    //         .point
                    //         .point()
                    // );
                    let intersecting_edge =
                        unsafe { self.triangulation.qeds.edge_a_ref(intersecting_edge) };
                    // Now we iterate through all of the triangles.
                    // If the final point is to the left of the intersecting edge, it is
                    // within the triangle, so we don't include the intersection and can
                    // break.
                    match left_or_right(
                        intersecting_edge.edge().point.point(),
                        intersecting_edge.sym().edge().point.point(),
                        self.b,
                    ) {
                        Direction::Left => {
                            return None;
                        }
                        _ => (),
                    }
                    // intersections.push(Intersection::Edge(intersecting_edge.target()));
                    // To get the same edge on the other triangle, we simply need the
                    // edge in the opposite direction.
                    let e = intersecting_edge.sym();
                    // Get the point opposite.
                    let c_edge = e.l_next().sym();
                    let c = c_edge.edge().point.point();
                    match left_or_right(self.a, c, self.b) {
                        Direction::Left => {
                            let intersection = Intersection::Edge(e.l_next().l_next().target());
                            self.current_intersection = Some(intersection);
                            return Some(intersection);
                        }
                        Direction::Straight => {
                            let intersection = Intersection::Point(c_edge.target());
                            self.current_intersection = Some(intersection);
                            return Some(intersection);
                        }
                        Direction::Right => {
                            let intersection = Intersection::Edge(e.l_next().target());
                            self.current_intersection = Some(intersection);
                            return Some(intersection);
                        }
                    }
                }
            }
        } else {
            // If we are not already starting at an intersection, we need to
            // find the first one. This is the "start" of the iter.
            // [`start`] is an edge of the triangle in which the fist point ([`a`])
            // is located OR [`a`] lies on this edge.
            let start = self.triangulation.locate(self.a).unwrap();
            // End condition. If b also lies on this edge we are done OR if the
            // two points are the same we are done.
            if start.edge().point.point() == self.b || self.a == self.b {
                return None;
            }
            if self.a == start.edge().point.point() {
                // Special case.
                self.current_intersection = Some(Intersection::Point(start.target()));
                return self.next();
            }

            // TODO: not quite correct either, doesn't accout for starting on an
            // edge.
            if start.in_left_face(self.b) {
                return None;
            }

            // Given a ray originating at [`a`]
            // and passing through [`b`], this ray must pass through one of the
            // edges of the triangle (or originate on one of the edges, which
            // counts as an intersection).

            // We must test agains three lines made by combining [`a`] with each
            // of the vertices of the triangle. These lines are ax,ay, and az in
            // CCW order.
            let x_edge = start;
            let x = x_edge.edge().point.point();
            let y_edge = start.l_next();
            let y = y_edge.edge().point.point();
            let z_edge = start.l_next().l_next();
            let z = z_edge.edge().point.point();
            let ax = (self.a, x);
            let ay = (self.a, y);
            let az = (self.a, z);

            let right_of_ax = left_or_right(ax.0, ax.1, self.b);
            let right_of_ay = left_or_right(ay.0, ay.1, self.b);
            let right_of_az = left_or_right(az.0, az.1, self.b);

            // Find the point at which the direction changes from left to right.
            let intersection = match (right_of_ax, right_of_ay, right_of_az) {
                (Direction::Straight, _, _) => Intersection::Point(x_edge.target()),
                (_, Direction::Straight, _) => Intersection::Point(y_edge.target()),
                (_, _, Direction::Straight) => Intersection::Point(z_edge.target()),
                (Left, Right, _) => Intersection::Edge(x_edge.target()),
                (_, Left, Right) => Intersection::Edge(y_edge.target()),
                (_, Right, Left) => Intersection::Edge(z_edge.target()),
                _ => unreachable!(),
            };
            self.current_intersection = Some(intersection);
            Some(intersection)
        }
    }
}

impl<'a, T> EdgeRefA<'a, Segment<T>, Space> {
    pub fn triangle_across(&self) -> Self {
        let edge_on_next_tri = self.sym();
        let canonical_edge = edge_on_next_tri.get_tri_canonical();
        canonical_edge
    }
}

#[derive(Clone, Copy)]
pub struct RawTriangleIter<'a, T> {
    triangulation: &'a SurfaceTriangulation<T>,
    next: EdgeTarget,
}

impl<'a, T> RawTriangleIter<'a, T> {
    pub fn new(triangulation: &'a SurfaceTriangulation<T>) -> Self {
        Self {
            triangulation,
            next: EdgeTarget::new(0, 0, 0),
        }
    }

    fn inc(&mut self) {
        if self.next.r == 0 {
            self.next.r = 2;
        } else {
            self.next = EdgeTarget::new(self.next.e + 1, 0, 0);
        }
    }
}

impl<'a, T> Iterator for RawTriangleIter<'a, T> {
    type Item = EdgeRefA<'a, Segment<T>, Space>;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // First we check that the edge actually exists for this given
            // EdgeTarget.
            if self.triangulation.qeds.quads.contains(self.next.e) {
                let edge_ref = unsafe { self.triangulation.qeds.edge_a_ref(self.next) };
                self.inc();
                if edge_ref.is_tri_canonical() && edge_ref.is_tri_real() {
                    break Some(edge_ref);
                }
            } else {
                self.inc();
                if self.next.e >= self.triangulation.qeds.quads.len() {
                    break None;
                }
            }
        }
    }
}

// #[derive(Clone)]
// pub struct TriangleIterMut<'a, T> {
//     triangulation: &'a mut SurfaceTriangulation<T>,
//     next: EdgeTarget,
// }

// impl<'a, T> TriangleIterMut<'a, T> {
//     pub fn new(triangulation: &'a mut SurfaceTriangulation<T>) -> Self {
//         Self {
//             triangulation,
//             // We skip zero because that is the edge for the infinite triangle. (TODO: currently it is e0r2f0);
//             next: EdgeTarget::new(0, 2, 0),
//         }
//     }

//     fn inc(&mut self) {
//         if self.next.r == 0 {
//             self.next.r = 2;
//         } else {
//             self.next = EdgeTarget::new(self.next.e + 1, 0, 0);
//         }
//     }
// }

// impl<'a, T> Iterator for TriangleIterMut<'a, T> {
//     type Item = (&'a mut Segment<T>, &'a mut Segment<T>, &'a mut Segment<T>);
//     fn next(&mut self) -> Option<Self::Item> {
//         loop {
//             // First we check that the edge actually exists for this given
//             // EdgeTarget.
//             if self.triangulation.qeds.quads.contains(self.next.e) {
//                 // self.triangulation.qeds.
//                 let tri_a: &mut Edge<Segment<T>> = unsafe {self.triangulation.qeds.edge_a_mut(self.next)};
//                 // let tri_a = unsafe { self.triangulation.qeds.edge_a_ref(self.next).edge_mut() };
//                 let tri_b = unsafe { self.triangulation.qeds.edge_a_ref(self.next).l_next().edge() };
//                 let tri_c = unsafe { self.triangulation.qeds.edge_a_ref(self.next).l_next().l_next().edge_mut() };
//                 self.inc();
//                 let is_tri_canonical: bool = unsafe { self.triangulation.qeds.edge_a_ref(self.next).is_tri_canonical() };
//                 if is_tri_canonical {
//                     break Some((&mut tri_a.point,&mut tri_b.point,&mut tri_c.point));
//                 }
//             } else {
//                 self.inc();
//                 if self.next.e >= self.triangulation.qeds.quads.len() {
//                     break None;
//                 }
//             }
//         }
//     }
// }

#[derive(Clone, Copy)]
pub struct TriangleIter<'a, T> {
    raw_triangles: RawTriangleIter<'a, T>,
}

impl<'a, T> TriangleIter<'a, T> {
    pub fn new(triangulation: &'a SurfaceTriangulation<T>) -> Self {
        Self {
            raw_triangles: triangulation.raw_triangles(),
        }
    }
}

impl<'a, T: Clone> Iterator for TriangleIter<'a, T> {
    type Item = (Segment<T>, Segment<T>, Segment<T>);
    fn next(&mut self) -> Option<Self::Item> {
        let tri_a = self.raw_triangles.next()?;
        let tri_b = tri_a.l_next();
        let tri_c = tri_b.l_next();
        let ccw = is_ccw(
            tri_a.edge().point.point,
            tri_b.edge().point.point,
            tri_c.edge().point.point,
        );
        if !ccw {
            println!("a: {}",tri_a.edge().point.point);
            println!("b: {}",tri_b.edge().point.point);
            println!("c: {}",tri_c.edge().point.point);
        }
        assert!(ccw);
        Some((
            tri_a.edge().point.clone(),
            tri_b.edge().point.clone(),
            tri_c.edge().point.clone(),
        ))
    }
}

#[derive(Clone, Copy)]
pub struct TriangleEdgeIter<'a, T> {
    raw_triangles: RawTriangleIter<'a, T>,
}

impl<'a, T> TriangleEdgeIter<'a, T> {
    pub fn new(triangulation: &'a SurfaceTriangulation<T>) -> Self {
        Self {
            raw_triangles: triangulation.raw_triangles(),
        }
    }
}

impl<'a, T: Clone> Iterator for TriangleEdgeIter<'a, T> {
    type Item = (
        (EdgeTarget, Segment<T>),
        (EdgeTarget, Segment<T>),
        (EdgeTarget, Segment<T>),
    );
    fn next(&mut self) -> Option<Self::Item> {
        let tri_a = self.raw_triangles.next()?;
        let tri_b = tri_a.l_next();
        let tri_c = tri_b.l_next();
        assert_eq!(tri_c.l_next().target, tri_a.target);
        Some((
            (tri_a.target, tri_a.edge().point.clone()),
            (tri_b.target, tri_b.edge().point.clone()),
            (tri_c.target, tri_c.edge().point.clone()),
        ))
    }
}

// TODO: does this actuall return all the nodes?
#[derive(Clone, Copy)]
pub struct NodeIter<'a, T> {
    triangulation: &'a SurfaceTriangulation<T>,
    next: EdgeTarget,
}

impl<'a, T> NodeIter<'a, T> {
    pub fn new(triangulation: &'a SurfaceTriangulation<T>) -> Self {
        Self {
            triangulation,
            // We skip zero because that is the edge for the infinite triangle. (TODO: currently it is e0r2f0);
            next: EdgeTarget::new(0, 2, 0),
        }
    }

    fn inc(&mut self) {
        if self.next.r == 0 {
            self.next.r = 2;
        } else {
            self.next = EdgeTarget::new(self.next.e + 1, 0, 0);
        }
    }
}

impl<'a, T> Iterator for NodeIter<'a, T> {
    type Item = NodeTarget;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // First we check that the edge actually exists for this given
            // EdgeTarget.
            if self.triangulation.qeds.quads.contains(self.next.e) {
                let edge_ref = unsafe { self.triangulation.qeds.edge_a_ref(self.next) };
                self.inc();
                let is_tri_canonical: bool = edge_ref.is_tri_canonical();
                if is_tri_canonical {
                    break Some(NodeTarget::new_unchecked(edge_ref.target()));
                }
            } else {
                self.inc();
                if self.next.e >= self.triangulation.qeds.quads.len() {
                    break None;
                }
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Space {
    In,
    Out,
}

impl Default for Space {
    fn default() -> Self {
        Self::In
    }
}

impl<'a, T> EdgeRefA<'a, Segment<T>, Space> {
    /// An edge is the canonical edge for a tri iff it has the lowest e value
    /// for the tri.
    fn is_tri_canonical(&self) -> bool {
        let current_e = self.target().e;
        let next_edge = self.sym_rot().onext();
        let next_e = next_edge.target().e;
        let second_e = next_edge.onext().target().e;
        (current_e < next_e) && (current_e < second_e)
    }

    fn is_tri_real(&self) -> bool {
        let first = self.target().sym().rot();
        let mut current = first;
        loop {
            // TODO: this is very inefficient, fix BData
            let edge = unsafe { self.qeds.edge_b_ref(current) };
            if edge.edge().point == Space::Out {
                return false;
            }
            current = edge.oprev().target();
            if current == first {
                break;
            }
        }
        // for edge in self.l_face().edges {
        // }
        true
    }

    pub fn get_tri_canonical(&self) -> Self {
        let second_edge = self.l_next();
        let third_edge = second_edge.l_next();
        let mut canonical_edge = *self;
        if second_edge.target().e < canonical_edge.target().e {
            canonical_edge = second_edge;
        }
        if third_edge.target().e < canonical_edge.target().e {
            canonical_edge = third_edge;
        }
        canonical_edge
    }
}

// impl EdgeTarget {
//     /// An edge is the canonical edge for a tri iff it has the lowest e value
//     /// for the tri.
//     fn is_tri_canonical(self, triangulation: &SurfaceTriangulation) -> bool {
//         let edge_ref = unsafe { triangulation.qeds.edge_a_ref(self) };
//         edge_ref.is_tri_canonical()
//     }

//     fn is_tri_real(self, triangulation: &SurfaceTriangulation) -> bool {
//         let edge_ref = unsafe { triangulation.qeds.edge_a_ref(self) };
//         edge_ref.is_tri_real()
//     }

//     pub fn get_tri_canonical(self, triangulation: &SurfaceTriangulation) -> NodeTarget {
//         let edge_ref = unsafe { triangulation.qeds.edge_a_ref(self) };
//         let canonical_edge = edge_ref.get_tri_canonical();
//         NodeTarget(canonical_edge.target())
//     }

//     pub fn triangle_across(self, triangulation: &SurfaceTriangulation) -> NodeTarget {
//         let edge_ref = unsafe { triangulation.qeds.edge_a_ref(self) };
//         let triangle_across = edge_ref.triangle_across();
//         NodeTarget(triangle_across.target())
//     }

//     pub fn edge(self, triangulation: &SurfaceTriangulation) -> Edge<Segment> {
//         let edge_ref = unsafe { triangulation.qeds.edge_a_ref(self) };
//         edge_ref.edge().clone()
//     }
// }

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
    // TODO: hard-code this
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
    let det = matrix.determinant();
    // TODO: It seems the LU algorithm from nalgebra is better
    // let new_det = determinant_4x4(a.x,
    //     a.y,
    //     a.x.powi(2) + a.y.powi(2),
    //     1.0,
    //     b.x,
    //     b.y,
    //     b.x.powi(2) + b.y.powi(2),
    //     1.0,
    //     c.x,
    //     c.y,
    //     c.x.powi(2) + c.y.powi(2),
    //     1.0,
    //     d.x,
    //     d.y,
    //     d.x.powi(2) + d.y.powi(2),
    //     1.0,);
    // assert_eq!(det, new_det);
    det > 0.0
}
fn determinant_4x4(
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    f: f64,
    g: f64,
    h: f64,
    i: f64,
    j: f64,
    k: f64,
    l: f64,
    m: f64,
    n: f64,
    o: f64,
    p: f64,
) -> f64 {
    a * determinant_3x3(f, g, h, j, k, l, n, o, p) - b * determinant_3x3(e, g, h, i, k, l, m, o, p)
        + c * determinant_3x3(e, f, h, i, j, l, m, n, p)
        - d * determinant_3x3(e, f, g, i, j, k, m, n, o)
}

fn determinant_3x3(a: f64, b: f64, c: f64, d: f64, e: f64, f: f64, g: f64, h: f64, i: f64) -> f64 {
    a * determinant_2x2(e, f, h, i) - b * determinant_2x2(d, f, g, i)
        + c * determinant_2x2(d, e, g, h)
}

pub fn is_ccw(p1: Point, p2: Point, p3: Point) -> bool {
    determinant_3x3(p1.x, p1.y, 1.0, p2.x, p2.y, 1.0, p3.x, p3.y, 1.0) > 0.0
}

fn determinant_2x2(a: f64, b: f64, c: f64, d: f64) -> f64 {
    a * d - b * c
}

pub fn has_edge<T>(triangulation: &SurfaceTriangulation<T>, pa: Point, pb: Point) -> bool {
    get_edge(triangulation, pa, pb).is_some()
}

pub fn get_edge<T>(
    triangulation: &SurfaceTriangulation<T>,
    pa: Point,
    pb: Point,
) -> Option<EdgeRefA<'_, Segment<T>, Space>> {
    for (i, quad) in triangulation.qeds().unwrap().quads.iter() {
        let edge1 = quad.edges_a[0].point.point;
        let edge2 = quad.edges_a[1].point.point;
        let edge_target = if edge1 == pa && edge2 == pb {
            EdgeTarget::new(i, 0, 0)
        } else if edge1 == pb && edge2 == pa {
            EdgeTarget::new(i, 2, 0)
        } else {
            continue;
        };
        return Some(unsafe { triangulation.qeds().unwrap().edge_a_ref(edge_target) });
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn surface_one_point_triangulation_only() {
        let v1 = Point::new(4.19, 0.0);
        let v2 = Point::new(4.189999, -15.0);
        let v3 = Point::new(17.99, -0.000002);
        let mut triangulation =
            SurfaceTriangulation::new((v1, (0.0, 0.0)), (v2, (0.0, 0.0)), (v3, (0.0, 0.0)));
        assert_eq!(triangulation.triangles().count(), 1);
        debug_assert_spaces(&triangulation);
        {
            let centroids: Vec<_> = triangulation
                .triangles()
                .map(|(a, b, c)| {
                    println!("a {} b {} c {}", a.point, b.point, c.point);
                    let ax = (a.point.x + b.point.x + c.point.x) / 3.0;
                    let ay = (a.point.y + b.point.y + c.point.y) / 3.0;
                    Point::new(ax, ay)
                })
                .collect();
            for centroid in centroids {
                println!("adding point: {}", centroid);
                for (i, quad) in triangulation.qeds.quads.iter() {
                    println!(
                        "BeforeQuad[{}]: {}-{} {:?}-{:?}",
                        i,
                        quad.edges_a[0].point.point,
                        quad.edges_a[1].point.point,
                        quad.edges_b[0].point,
                        quad.edges_b[1].point
                    );
                }
                triangulation.add_point(centroid, (0.0, 0.0), true);
                for (i, quad) in triangulation.qeds.quads.iter() {
                    println!(
                        "AfterQuad[{}]: {}-{} {:?}-{:?}",
                        i,
                        quad.edges_a[0].point.point,
                        quad.edges_a[1].point.point,
                        quad.edges_b[0].point,
                        quad.edges_b[1].point
                    );
                }
            }
        }
        assert_eq!(triangulation.qeds.base_edges().count(), 6);
        // let p1 = Point::new(4.0, 1.0);
        // triangulation.add_point(p1, (0.0, 0.0));
        debug_assert_spaces(&triangulation);
        assert_eq!(triangulation.triangles().count(), 3);
    }

    #[test]
    fn surface_split_edge_triangulation_only() {
        let v1 = Point::new(0.0, 0.0);
        let v2 = Point::new(5.0, 0.0);
        let v3 = Point::new(5.0, 5.0);
        let mut triangulation =
            SurfaceTriangulation::new((v1, (0.0, 0.0)), (v2, (0.0, 0.0)), (v3, (0.0, 0.0)));
        assert_eq!(triangulation.triangles().count(), 1);
        debug_assert_spaces(&triangulation);
        {
            let mut line_midpoints = Vec::with_capacity(triangulation.qeds.quads.len() * 3);
            for (a, b, c) in
                triangulation
                    .triangle_edges()
                    .map(|((ta, edge_a), (tb, edge_b), (tc, edge_c))| {
                        let ax = (edge_a.point.x + edge_b.point.x) / 2.0;
                        let ay = (edge_a.point.y + edge_b.point.y) / 2.0;

                        let bx = (edge_b.point.x + edge_c.point.x) / 2.0;
                        let by = (edge_b.point.y + edge_c.point.y) / 2.0;

                        let cx = (edge_c.point.x + edge_a.point.x) / 2.0;
                        let cy = (edge_c.point.y + edge_a.point.y) / 2.0;
                        let a = Point::new(ax, ay);
                        let b = Point::new(bx, by);
                        let c = Point::new(cx, cy);
                        ((ta, a), (tb, b), (tc, c))
                    })
            {
                line_midpoints.push(a);
            }
            let triangles: Vec<_> = triangulation
                .triangles()
                .map(|(a, b, c)| (a.point, b.point, c.point))
                .collect();
            println!("###Triangles Before");
            for t in triangles {
                println!("Triangle: {}-{}-{}", t.0, t.1, t.2);
            }
            for (i, quad) in triangulation.qeds.quads.iter() {
                println!("BeforeQuad[{}]: {:?}", i, quad.edges_b);
            }
            let (edge_target, midpoint) = line_midpoints.first().unwrap();
            let edge_target = *edge_target;
            let point = *midpoint;
            let n1 = triangulation.triangle_edges().count();
            // let e = triangulation.add_point_to_edge(edge_target, midpoint, (0.0, 0.0));
            unsafe {
                let oprev = triangulation.qeds.edge_a_ref(edge_target).oprev().target();
                // drop(oprev_ref);

                triangulation.qeds.delete(edge_target);
                let mut edge_target = oprev;
                let first = triangulation
                    .qeds
                    .edge_a_ref(edge_target)
                    .edge()
                    .point
                    .point;
                let mut base = triangulation
                    .qeds
                    .make_edge_with_a(
                        Segment::new(first, (0.0, 0.0)),
                        Segment::new(point, (0.0, 0.0)),
                    )
                    .target();
                let return_value = base.sym();
                triangulation.qeds.splice(base, edge_target);

                // TODO: not 100% sure this works
                let left_b = triangulation
                    .qeds
                    .edge_a_ref(base)
                    .onext()
                    .rot()
                    .edge()
                    .point;
                let right_b = triangulation
                    .qeds
                    .edge_a_ref(base)
                    .oprev()
                    .sym()
                    .rot()
                    .edge()
                    .point;

                triangulation.qeds.edge_b_mut(base.sym().rot()).point = left_b;
                triangulation.qeds.edge_b_mut(base.rot()).point = right_b;
                loop {
                    let base_ref = triangulation.connect(edge_target, base.sym());
                    edge_target = base_ref.oprev().target();
                    base = base_ref.target();
                    // TODO: the below BData settings need checking
                    let left_b = base_ref.onext().rot().edge().point;
                    let right_b = base_ref.oprev().sym().rot().edge().point;
                    triangulation.qeds.edge_b_mut(base.sym().rot()).point = left_b;
                    triangulation.qeds.edge_b_mut(base.rot()).point = right_b;
                    if triangulation.qeds.edge_a(edge_target.sym()).point.point == first {
                        break;
                    }
                }
                edge_target = triangulation.qeds.edge_a_ref(base).oprev().target();
                let e = edge_target;
                triangulation.retriangulate_suspect_edges(e, point, first);
                debug_assert_eq!(
                    triangulation
                        .qeds
                        .edge_a_ref(return_value)
                        .edge()
                        .point
                        .point(),
                    point
                );
            }

            for (i, quad) in triangulation.qeds.quads.iter() {
                println!("AfterQuad[{}]: {:?}", i, quad.edges_b);
            }
            let n2 = triangulation.triangle_edges().count();
            let triangles: Vec<_> = triangulation
                .triangles()
                .map(|(a, b, c)| (a.point, b.point, c.point))
                .collect();
            println!("###Triangles After");
            for t in triangles {
                println!("Triangle: {}-{}-{}", t.0, t.1, t.2);
            }
            println!("n1: {} n2: {}", n1, n2);
        }
        assert_eq!(triangulation.qeds.base_edges().count(), 6);
        // let p1 = Point::new(4.0, 1.0);
        // triangulation.add_point(p1, (0.0, 0.0));
        debug_assert_spaces(&triangulation);
        assert_eq!(triangulation.triangles().count(), 2);
    }

    //     #[test]
    //     fn empty_triangulation() {
    //         let triangulation = SurfaceTriangulation::new();
    //         let tris: Vec<EdgeRefA<Segment, Space>> = triangulation.triangles().collect();
    //         assert_eq!(tris.len(), 2);
    //         let tri_info = triangulation.classify_triangles();
    //         for (target, info) in tri_info.0 {
    //             println!("EdgeTarget: {:?} Info: {:?}", target, info);
    //         }
    //     }

    //     #[test]
    //     fn one_point_triangulation_location() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         // assert_eq!(triangulation.qeds().unwrap().quads.len(), 5);
    //         let south = unsafe { triangulation.qeds.edge_a_ref(EdgeTarget::new(0, 2, 0)) };
    //         assert_eq!(
    //             south.edge().point.point,
    //             Point::new(SafeFloat::MIN.0, SafeFloat::MIN.0)
    //         );
    //         assert_eq!(
    //             south.sym().edge().point.point,
    //             Point::new(SafeFloat::MAX.0, SafeFloat::MIN.0)
    //         );

    //         assert_eq!(south.target(), south.l_next().l_next().l_next().target());
    //         assert_eq!(
    //             south.l_next().l_next().oprev().target(),
    //             EdgeTarget::new(2, 0, 0)
    //         );
    //         assert_eq!(
    //             south.l_next().l_next().oprev().l_next().target(),
    //             EdgeTarget::new(4, 0, 0)
    //         );

    //         let east = unsafe { triangulation.qeds.edge_a_ref(EdgeTarget::new(1, 0, 0)) };
    //         assert_eq!(
    //             east.edge().point.point,
    //             Point::new(SafeFloat::MAX.0, SafeFloat::MIN.0)
    //         );
    //         assert_eq!(
    //             east.sym().edge().point.point,
    //             Point::new(SafeFloat::MAX.0, SafeFloat::MAX.0)
    //         );

    //         let centre = unsafe { triangulation.qeds.edge_a_ref(EdgeTarget::new(3, 0, 0)) };
    //         assert_eq!(
    //             centre.edge().point.point,
    //             Point::new(SafeFloat::MAX.0, SafeFloat::MAX.0)
    //         );
    //         assert_eq!(
    //             centre.sym().edge().point.point,
    //             Point::new(SafeFloat::MIN.0, SafeFloat::MIN.0)
    //         );

    //         let north = unsafe { triangulation.qeds.edge_a_ref(EdgeTarget::new(2, 0, 0)) };
    //         assert_eq!(
    //             north.edge().point.point,
    //             Point::new(SafeFloat::MAX.0, SafeFloat::MAX.0)
    //         );
    //         assert_eq!(
    //             north.sym().edge().point.point,
    //             Point::new(SafeFloat::MIN.0, SafeFloat::MAX.0)
    //         );

    //         let west = unsafe { triangulation.qeds.edge_a_ref(EdgeTarget::new(4, 0, 0)) };
    //         assert_eq!(
    //             west.edge().point.point,
    //             Point::new(SafeFloat::MIN.0, SafeFloat::MAX.0)
    //         );
    //         assert_eq!(
    //             west.sym().edge().point.point,
    //             Point::new(SafeFloat::MIN.0, SafeFloat::MIN.0)
    //         );

    //         let p1 = Point::new(1.0, 0.0);
    //         triangulation.add_point(p1);
    //         assert_eq!(triangulation.qeds().unwrap().quads.len(), 8);
    //         assert_eq!(triangulation.retriangulate_all(), 0);
    //     }

    //     #[test]
    //     fn two_point_triangulation() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let p1 = Point::new(0.0, 0.0);
    //         let p2 = Point::new(1.0, 1.0);
    //         triangulation.add_point(p1);
    //         triangulation.add_point(p2);
    //         for (i, edge) in triangulation.qeds().unwrap().base_edges().enumerate() {
    //             println!(
    //                 "Edge[{}]: {}-{}",
    //                 i,
    //                 edge.edge().point.point,
    //                 edge.sym().edge().point.point
    //             );
    //         }
    //         assert_eq!(triangulation.qeds().unwrap().quads.len(), 11);
    //         assert!(has_edge(&triangulation, p1, p2), "missing edge");
    //     }

    //     #[test]
    //     fn triangle_triangulation() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let p1 = Point::new(0.0, 0.0);
    //         let p2 = Point::new(5.0, 0.10);
    //         let p3 = Point::new(2.5, 5.1);
    //         triangulation.add_constraint(p1, p2);
    //         triangulation.add_constraint(p2, p3);
    //         triangulation.add_constraint(p3, p1);
    //         assert_eq!(triangulation.triangles().count(), 8);
    //         let tri_info = triangulation.classify_triangles();
    //         assert_eq!(tri_info.ns(), (1, 0, 7, 0));
    //         for (target, info) in tri_info.0 {
    //             println!("EdgeTarget: {:?} Info: {:?}", target, info);
    //         }
    //     }

    //     #[test]
    //     fn quad_triangulation() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let p1 = Point::new(0.0, 0.0);
    //         let p2 = Point::new(5.0, 0.10);
    //         let p3 = Point::new(2.5, 5.1);
    //         let p4 = Point::new(5.0, 5.0);
    //         triangulation.add_constraint(p1, p2);
    //         triangulation.add_constraint(p3, p1);
    //         triangulation.add_constraint(p4, p2);
    //         triangulation.add_constraint(p4, p3);
    //         // 10 triangles, 2 from the quad and 8 from the bounding box.
    //         assert_eq!(triangulation.triangles().count(), 10);
    //         let tri_info = triangulation.classify_triangles();
    //         assert_eq!(tri_info.ns(), (0, 2, 8, 0));
    //         for (target, info) in tri_info.0 {
    //             println!("EdgeTarget: {:?} Info: {:?}", target, info);
    //         }
    //     }

    //     #[test]
    //     fn star_triangulation() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let p1 = Point::new(0.0, 0.0);
    //         let p2 = Point::new(0.5, 0.5);
    //         let p3 = Point::new(-0.5, 0.5);
    //         let p4 = Point::new(-2.0, -1.0);
    //         let p5 = Point::new(2.0, -1.0);
    //         let p6 = Point::new(0.0, 2.0);
    //         triangulation.add_constraint(p4, p1);
    //         triangulation.add_constraint(p1, p5);
    //         triangulation.add_constraint(p5, p2);
    //         triangulation.add_constraint(p2, p6);
    //         triangulation.add_constraint(p6, p3);
    //         triangulation.add_constraint(p3, p4);

    //         // 14 triangles, 4 from the quad and 10 from the bounding box.
    //         assert_eq!(triangulation.triangles().count(), 14);
    //         let tri_info = triangulation.classify_triangles();
    //         assert_eq!(tri_info.ns(), (0, 7, 7, 0));
    //     }

    //     #[test]
    //     fn ring_triangulation() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let p1 = Point::new(-1.0, -1.0);
    //         let p2 = Point::new(1.0, -1.0);
    //         let p3 = Point::new(0.0, 1.0);
    //         let p4 = Point::new(-2.0, -1.5);
    //         let p5 = Point::new(2.0, -1.5);
    //         let p6 = Point::new(0.0, 2.0);
    //         triangulation.add_constraint(p1, p2);
    //         triangulation.add_constraint(p2, p3);
    //         triangulation.add_constraint(p3, p1);
    //         triangulation.add_constraint(p4, p5);
    //         triangulation.add_constraint(p5, p6);
    //         triangulation.add_constraint(p6, p4);

    //         // 14 triangles, 7 from the quad and 7 from the bounding box.
    //         assert_eq!(triangulation.triangles().count(), 14);
    //         let tri_info = triangulation.classify_triangles();
    //         assert_eq!(tri_info.ns(), (1, 0, 13, 0));
    //     }

    //     #[test]
    //     fn figure8_triangulation() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let p1 = Point::new(0.0, 0.0);
    //         let p2 = Point::new(5.0, 0.0);
    //         let p3 = Point::new(5.0, 3.0);
    //         let p4 = Point::new(0.0, 3.0);

    //         let q1 = Point::new(1.0, 1.0);
    //         let q2 = Point::new(2.0, 1.0);
    //         let q3 = Point::new(2.0, 2.0);
    //         let q4 = Point::new(1.0, 2.0);

    //         let r = Point::new(2.0, 0.0);

    //         let r1 = r + q1;
    //         let r2 = r + q2;
    //         let r3 = r + q3;
    //         let r4 = r + q4;

    //         triangulation.add_constraint(p1, p2);
    //         triangulation.add_constraint(p2, p3);
    //         triangulation.add_constraint(p3, p4);
    //         triangulation.add_constraint(p4, p1);

    //         triangulation.add_constraint(q1, q2);
    //         triangulation.add_constraint(q2, q3);
    //         triangulation.add_constraint(q3, q4);
    //         triangulation.add_constraint(q4, q1);

    //         triangulation.add_constraint(r1, r2);
    //         triangulation.add_constraint(r2, r3);
    //         triangulation.add_constraint(r3, r4);
    //         triangulation.add_constraint(r4, r1);

    //         // 26 triangles, 18 from the quad and 8 from the bounding box.
    //         assert_eq!(triangulation.triangles().count(), 26);
    //         let _tri_info = triangulation.classify_triangles();
    //         // assert_eq!(tri_info.ns(),(1,0,13,0));
    //     }

    //     #[test]
    //     fn l_shaped_triangulation() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         // triangulation.add_point(Point::new(-5.5, 5.5));
    //         // triangulation.add_point(Point::new(-4.5, 5.5));

    //         triangulation.add_constraint(Point::new(-5.5, 5.5), Point::new(-4.5, 5.5));
    //         triangulation.add_constraint(Point::new(-4.5, 5.5), Point::new(-4.5, 4.5));

    //         // 8 triangles, 1 from the l-shape and 7 from the bounding box.
    //         assert_eq!(triangulation.triangles().count(), 8);
    //         // let _tri_info = triangulation.classify_triangles();
    //         // assert_eq!(tri_info.ns(),(1,0,13,0));
    //     }

    //     #[test]
    //     fn u_shaped_triangulation() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         triangulation.add_constraint(Point::new(-5.5, 5.5), Point::new(-4.5, 5.5));
    //         triangulation.add_constraint(Point::new(-4.5, 5.5), Point::new(-4.5, 4.5));
    //         triangulation.add_constraint(Point::new(-5.5, 5.5), Point::new(-5.5, 4.5));

    //         // 26 triangles, 18 from the quad and 8 from the bounding box.
    //         // assert_eq!(triangulation.triangles().count(), 26);
    //         let _tri_info = triangulation.classify_triangles();
    //         // assert_eq!(tri_info.ns(),(1,0,13,0));
    //     }

    //     #[test]
    //     fn u2_shaped_triangulation() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         triangulation.add_constraint(Point::new(-5.5, 2.5), Point::new(-4.5, 2.5));
    //         triangulation.add_constraint(Point::new(-5.5, 1.5), Point::new(-4.5, 1.5));
    //         triangulation.add_constraint(Point::new(-5.5, 2.5), Point::new(-5.5, 1.5));
    //         triangulation.add_constraint(Point::new(-3.5, 2.5), Point::new(-3.5, 1.5));

    //         // 26 triangles, 18 from the quad and 8 from the bounding box.
    //         // assert_eq!(triangulation.triangles().count(), 26);
    //         let _tri_info = triangulation.classify_triangles();
    //         // assert_eq!(tri_info.ns(),(1,0,13,0));
    //     }

    //     #[test]
    //     fn split_edge() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         assert_eq!(triangulation.qeds().unwrap().quads.len(), 5);
    //         let p1 = Point::new(0.0, 0.0);
    //         let p2 = Point::new(5.0, 0.0);
    //         let p3 = Point::new(2.5, 5.1);
    //         triangulation.add_point(p1);
    //         triangulation.add_point(p2);
    //         triangulation.add_point(p3);
    //         assert_eq!(triangulation.qeds().unwrap().quads.len(), 14);
    //         triangulation.add_point(Point::new(2.5, 0.0));
    //         assert_eq!(triangulation.qeds().unwrap().quads.len(), 17);
    //     }

    //     fn valid_triangulation(triangulation: &SurfaceTriangulation) {
    //         for (i, _quad) in triangulation.qeds().unwrap().quads.iter() {
    //             let target1 = EdgeTarget::new(i, 0, 0);
    //             let target2 = EdgeTarget::new(i, 2, 0);
    //             unsafe {
    //                 assert_ne!(
    //                     triangulation
    //                         .qeds()
    //                         .unwrap()
    //                         .edge_a_ref(target1)
    //                         .onext()
    //                         .target(),
    //                     target1,
    //                     "Edge {:?} with points {} - {} has a floating end",
    //                     target1,
    //                     triangulation
    //                         .qeds()
    //                         .unwrap()
    //                         .edge_a_ref(target1)
    //                         .edge()
    //                         .point
    //                         .point,
    //                     triangulation
    //                         .qeds()
    //                         .unwrap()
    //                         .edge_a_ref(target1)
    //                         .sym()
    //                         .edge()
    //                         .point
    //                         .point
    //                 );
    //                 assert_ne!(
    //                     triangulation
    //                         .qeds()
    //                         .unwrap()
    //                         .edge_a_ref(target2)
    //                         .onext()
    //                         .target(),
    //                     target2,
    //                     "Edge {:?} with points {} - {} has a floating end",
    //                     target2,
    //                     triangulation
    //                         .qeds()
    //                         .unwrap()
    //                         .edge_a_ref(target2)
    //                         .edge()
    //                         .point
    //                         .point,
    //                     triangulation
    //                         .qeds()
    //                         .unwrap()
    //                         .edge_a_ref(target2)
    //                         .sym()
    //                         .edge()
    //                         .point
    //                         .point
    //                 );
    //             }
    //         }
    //     }

    //     #[test]
    //     fn detailed_triangulation() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let p1 = Point::new(0.0, 0.0);
    //         let p2 = Point::new(10.0, 1.0);
    //         let p3 = Point::new(5.0, 10.0);
    //         // let p4 = Point::new(11.0, 11.0);
    //         // let p5 = Point::new(15.0, 2.5);
    //         // let p6 = Point::new(12.0, 12.5);
    //         // let p7 = Point::new(15.0, 12.5);
    //         triangulation.add_point(p1);
    //         triangulation.add_point(p2);
    //         triangulation.add_point(p3);
    //         valid_triangulation(&triangulation);
    //         // triangulation.add_point(p4);
    //     }

    //     #[test]
    //     fn big_spiral_triangulation() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let a = 1.0;
    //         let n_max = 100;
    //         for n in 0..n_max {
    //             let t = (n as f64 / 6.0) * std::f64::consts::PI;
    //             let x = a * t * t.cos();
    //             let y = a * t * t.sin();
    //             triangulation.add_point(Point::new(x, y));
    //         }
    //         // We shouldn't need to retriangulate any more edges.
    //         assert_eq!(triangulation.retriangulate_all(), 0);
    //         valid_triangulation(&triangulation);
    //     }

    //     #[test]
    //     fn small_spiral_triangulation() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let p1 = Point::new(0.0, 0.0);
    //         let p2 = Point::new(0.45344999999999996, 0.261799);
    //         let p3 = Point::new(0.5235989999999999, 0.9068999999999999);
    //         let p4 = Point::new(0.408105, 0.235619);
    //         let p5 = Point::new(0.47123899999999996, 0.81621);
    //         let p6 = Point::new(0.0, 1.413717);
    //         triangulation.add_point(p1);
    //         triangulation.add_point(p2);
    //         triangulation.add_point(p3);
    //         assert!(has_edge(&triangulation, p1, p2));
    //         assert!(has_edge(&triangulation, p2, p3));
    //         assert!(has_edge(&triangulation, p3, p1));
    //         triangulation.add_point(p4);
    //         assert!(has_edge(&triangulation, p1, p4));
    //         assert!(has_edge(&triangulation, p4, p2));
    //         triangulation.add_point(p5);
    //         assert!(has_edge(&triangulation, p1, p5));
    //         assert!(has_edge(&triangulation, p5, p3));
    //         assert!(has_edge(&triangulation, p2, p3));
    //         assert!(has_edge(&triangulation, p2, p3));
    //         triangulation.add_point(p6);

    //         // We shouldn't need to retriangulate any more edges.
    //         assert_eq!(triangulation.retriangulate_all(), 0);
    //         valid_triangulation(&triangulation);
    //     }

    //     #[test]
    //     fn line_triangulation() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let p1 = Point::new(0.0, 0.0);
    //         let p2 = Point::new(5.0, 0.0);
    //         let p3 = Point::new(7.0, 0.0);
    //         triangulation.add_point(p1);
    //         triangulation.add_point(p2);
    //         triangulation.add_point(p3);
    //         assert_eq!(triangulation.qeds().unwrap().quads.len(), 14);
    //     }

    //     #[test]
    //     fn no_intersections() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let p1 = Point::new(0.0, 0.0);
    //         let p2 = Point::new(10.0, 0.0);
    //         let p3 = Point::new(5.0, 10.0);
    //         triangulation.add_point(p1);
    //         triangulation.add_point(p2);
    //         triangulation.add_point(p3);
    //         valid_triangulation(&triangulation);
    //         assert_eq!(
    //             triangulation
    //                 .find_intersections_between_points(Point::new(2.0, 2.0), Point::new(4.0, 3.0))
    //                 .map(|x| {
    //                     if let Intersection::Point(p) = x {
    //                         println!(
    //                             "IntersectionPoint: {}",
    //                             p.edge(&triangulation).point.point()
    //                         );
    //                     }
    //                     x
    //                 })
    //                 .collect::<Vec<_>>(),
    //             vec![]
    //         );
    //     }

    //     #[test]
    //     fn single_intersection() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let p1 = Point::new(0.0, 0.0);
    //         let p2 = Point::new(5.0, -1.0);
    //         let p3 = Point::new(10.0, 0.0);
    //         let p4 = Point::new(5.0, 1.0);
    //         triangulation.add_constraint(p1, p2);
    //         triangulation.add_constraint(p2, p3);
    //         triangulation.add_constraint(p3, p4);
    //         triangulation.add_constraint(p4, p1);
    //         triangulation.add_constraint(p1, p3);
    //     }

    //     #[test]
    //     fn double_intersection() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let p1 = Point::new(0.0, 0.0);
    //         let p2 = Point::new(5.0, -1.0);
    //         let p3 = Point::new(10.0, 0.0);
    //         let p4 = Point::new(5.5, 1.0);
    //         let p5 = Point::new(4.5, 1.0);
    //         triangulation.add_constraint(p1, p2);
    //         triangulation.add_constraint(p2, p3);
    //         triangulation.add_constraint(p3, p4);
    //         triangulation.add_constraint(p4, p5);
    //         triangulation.add_constraint(p5, p1);
    //         triangulation.add_constraint(p1, p3);
    //     }

    //     #[test]
    //     fn regress1() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let segments = vec![
    //             (
    //                 Point::new(-75.3645612469715, 10.649471266802962),
    //                 Point::new(11.264242662961536, -74.1427143080474),
    //             ),
    //             (
    //                 Point::new(-65.08752072884896, 85.23877185897558),
    //                 Point::new(-36.55529677285707, 25.159802183655742),
    //             ),
    //             (
    //                 Point::new(-71.01084733411076, 30.660749902036656),
    //                 Point::new(67.62855075915658, -89.10279376500583),
    //             ),
    //         ];
    //         for (p1, p2) in segments.iter() {
    //             triangulation.add_constraint(*p1, *p2);
    //         }
    //         valid_triangulation(&triangulation);
    //     }

    //     #[test]
    //     fn regress2() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let segments = vec![
    //             (
    //                 Point {
    //                     x: -40.292118735040106,
    //                     y: 82.0073097016039,
    //                 },
    //                 Point {
    //                     x: -42.07946675183436,
    //                     y: -44.51802168917966,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: -52.00978628295401,
    //                     y: 70.99340776302836,
    //                 },
    //                 Point {
    //                     x: 60.41164157780554,
    //                     y: 77.45771825005286,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: 87.60492632142518,
    //                     y: 13.83554980163639,
    //                 },
    //                 Point {
    //                     x: -64.2901267321343,
    //                     y: 68.80564964035884,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: 40.89557780924511,
    //                     y: 34.30455081240521,
    //                 },
    //                 Point {
    //                     x: -96.09647543301776,
    //                     y: 93.8999414141874,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: -45.63618305392634,
    //                     y: 73.74453797959046,
    //                 },
    //                 Point {
    //                     x: -44.04842079916467,
    //                     y: 19.50193303235106,
    //                 },
    //             ),
    //         ];
    //         for (p1, p2) in segments.iter() {
    //             triangulation.add_constraint(*p1, *p2);
    //         }
    //         valid_triangulation(&triangulation);
    //     }

    //     #[test]
    //     fn regress3() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let segments = vec![
    //             (
    //                 Point {
    //                     x: -0.4214359290103573,
    //                     y: 95.82816145300376,
    //                 },
    //                 Point {
    //                     x: -9.699326930357884,
    //                     y: 26.125595226387333,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: -13.147969911139384,
    //                     y: -54.45197788876022,
    //                 },
    //                 Point {
    //                     x: -44.33529262658724,
    //                     y: -75.64912335744887,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: -16.705404907959036,
    //                     y: 35.1977871928658,
    //                 },
    //                 Point {
    //                     x: 60.46183760891148,
    //                     y: 10.128488435745282,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: -59.73082630206363,
    //                     y: -57.035655725149056,
    //                 },
    //                 Point {
    //                     x: 49.69180679531209,
    //                     y: 96.08879129160505,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: 12.685519199969605,
    //                     y: 78.74650844233386,
    //                 },
    //                 Point {
    //                     x: -0.5177468902998896,
    //                     y: -41.830391168991895,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: -91.33303213105921,
    //                     y: 42.425551365690154,
    //                 },
    //                 Point {
    //                     x: 67.18820631566183,
    //                     y: 49.651436020714186,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: 19.547046919264545,
    //                     y: 21.14100619435102,
    //                 },
    //                 Point {
    //                     x: 43.67012362837994,
    //                     y: -56.81803213245602,
    //                 },
    //             ),
    //         ];
    //         for (p1, p2) in segments.iter() {
    //             triangulation.add_constraint(*p1, *p2);
    //         }
    //         valid_triangulation(&triangulation);
    //     }

    //     #[test]
    //     fn regress4() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let segments = vec![
    //             (
    //                 Point {
    //                     x: 49.354794208915905,
    //                     y: 77.9312265303424,
    //                 },
    //                 Point {
    //                     x: -11.492740263412088,
    //                     y: 8.956223279493656,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: 90.35804944791943,
    //                     y: -64.48450750858385,
    //                 },
    //                 Point {
    //                     x: 25.29536887506309,
    //                     y: 8.406416670169662,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: 46.09321063375256,
    //                     y: 61.22935053925707,
    //                 },
    //                 Point {
    //                     x: -26.439979432037305,
    //                     y: 51.94522246412245,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: 96.6694613097255,
    //                     y: 45.14085139658687,
    //                 },
    //                 Point {
    //                     x: 58.35546494466675,
    //                     y: 53.009402448096495,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: 58.99661219941143,
    //                     y: 65.20711681694809,
    //                 },
    //                 Point {
    //                     x: -63.66585352841398,
    //                     y: -59.20567981731186,
    //                 },
    //             ),
    //         ];
    //         for (p1, p2) in segments.iter() {
    //             triangulation.add_constraint(*p1, *p2);
    //         }
    //         valid_triangulation(&triangulation);
    //     }

    //     #[test]
    //     fn regress5() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let segments = vec![
    //             (
    //                 Point {
    //                     x: -51.400285762967044,
    //                     y: 78.08416289394077,
    //                 },
    //                 Point {
    //                     x: 73.94719688650767,
    //                     y: -9.952271890507092,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: -80.05915835704234,
    //                     y: 97.55389372322861,
    //                 },
    //                 Point {
    //                     x: 75.64795651512509,
    //                     y: -27.169794412596275,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: 81.99617838589236,
    //                     y: 8.336094155253178,
    //                 },
    //                 Point {
    //                     x: 99.33459619053124,
    //                     y: -28.76119989156618,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: 52.06631213621819,
    //                     y: 91.79642677407745,
    //                 },
    //                 Point {
    //                     x: 98.04001031185186,
    //                     y: -81.61857001089774,
    //                 },
    //             ),
    //         ];
    //         for (p1, p2) in segments.iter() {
    //             triangulation.add_constraint(*p1, *p2);
    //         }
    //         valid_triangulation(&triangulation);
    //     }

    //     #[test]
    //     fn regress6() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let segments = vec![
    //             (
    //                 Point {
    //                     x: -71.48547038818833,
    //                     y: -74.39322987310905,
    //                 },
    //                 Point {
    //                     x: -89.94096477395534,
    //                     y: 19.301100606945766,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: -86.9828479383787,
    //                     y: 0.05961466115694236,
    //                 },
    //                 Point {
    //                     x: 43.25946276642907,
    //                     y: -77.91468968915898,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: -4.142054219662228,
    //                     y: 70.85270987086068,
    //                 },
    //                 Point {
    //                     x: -90.96995169850732,
    //                     y: 36.77453032058807,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: -73.374254157982,
    //                     y: 87.53245031954737,
    //                 },
    //                 Point {
    //                     x: -46.959018139029226,
    //                     y: -48.21755053303174,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: 97.61862767011499,
    //                     y: 56.159609123931716,
    //                 },
    //                 Point {
    //                     x: -60.03876659083822,
    //                     y: -4.969101911120703,
    //                 },
    //             ),
    //             (
    //                 Point {
    //                     x: -99.10800691452306,
    //                     y: 43.80328990194849,
    //                 },
    //                 Point {
    //                     x: -58.187890315932435,
    //                     y: -11.073968979166835,
    //                 },
    //             ),
    //         ];
    //         for (p1, p2) in segments.iter() {
    //             triangulation.add_constraint(*p1, *p2);
    //         }
    //         valid_triangulation(&triangulation);
    //     }

    //     #[ignore]
    //     #[test]
    //     fn simple_intersections() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let p1 = Point::new(0.0, 0.0);
    //         let p2 = Point::new(10.0, 0.0);
    //         let p3 = Point::new(5.0, 10.0);
    //         let p4 = Point::new(20.0, 20.0);
    //         triangulation.add_point(p1);
    //         triangulation.add_point(p2);
    //         triangulation.add_point(p3);
    //         triangulation.add_point(p4);
    //         valid_triangulation(&triangulation);
    //         {
    //             // single intersection
    //             let intersections: Vec<_> = triangulation
    //                 .find_intersections_between_points(Point::new(2.0, 2.0), Point::new(9.0, 8.0))
    //                 .collect();
    //             println!("intersections: {:?}", intersections);
    //             assert_eq!(intersections.len(), 1);
    //             let first_intersection =
    //                 unsafe { triangulation.qeds.edge_a_ref(intersections[0].target()) };
    //             // The edges found have a specified direction.
    //             assert_eq!(
    //                 (
    //                     first_intersection.edge().point.point(),
    //                     first_intersection.sym().edge().point.point()
    //                 ),
    //                 (p2, p3)
    //             );
    //         }
    //         {
    //             // two intersections
    //             let intersections = triangulation
    //                 .find_intersections_between_points(Point::new(2.0, 2.0), Point::new(18.0, 6.0));
    //             assert_eq!(intersections.clone().count(), 2);
    //             let first_intersection = unsafe {
    //                 triangulation
    //                     .qeds
    //                     .edge_a_ref(intersections.clone().collect::<Vec<_>>()[0].target())
    //             };
    //             // The edges found have a specified direction.
    //             assert_eq!(
    //                 (
    //                     first_intersection.edge().point.point(),
    //                     first_intersection.sym().edge().point.point()
    //                 ),
    //                 (p2, p3)
    //             );
    //             let second_intersection = unsafe {
    //                 triangulation
    //                     .qeds
    //                     .edge_a_ref(intersections.collect::<Vec<_>>()[1].target())
    //             };
    //             assert_eq!(
    //                 (
    //                     second_intersection.edge().point.point(),
    //                     second_intersection.sym().edge().point.point()
    //                 ),
    //                 (p2, p4)
    //             );
    //         }
    //         let p5 = Point::new(7.5, 5.0);
    //         triangulation.add_point(p5);
    //         {
    //             // through a vertex (but starting within a triangle).
    //             println!("Testing through a vertex");
    //             let intersections = triangulation
    //                 .find_intersections_between_points(Point::new(6.0, 5.0), Point::new(9.0, 5.0));
    //             assert_eq!(intersections.count(), 0);
    //             // let first_intersection = unsafe {triangulation.qeds.edge_a_ref(intersections[0])};
    //             // // The edges found have a specified direction.
    //             // assert_eq!((first_intersection.edge().point.point(),first_intersection.sym().edge().point.point()),(p2,p3));
    //             // let second_intersection = unsafe {triangulation.qeds.edge_a_ref(intersections[1])};
    //             // assert_eq!((second_intersection.edge().point.point(),second_intersection.sym().edge().point.point()),(p2,p4));
    //         }
    //         {
    //             println!("Testing from a vertex (without leaving initial triangle)");
    //             let intersections = triangulation
    //                 .find_intersections_between_points(Point::new(0.0, 0.0), Point::new(1.0, 1.0));
    //             assert_eq!(intersections.count(), 0);
    //         }
    //         {
    //             // through a vertex (but starting within a triangle).
    //             println!("Testing from a vertex (with 1 intersection)");
    //             let intersections = triangulation
    //                 .find_intersections_between_points(Point::new(0.0, 0.0), Point::new(8.0, 10.0));
    //             assert_eq!(intersections.clone().count(), 1);

    //             let first_intersection = unsafe {
    //                 triangulation
    //                     .qeds
    //                     .edge_a_ref(intersections.collect::<Vec<_>>()[0].target())
    //             };
    //             assert_eq!(
    //                 (
    //                     first_intersection.edge().point.point(),
    //                     first_intersection.sym().edge().point.point()
    //                 ),
    //                 (p5, p3)
    //             );
    //         }
    //     }

    //     #[test]
    //     fn simple_constraints() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let p1 = Point::new(-10.0, 0.0);
    //         let p2 = Point::new(10.0, 0.0);
    //         // let p3 = Point::new(5.0, 10.0);
    //         triangulation.add_constraint(p1, p2);
    //         assert!(has_constraint(&triangulation, p1, p2));
    //         // triangulation.add_point(p2);
    //         // triangulation.add_point(p3);
    //         // valid_triangulation(&triangulation);
    //         // assert_eq!(
    //         // triangulation.find_intersections_between_points(Point::new(2.0, 2.0), Point::new(4.0, 3.0)),
    //         // vec![]
    //         // );
    //     }

    //     #[test]
    //     fn cross_constraints() {
    //         let mut triangulation = SurfaceTriangulation::new();
    //         let p1 = Point::new(-10.0, 0.0);
    //         let p2 = Point::new(10.0, 0.0);
    //         let p3 = Point::new(0.0, -10.0);
    //         let p4 = Point::new(0.0, 10.0);
    //         triangulation.add_constraint(p1, p2);
    //         triangulation.add_constraint(p3, p4);
    //         let p0 = Point::new(0.0, 0.0);
    //         assert!(has_constraint(&triangulation, p0, p1));
    //         assert!(has_constraint(&triangulation, p0, p2));
    //         assert!(has_constraint(&triangulation, p0, p3));
    //         assert!(has_constraint(&triangulation, p0, p4));
    //         let p5 = Point::new(3.0, 3.0);
    //         let p6 = Point::new(6.0, 6.0);
    //         triangulation.add_constraint(p5, p6);
    //         // assert!(has_constraint(&triangulation, p5, p6));
    //         // let n_constraints = triangulation
    //         //     .qeds
    //         //     .quads
    //         //     .iter()
    //         //     .filter(|(_i, q)| q.edges_a[0].point.constraint)
    //         //     .count();
    //         // assert_eq!(n_constraints, 9);
    //     }
}
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub enum IntersectionResult {
    Parallel,
    LineIntersection(f64, f64),
}

// fn segment_intersection(s1: (Point, Point), s2: (Point, Point)) -> IntersectionResult {
//     let p = s1.0;
//     let r_abs = s1.1;
//     let q = s2.0;
//     let s_abs = s2.1;
//     // line1 t = p + t*r
//     // line2 r = q + u*s
//     let r = r_abs - p;
//     let s = s_abs - q;

//     // Sometimes we get numerical errors when the points coincide, so we
//     // check for that first
//     let t = if p == q {
//         0.0
//     } else if p == s_abs {
//         0.0
//     } else if r_abs == q {
//         1.0
//     } else if r_abs == s_abs {
//         1.0
//     } else {
//         cross(q - p, s) / cross(r, s)
//     };

//     let u = if q == p {
//         0.0
//     } else if q == r_abs {
//         0.0
//     } else if s_abs == p {
//         1.0
//     } else if s_abs == r_abs {
//         1.0
//     } else {
//         cross(q - p, r) / cross(r, s)
//     };

//     // If r `cross` s is 0 then the lines are parallel and do not intersect, so
//     // return nothing.
//     let cross_val = cross(r, s);
//     if cross_val == 0.0 {
//         IntersectionResult::Parallel
//     } else {
//         IntersectionResult::LineIntersection(t, u)
//     }
// }

fn cross(pa: Point, pb: Point) -> f64 {
    pa.x * pb.y - pb.x * pa.y
}

fn show_triangle<T>(triangle: EdgeRefA<Segment<T>, Space>) -> String {
    format!(
        "{}-{}-{}",
        triangle.edge().point.point(),
        triangle.l_next().edge().point.point(),
        triangle.l_next().l_next().edge().point.point()
    )
}

fn show_edge<T>(edge: EdgeRefA<Segment<T>, Space>) -> String {
    format!(
        "{}-{}",
        edge.edge().point.point(),
        edge.sym().edge().point.point(),
    )
}

/// Assert that every node as a component.
pub fn debug_assert_spaces<T>(triangulation: &SurfaceTriangulation<T>) {
    println!("Debugging Spokes");
    for (i, quad) in triangulation.qeds.quads.iter() {
        println!(
            "DebugBefore[{}]: {}-{} {:?}-{:?}",
            i,
            quad.edges_a[0].point.point,
            quad.edges_a[1].point.point,
            quad.edges_b[1].point,
            quad.edges_b[0].point,
        );
    }
    // Get all B edges.
    let mut b_targets = Vec::new();
    for (i, quad) in triangulation.qeds.quads.iter() {
        b_targets.push(EdgeTarget::new(i, 1, 0));
        b_targets.push(EdgeTarget::new(i, 3, 0));
    }
    for b_target in b_targets {
        let mut current = unsafe { triangulation.qeds.edge_b_ref(b_target) };
        let first = current.target();
        let space = current.edge().point;
        loop {
            if current.edge().point != space {
                panic!("invalid spoke");
            }
            current = current.oprev();
            if current.target() == first {
                break;
            }
        }
    }
}
