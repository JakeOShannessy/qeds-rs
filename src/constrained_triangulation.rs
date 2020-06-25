use crate::point::*;
use crate::qeds::*;
use nalgebra::Matrix4;

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Segment {
    /// The point at the origin of this edge.
    pub point: Point,
    /// Indicates whether this segment is a constraint. If it is a constraint it
    /// cannot be swapped.
    pub constraint: bool,
}

impl Segment {
    pub fn new(point: Point) -> Self {
        Self {
            point,
            constraint: false,
        }
    }

    pub fn new_constraint(point: Point) -> Self {
        Self {
            point,
            constraint: true,
        }
    }
}

impl HasPoint for Segment {
    fn point(&self) -> Point {
        self.point
    }
}

#[derive(Clone, Debug)]
/// A Qeds data structure specialised to a 2d triangulation.
pub struct ConstrainedTriangulation {
    /// The quad-edge data structure we use as the basis for the triangulation.
    pub qeds: Qeds<Segment, ()>,
    pub boundary_edge: EdgeTarget,
    pub bounds: Option<(Point, Point)>,
}

impl ConstrainedTriangulation {
    pub fn new() -> Self {
        let mut qeds = Qeds::new();
        let sw = Point::new(SafeFloat::MIN.0, SafeFloat::MIN.0);
        let se = Point::new(SafeFloat::MAX.0, SafeFloat::MIN.0);
        let ne = Point::new(SafeFloat::MAX.0, SafeFloat::MAX.0);
        let nw = Point::new(SafeFloat::MIN.0, SafeFloat::MAX.0);
        let south = qeds
            .make_edge_with_a(Segment::new(sw), Segment::new(se))
            .target();
        let east = qeds
            .make_edge_with_a(Segment::new(se), Segment::new(ne))
            .target();
        unsafe {
            qeds.splice(east, south.sym());
            let centre = qeds.connect(east, south).target();

            let north = qeds
                .make_edge_with_a(Segment::new(ne), Segment::new(nw))
                .target();
            qeds.splice(north, east.sym());
            qeds.connect(north, centre.sym()).target();

            ConstrainedTriangulation {
                qeds,
                boundary_edge: south,
                bounds: None,
            }
        }
    }

    pub fn qeds(&self) -> Option<&Qeds<Segment, ()>> {
        Some(&self.qeds)
    }

    pub fn some_edge_a(&self) -> Option<EdgeRefA<Segment, ()>> {
        let (i, _) = self.qeds.quads.iter().next()?;
        unsafe { Some(self.qeds.edge_a_ref(EdgeTarget::new(i, 0, 0))) }
    }

    // TODO: This has not yet been proved to be stable. It may also loop
    // inifintely, particularly with constrianed triangulations.
    pub fn locate(&self, point: Point) -> Option<EdgeRefA<Segment, ()>> {
        use rand::Rng;
        let mut e = self.some_edge_a().unwrap();
        let mut rng = rand::thread_rng();
        let mut current_iterations = 0;
        loop {
            current_iterations += 1;
            if current_iterations > 2000 {
                panic!("locating failed for: {}", point);
            }
            // A random variable to determine if onext is tested first. If not
            // it is tested second.
            let onext_first: bool = rng.gen();
            if point == e.edge().point.point {
                return Some(e);
            } else if point == e.sym().edge().point.point {
                // If we lie on the eDest, return eSym(). This is a departure
                // from the paper which just returns e.
                return Some(e.sym());
            } else if e.lies_right_strict(point) {
                e = e.sym();
            } else {
                if onext_first {
                    if !e.onext().lies_right_strict(point) {
                        e = e.onext();
                    } else if !e.d_prev().lies_right_strict(point) {
                        e = e.d_prev();
                    } else {
                        return Some(e);
                    }
                } else {
                    if !e.d_prev().lies_right_strict(point) {
                        e = e.d_prev();
                    } else if !e.onext().lies_right_strict(point) {
                        e = e.onext();
                    } else {
                        return Some(e);
                    }
                }
            }
        }
    }

    pub fn add_point(&mut self, mut point: Point) -> Option<EdgeTarget> {
        point.snap();
        self.update_bounds(point);
        if let Some(edge_ref) = self.locate(point) {
            if edge_ref.edge().point.point == point || edge_ref.sym().edge().point.point == point {
                // If point lies on another point, just return the located edge.
                Some(edge_ref.target())
            // } else if point_lies_on_edge {
            //     // Delete that edge and replace it with 2 edges
            } else {
                let edge_target = edge_ref.target();
                Some(self.add_to_l_face(edge_target, point))
            }
        } else {
            // Point was out of bounds (probably)
            None
        }
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

    // TODO: remove the recursion from this function as we don't have tail
    // recursions.
    /// Edges are returned with each edge pointing from right to left across the
    /// intersection line.
    fn find_intersections(&self, a: Point, b: Point) -> Vec<EdgeTarget> {
        use Direction::*;
        // TODO: this belongs in the base trianglulation.
        // [`start`] is an edge of the triangle in which the fist point ([`a`])
        // is located.
        let mut start = self.locate(a).unwrap();
        println!(
            "Start is: {}-{}",
            start.edge().point.point(),
            start.sym().edge().point.point()
        );
        let mut intersecting_edge = if a == start.edge().point.point() {
            // TODO: there is an issue here as we might not be looking at the
            // right triangle. This applies both to the initial point and
            // subsequent "restarts".

            // First we need to place ourselves on the correct triangle. To do
            // that we iterate around the point until we are sure that the line
            // passes through a given triangle. This is indicated by subsequent
            // "spokes" switching from left to right. If the line passes through
            // a vertex we can deal with it here.
            let initial_start = start;
            loop {
                // End condition.
                if start.edge().point.point() == b {return vec![];}
                println!("looping top: {}-{}",start.edge().point.point(), start.sym().edge().point.point(),);
                let start_dir = left_or_right(
                    start.edge().point.point(),
                    start.sym().edge().point.point(),
                    b,
                );
                match start_dir {
                    Direction::Left => {
                        println!("It is to the left");
                        loop {
                            let next_edge = start.onext();
                            println!("looping next: {}-{}",next_edge.edge().point.point(), next_edge.sym().edge().point.point(),);
                            let next_dir = left_or_right(
                                next_edge.edge().point.point(),
                                next_edge.sym().edge().point.point(),
                                b,
                            );
                            match next_dir {
                                Straight => {
                                    return self.find_intersections(
                                        next_edge.sym().edge().point.point(),
                                        b,
                                    )
                                }
                                Left => {
                                    start = next_edge;
                                }
                                Right => break,
                            }
                        }
                        break;
                    }
                    _ => {
                        println!("It is straight or to the right");
                        start = start.onext();
                        if start == initial_start {
                            panic!("looped around ring");
                        }
                    }
                }
            }

            //Point a is on startOrg. This changes our initial question
            // somewhat. We know that the line either intersects the edge
            // opposite a, or also passes through one of the other vertices.
            let opposite_edge = start.l_next();
            // Does the line pass through the vertex to the right?
            let right_point = opposite_edge.edge().point.point();
            let passes_through_right = left_or_right(a, right_point, b) == Direction::Straight;
            if passes_through_right {
                // If so we need to start the calculation from that vertex.
                return self.find_intersections(a, right_point);
            }
            // Does the line pass through the vertex to the left?
            let left_point = opposite_edge.sym().edge().point.point();
            let passes_through_left = left_or_right(a, left_point, b) == Direction::Straight;
            if passes_through_left {
                // If so we need to start the calculation from that vertex.
                return self.find_intersections(a, left_point);
            }
            // If it passes through neither, then the intersecting edges is
            // simply the opposite edge.
            opposite_edge
        } else {
            // In the case of the first triangle we must test agains three lines
            // made by combining a with each of the vertices of the triangle. These
            // lines are ax,ay, and az in CCW order.
            let x = start.edge().point.point();
            let y = start.l_next().edge().point.point();
            let z = start.l_next().l_next().edge().point.point();
            println!("x is: {}", x);
            println!("y is: {}", y);
            println!("z is: {}", z);
            let ax = (a, x);
            let ay = (a, y);
            let az = (a, z);

            let right_of_ax = match left_or_right(ax.0, ax.1, b) {
                Direction::Straight => {
                    return self.find_intersections(x, b);
                }
                other => other,
            };
            let right_of_ay = match left_or_right(ay.0, ay.1, b) {
                Direction::Straight => {
                    return self.find_intersections(y, b);
                }
                other => other,
            };
            let right_of_az = match left_or_right(az.0, az.1, b) {
                Direction::Straight => {
                    return self.find_intersections(z, b);
                }
                other => other,
            };
            println!(
                "{} is to the {:?} of {}-{}",
                b,
                left_or_right(ax.0, ax.1, b),
                ax.0,
                ax.1
            );
            println!(
                "{} is to the {:?} of {}-{}",
                b,
                left_or_right(ay.0, ay.1, b),
                ay.0,
                ay.1
            );
            println!(
                "{} is to the {:?} of {}-{}",
                b,
                left_or_right(az.0, az.1, b),
                az.0,
                az.1
            );

            // Find the point at which the direction changes from left to right.
            if (right_of_ax, right_of_ay) == (Left, Right) {
                start
            } else if (right_of_ay, right_of_az) == (Left, Right) {
                start.l_next()
            } else if (right_of_az, right_of_ax) == (Left, Right) {
                start.l_next().l_next()
            } else {
                unreachable!()
            }
        };
        let mut edges: Vec<EdgeTarget> = Vec::new();
        // Now we iterate through all of the triangles.
        loop {
            // If the final point is to the left of the intersecting edge, it is
            // within the triangle, so we don't include the intersection and can
            // break.
            match left_or_right(
                /// TODO: b should be last not first.
                b,
                intersecting_edge.edge().point.point(),
                intersecting_edge.sym().edge().point.point(),
            ) {
                Direction::Left => {
                    println!(
                        "Internal: {} is to the left of {}-{}",
                        b,
                        intersecting_edge.edge().point.point(),
                        intersecting_edge.sym().edge().point.point()
                    );
                    // TODO: should probably be continue
                    break;
                }
                _ => (),
            }
            edges.push(intersecting_edge.target());
            // To get the same edge on the other triangle, we simply need the
            // edge in the opposite direction.
            let e = intersecting_edge.sym();
            // Get the point opposite.
            let c = e.l_next().sym().edge().point.point();
            match left_or_right(a, c, b) {
                Direction::Left => {
                    intersecting_edge = e.l_next().l_next();
                }
                Direction::Straight => {
                    let mut other_intersections = self.find_intersections(c, b);
                    edges.append(&mut other_intersections);
                    return edges;
                }
                Direction::Right => {
                    intersecting_edge = e.l_next();
                }
            }
        }
        edges
    }

    /// A poorly named function, but given a point described by an edge leading
    /// away from it and another point find the edge that lies to the right of a
    /// line from point A to point B.
    fn find_edge_to_right<'a>(
        &self,
        mut edge: EdgeRefA<'a, Segment, ()>,
        b: Point,
    ) -> EdgeRefA<'a, Segment, ()> {
        let initial_edge = edge;
        loop {
            println!("looping top");
            let start_dir = left_or_right(
                edge.edge().point.point(),
                edge.sym().edge().point.point(),
                b,
            );
            match start_dir {
                Direction::Left => {
                    println!("It is to the left");
                    loop {
                        println!("looping next");
                        let next_edge = edge.onext();
                        let next_dir = left_or_right(
                            next_edge.edge().point.point(),
                            next_edge.sym().edge().point.point(),
                            b,
                        );
                        match next_dir {
                            Direction::Straight => todo!("can't handle straight yet"),
                            Direction::Left => {
                                edge = next_edge;
                            }
                            Direction::Right => break,
                        }
                    }
                    break;
                }
                _ => {
                    println!("It is straight or to the right");
                    edge = edge.onext();
                    if edge == initial_edge {
                        panic!("looped around ring");
                    }
                }
            }
        }
        edge
    }

    pub fn add_constraint(&mut self, mut pa: Point, mut pb: Point) -> Option<()> {
        let pa_edge = self.add_point(pa)?;
        let pb_edge = self.add_point(pb)?;

        unsafe {
            pa = self.qeds.edge_a_ref(pa_edge).edge().point.point();
            pb = self.qeds.edge_a_ref(pb_edge).edge().point.point();
        }

        self.update_bounds(pa);
        self.update_bounds(pb);

        // As part of this process, find all the edges that intersect with our
        // new segment is an important part. This can be done relatively
        // efficiently.
        let intersecting_edges = self.find_intersections(pa, pb);

        // Step 1. Wherever the new segment crosses an existing constrained
        // segment add a new point. Also update the edges accordingly. As part
        // of this process, the constraint will be split into a series of
        // smaller constraints.
        println!("Inserting constraint from {} to {}", pa, pb);
        loop {
            let new_pa = loop {
                if let Some(&edge_target) = intersecting_edges.get(0) {
                    let e = unsafe { self.qeds.edge_a_ref(edge_target) };
                    if e.edge().point.constraint {
                        // It is critical that whenever we do this, we also account
                        // for it in the constraint we're inserting, as it will
                        // subtly modify it. Essentially it will split it into many
                        // smaller constraints. [`px`] is the new point at the
                        // intersection of the constraints. We need to be able to
                        // split a constraint too.

                        // First find/choose the intersection point.
                        let mut px: Point = match segment_intersection((e.edge().point.point(),e.sym().edge().point.point()), (pa,pb)) {
                            IntersectionResult::Parallel => panic!("unexpected parallel"),
                            IntersectionResult::LineIntersection(p) => p,
                        };
                        drop(e);
                        // Then insert the point.
                        px = unsafe {
                            let target = self.add_point(px).unwrap();
                            self.qeds
                                .edge_a_ref(target)
                                .edge()
                                .point
                                .point()
                        };

                        // We have now created a smaller constraint that we know
                        // crosses no other constraints (although it may of course
                        // cross other unconstrained edges). Wherever a an existing
                        // unconstrained segment crosses the new segment, flip it so
                        // that it aligns with the constraint.
                        let s_intersecting = self.find_intersections(px, pb);
                        // TODO: this part will require an understanding of the QEDS.
                        for s_e in s_intersecting {
                            let s_e_ref = unsafe { self.qeds.edge_a_ref(s_e) };
                            let s_e_constraint = s_e_ref.edge().point.constraint;
                            drop(s_e_ref);
                            if !s_e_constraint {
                                unsafe {
                                    // Swap the unconstrained edge.
                                    {
                                        self.swap(s_e);
                                    }
                                }

                            } else {
                                panic!("we should not be intersecting with any unconstrained edges");
                            }
                        }
                        // unreachable!("we should have made our constraint and exited already.");
                        // Step 3. Make and connect the constraint piece to the triangulation.
                        // let constraint_piece = self
                        //     .qeds
                        //     .make_edge_with_a(Segment::new_constraint(pa), Segment::new_constraint(px))
                        //     .target();
                        // Get the edge that will become immediately to the right of the
                        // constraint around pa.
                        // unsafe {
                        //     let ea = self
                        //         .find_edge_to_right(self.qeds.edge_a_ref(pa_edge), pb)
                        //         .target();
                        //     let eb = self
                        //         .find_edge_to_right(self.qeds.edge_a_ref(pb_edge), pa)
                        //         .target();
                        //     // TODO: where to splice. This seems like a good start but
                        //     // may be completely incorrect.
                        //     self.qeds.splice(constraint_piece, ea);
                        //     self.qeds.splice(constraint_piece.sym(), eb);
                        // }
                        break Some(px);
                    }

                } else {
                    println!("No intersecting edges");
                    break None;
                };
            };
            // After this process there should be an unconstrained edge
            // between our two points. We need to change it to be a
            // constraint. We know an existing edge for pa, we just need
            // to loop around pa until we find one with the destination
            // px.
            unsafe {
                let initial_edge = self.qeds.edge_a_ref(pa_edge).target();
                let mut edge = initial_edge;
                loop {
                    println!("Looking at edge {}-{}", pa, self.qeds.edge_a_ref(edge).sym().edge().point.point());
                    // TODO: this is not quite right
                    if self.qeds.edge_a_ref(edge).sym().edge().point.point() == pb {
                        self.qeds.edge_a_mut(edge).point.constraint = true;
                        break;
                    } else {
                        edge = self.qeds.edge_a_ref(edge).onext().target();
                        if edge == initial_edge {
                            panic!("looped around ring");
                        }
                    }
                }
            }
            if let Some(new_pa) = new_pa {
                pa = new_pa;
            } else {
                break;
            }
        }
        Some(())
    }

    /// This function accounts for the point lying on an existing edge.
    fn add_to_l_face(&mut self, mut edge_target: EdgeTarget, point: Point) -> EdgeTarget {
        unsafe {
            let point_a = self.qeds.edge_a_ref(edge_target).edge().point.point;
            let point_b = self.qeds.edge_a_ref(edge_target).sym().edge().point.point;
            let mut reinstate_as_constraint = false;
            // println!("Edge Target: {} - {}", self.qeds.edge_a_ref(edge_target).edge().point,self.qeds.edge_a_ref(edge_target).sym().edge().point);
            if point_a == point {
                return edge_target;
            } else if point_b == point {
                return edge_target.sym();
            } else if self.qeds.edge_a_ref(edge_target).lies_right(point) == Lies::On {
                // println!("Lies on edge: {} - {}", self.qeds.edge_a_ref(edge_target).edge().point,self.qeds.edge_a_ref(edge_target).sym().edge().point);
                let oprev_ref = self.qeds.edge_a_ref(edge_target).oprev();
                // We need to remember if this edge was constrained so
                // that we can reinstate it.
                reinstate_as_constraint = oprev_ref.edge().point.constraint;
                let oprev = oprev_ref.target();
                drop(oprev_ref);

                self.qeds.delete(edge_target);
                edge_target = oprev;
            }
            let first = self.qeds.edge_a_ref(edge_target).edge().point.point;
            let mut base = self
                .qeds
                .make_edge_with_a(Segment::new(first), Segment::new(point))
                .target();
            let return_value = base.sym();
            self.qeds.splice(base, edge_target);
            loop {
                let base_ref = self.qeds.connect(edge_target, base.sym());
                let this_is_constraint = reinstate_as_constraint
                    && (base_ref.sym().edge().point.point() == point_a
                        || base_ref.sym().edge().point.point() == point_b);
                edge_target = base_ref.oprev().target();
                base = base_ref.target();
                if this_is_constraint {
                    // TODO: this constraint setting code is not of great design
                    // and sorely needs checking.
                    self.qeds.edge_a_mut(base).point.constraint = true;
                    self.qeds.edge_a_mut(base.sym()).point.constraint = true;
                }
                if self.qeds.edge_a(edge_target.sym()).point.point == first {
                    break;
                }
            }
            edge_target = self.qeds.edge_a_ref(base).oprev().target();
            // The suspect edges are e(.Onext.Lprev)^k for k=0,1,2,3...
            let mut e = edge_target;
            // println!("Start Swap Loop: first: {}", first);
            loop {
                let t = self.qeds.edge_a_ref(e).oprev().target();
                let t_dest = self.qeds.edge_a_ref(t).sym().edge().point.point;
                let e_dest = self.qeds.edge_a_ref(e).sym().edge().point.point;
                let e_org = self.qeds.edge_a_ref(e).edge().point.point;
                // print!("Inspecting Edge for swap: {} {} : [{},{},{},{}] : ", e_org, e_dest, e_org, t_dest, e_dest, point);

                // TODO: we need to be cautious of infinite loops now that we're constrained.
                if self.qeds.edge_a_ref(e).lies_right_strict(t_dest)
                    && del_test_ccw(e_org, t_dest, e_dest, point)
                    && (self.qeds.edge_a_ref(e).edge().point.constraint == false)
                {
                    // println!("swap");
                    self.swap(e);
                    // This is different from the algorithm in the papaer
                    e = self.qeds.edge_a_ref(e).oprev().target();
                } else if e_org == first {
                    // println!("end");
                    break;
                } else {
                    // println!("don't swap");
                    e = self.qeds.edge_a_ref(e).onext().l_prev().target();
                }
            }
            return_value
        }
    }

    fn retriangulate_l_face(&mut self, edge_target: EdgeTarget) {
        unsafe { todo!() }
    }

    /// Assumes the EdgeTargets have already been put in the correct order.
    pub unsafe fn add_external_point_ordered(
        &mut self,
        boundary_edges: Vec<EdgeTarget>,
        point: Point,
    ) {
        let mut edges = boundary_edges.into_iter();
        let (g, e) = {
            let (a, b) = self.add_external_point(edges.next().unwrap(), point);
            (a.target(), b.target())
        };
        for f in edges {
            self.qeds.connect(e.sym(), f.sym());
        }
        self.boundary_edge = g;
    }

    pub unsafe fn add_external_point_unordered(
        &mut self,
        mut boundary_edges: Vec<EdgeTarget>,
        point: Point,
    ) {
        // Reorder the edges so that their base points are CCW around the point.
        boundary_edges.sort_by(|a, b| {
            let pa = self.qeds.edge_a_ref(*a).edge().point.point;
            let pb = self.qeds.edge_a_ref(*b).edge().point.point;
            match left_or_right(point, pb, pa) {
                Direction::Left => std::cmp::Ordering::Less,
                Direction::Straight => std::cmp::Ordering::Equal,
                Direction::Right => std::cmp::Ordering::Greater,
            }
        });
        self.add_external_point_ordered(boundary_edges, point);
    }

    pub fn boundary(&self) -> BoundaryIter<Segment, ()> {
        BoundaryIter::new(&self.qeds, self.boundary_edge)
    }

    pub unsafe fn add_external_point(
        &mut self,
        boundary_edge: EdgeTarget,
        p: Point,
    ) -> (EdgeRefA<Segment, ()>, EdgeRefA<Segment, ()>) {
        let dest = self
            .qeds
            .edge_a_ref(boundary_edge)
            .sym()
            .edge()
            .point
            .point
            .clone();
        let e = self
            .qeds
            .make_edge_with_a(Segment::new(p), Segment::new(dest))
            .target();
        self.qeds.splice(e.sym(), boundary_edge.sym());
        let f = self.qeds.connect(boundary_edge.sym(), e).target();
        (self.qeds.edge_a_ref(f), self.qeds.edge_a_ref(e))
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

    /// Perform Delaunay swapping on the entire triangulation until complete.
    /// Should not be necessary, mainly included for testing.
    pub fn retriangulate_all(&mut self) -> usize {
        let mut iterations = 0;
        let mut total_swaps = 0;
        loop {
            if iterations > 100 {
                panic!("too many triangulation iterations");
            }
            let swaps = self.retriangulate_all_single_pass();
            total_swaps += swaps;
            if swaps == 0 {
                break total_swaps;
            }
            iterations += 1;
        }
    }

    /// Warning: this is very inefficient and just for testing.
    fn retriangulate_all_single_pass(&mut self) -> usize {
        let mut swaps = 0;
        let edge_targets: Vec<EdgeTarget> =
            self.qeds.base_edges().map(|edge| edge.target()).collect();
        for e in edge_targets.into_iter() {
            unsafe {
                if self.del_test(e) && (self.qeds.edge_a_ref(e).edge().point.constraint == false) {
                    swaps += 1;
                    self.swap(e);
                }
            }
        }
        swaps
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
    // // println!("Point: {} {} {} {}",a,b,c,d);
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

fn determinant_2x2(a: f64, b: f64, c: f64, d: f64) -> f64 {
    a * d - b * c
}

pub fn has_edge(triangulation: &ConstrainedTriangulation, pa: Point, pb: Point) -> bool {
    get_edge(triangulation, pa, pb).is_some()
}

pub fn has_constraint(triangulation: &ConstrainedTriangulation, pa: Point, pb: Point) -> bool {
    if let Some(e) = get_edge(triangulation, pa, pb) {
        e.edge().point.constraint
    } else {
        false
    }
}

pub fn get_edge(
    triangulation: &ConstrainedTriangulation,
    pa: Point,
    pb: Point,
) -> Option<EdgeRefA<'_, Segment, ()>> {
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
    fn one_point_triangulation_only() {
        let mut triangulation = ConstrainedTriangulation::new();
        let p1 = Point::new(0.0, 0.0);
        triangulation.add_point(p1);
    }

    #[test]
    fn one_point_triangulation_location() {
        let mut triangulation = ConstrainedTriangulation::new();
        // assert_eq!(triangulation.qeds().unwrap().quads.len(), 5);
        let south = unsafe { triangulation.qeds.edge_a_ref(EdgeTarget::new(0, 0, 0)) };
        assert_eq!(
            south.edge().point.point,
            Point::new(SafeFloat::MIN.0, SafeFloat::MIN.0)
        );
        assert_eq!(
            south.sym().edge().point.point,
            Point::new(SafeFloat::MAX.0, SafeFloat::MIN.0)
        );

        assert_eq!(south.target(), south.l_next().l_next().l_next().target());
        assert_eq!(
            south.l_next().l_next().oprev().target(),
            EdgeTarget::new(3, 0, 0)
        );
        assert_eq!(
            south.l_next().l_next().oprev().l_next().target(),
            EdgeTarget::new(4, 0, 0)
        );

        let east = unsafe { triangulation.qeds.edge_a_ref(EdgeTarget::new(1, 0, 0)) };
        assert_eq!(
            east.edge().point.point,
            Point::new(SafeFloat::MAX.0, SafeFloat::MIN.0)
        );
        assert_eq!(
            east.sym().edge().point.point,
            Point::new(SafeFloat::MAX.0, SafeFloat::MAX.0)
        );

        let centre = unsafe { triangulation.qeds.edge_a_ref(EdgeTarget::new(2, 0, 0)) };
        assert_eq!(
            centre.edge().point.point,
            Point::new(SafeFloat::MAX.0, SafeFloat::MAX.0)
        );
        assert_eq!(
            centre.sym().edge().point.point,
            Point::new(SafeFloat::MIN.0, SafeFloat::MIN.0)
        );

        let north = unsafe { triangulation.qeds.edge_a_ref(EdgeTarget::new(3, 0, 0)) };
        assert_eq!(
            north.edge().point.point,
            Point::new(SafeFloat::MAX.0, SafeFloat::MAX.0)
        );
        assert_eq!(
            north.sym().edge().point.point,
            Point::new(SafeFloat::MIN.0, SafeFloat::MAX.0)
        );

        let west = unsafe { triangulation.qeds.edge_a_ref(EdgeTarget::new(4, 0, 0)) };
        assert_eq!(
            west.edge().point.point,
            Point::new(SafeFloat::MIN.0, SafeFloat::MAX.0)
        );
        assert_eq!(
            west.sym().edge().point.point,
            Point::new(SafeFloat::MIN.0, SafeFloat::MIN.0)
        );

        let p1 = Point::new(1.0, 0.0);
        triangulation.add_point(p1);
        assert_eq!(triangulation.qeds().unwrap().quads.len(), 8);
        assert_eq!(triangulation.retriangulate_all(), 0);
    }

    #[test]
    fn two_point_triangulation() {
        let mut triangulation = ConstrainedTriangulation::new();
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(1.0, 1.0);
        triangulation.add_point(p1);
        triangulation.add_point(p2);
        for (i, edge) in triangulation.qeds().unwrap().base_edges().enumerate() {
            println!(
                "Edge[{}]: {}-{}",
                i,
                edge.edge().point.point,
                edge.sym().edge().point.point
            );
        }
        assert_eq!(triangulation.qeds().unwrap().quads.len(), 11);
        assert!(has_edge(&triangulation, p1, p2), "missing edge");
    }

    #[test]
    fn triangle_triangulation() {
        let mut triangulation = ConstrainedTriangulation::new();
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(5.0, 0.10);
        let p3 = Point::new(2.5, 5.1);
        triangulation.add_point(p1);
        triangulation.add_point(p2);
        triangulation.add_point(p3);
    }

    #[test]
    fn split_edge() {
        let mut triangulation = ConstrainedTriangulation::new();
        assert_eq!(triangulation.qeds().unwrap().quads.len(), 5);
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(5.0, 0.0);
        let p3 = Point::new(2.5, 5.1);
        triangulation.add_point(p1);
        triangulation.add_point(p2);
        triangulation.add_point(p3);
        assert_eq!(triangulation.qeds().unwrap().quads.len(), 14);
        triangulation.add_point(Point::new(2.5, 0.0));
        assert_eq!(triangulation.qeds().unwrap().quads.len(), 17);
    }

    fn valid_triangulation(triangulation: &ConstrainedTriangulation) {
        for (i, quad) in triangulation.qeds().unwrap().quads.iter() {
            let target1 = EdgeTarget::new(i, 0, 0);
            let target2 = EdgeTarget::new(i, 2, 0);
            unsafe {
                assert_ne!(
                    triangulation
                        .qeds()
                        .unwrap()
                        .edge_a_ref(target1)
                        .onext()
                        .target(),
                    target1,
                    "Edge {:?} with points {} - {} has a floating end",
                    target1,
                    triangulation
                        .qeds()
                        .unwrap()
                        .edge_a_ref(target1)
                        .edge()
                        .point
                        .point,
                    triangulation
                        .qeds()
                        .unwrap()
                        .edge_a_ref(target1)
                        .sym()
                        .edge()
                        .point
                        .point
                );
                assert_ne!(
                    triangulation
                        .qeds()
                        .unwrap()
                        .edge_a_ref(target2)
                        .onext()
                        .target(),
                    target2,
                    "Edge {:?} with points {} - {} has a floating end",
                    target2,
                    triangulation
                        .qeds()
                        .unwrap()
                        .edge_a_ref(target2)
                        .edge()
                        .point
                        .point,
                    triangulation
                        .qeds()
                        .unwrap()
                        .edge_a_ref(target2)
                        .sym()
                        .edge()
                        .point
                        .point
                );
            }
        }
    }

    #[test]
    fn detailed_triangulation() {
        let mut triangulation = ConstrainedTriangulation::new();
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
    fn big_spiral_triangulation() {
        let mut triangulation = ConstrainedTriangulation::new();
        let a = 1.0;
        let n_max = 100;
        for n in 0..n_max {
            let t = (n as f64 / 6.0) * std::f64::consts::PI;
            let x = a * t * t.cos();
            let y = a * t * t.sin();
            triangulation.add_point(Point::new(x, y));
        }
        // We shouldn't need to retriangulate any more edges.
        assert_eq!(triangulation.retriangulate_all(), 0);
        valid_triangulation(&triangulation);
    }

    #[test]
    fn small_spiral_triangulation() {
        let mut triangulation = ConstrainedTriangulation::new();
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(0.45344999999999996, 0.261799);
        let p3 = Point::new(0.5235989999999999, 0.9068999999999999);
        let p4 = Point::new(0.408105, 0.235619);
        let p5 = Point::new(0.47123899999999996, 0.81621);
        let p6 = Point::new(0.0, 1.413717);
        triangulation.add_point(p1);
        triangulation.add_point(p2);
        triangulation.add_point(p3);
        assert!(has_edge(&triangulation, p1, p2));
        assert!(has_edge(&triangulation, p2, p3));
        assert!(has_edge(&triangulation, p3, p1));
        triangulation.add_point(p4);
        assert!(has_edge(&triangulation, p1, p4));
        assert!(has_edge(&triangulation, p4, p2));
        triangulation.add_point(p5);
        assert!(has_edge(&triangulation, p1, p5));
        assert!(has_edge(&triangulation, p5, p3));
        assert!(has_edge(&triangulation, p2, p3));
        assert!(has_edge(&triangulation, p2, p3));
        triangulation.add_point(p6);

        // We shouldn't need to retriangulate any more edges.
        assert_eq!(triangulation.retriangulate_all(), 0);
        valid_triangulation(&triangulation);
    }

    #[test]
    fn line_triangulation() {
        let mut triangulation = ConstrainedTriangulation::new();
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(5.0, 0.0);
        let p3 = Point::new(7.0, 0.0);
        triangulation.add_point(p1);
        triangulation.add_point(p2);
        triangulation.add_point(p3);
        assert_eq!(triangulation.qeds().unwrap().quads.len(), 14);
    }

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

    #[test]
    fn no_intersections() {
        let mut triangulation = ConstrainedTriangulation::new();
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(10.0, 0.0);
        let p3 = Point::new(5.0, 10.0);
        triangulation.add_point(p1);
        triangulation.add_point(p2);
        triangulation.add_point(p3);
        valid_triangulation(&triangulation);
        assert_eq!(
            triangulation.find_intersections(Point::new(2.0, 2.0), Point::new(4.0, 3.0)),
            vec![]
        );
    }

    #[test]
    fn simple_intersections() {
        let mut triangulation = ConstrainedTriangulation::new();
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(10.0, 0.0);
        let p3 = Point::new(5.0, 10.0);
        let p4 = Point::new(20.0, 20.0);
        triangulation.add_point(p1);
        triangulation.add_point(p2);
        triangulation.add_point(p3);
        triangulation.add_point(p4);
        valid_triangulation(&triangulation);
        {
            // single intersection
            let intersections =
                triangulation.find_intersections(Point::new(2.0, 2.0), Point::new(9.0, 8.0));
            assert_eq!(intersections.len(), 1);
            let first_intersection = unsafe { triangulation.qeds.edge_a_ref(intersections[0]) };
            // The edges found have a specified direction.
            assert_eq!(
                (
                    first_intersection.edge().point.point(),
                    first_intersection.sym().edge().point.point()
                ),
                (p2, p3)
            );
        }
        {
            // two intersections
            let intersections =
                triangulation.find_intersections(Point::new(2.0, 2.0), Point::new(18.0, 6.0));
            assert_eq!(intersections.len(), 2);
            let first_intersection = unsafe { triangulation.qeds.edge_a_ref(intersections[0]) };
            // The edges found have a specified direction.
            assert_eq!(
                (
                    first_intersection.edge().point.point(),
                    first_intersection.sym().edge().point.point()
                ),
                (p2, p3)
            );
            let second_intersection = unsafe { triangulation.qeds.edge_a_ref(intersections[1]) };
            assert_eq!(
                (
                    second_intersection.edge().point.point(),
                    second_intersection.sym().edge().point.point()
                ),
                (p2, p4)
            );
        }
        let p5 = Point::new(7.5, 5.0);
        triangulation.add_point(p5);
        {
            // through a vertex (but starting within a triangle).
            println!("Testing through a vertex");
            let intersections =
                triangulation.find_intersections(Point::new(6.0, 5.0), Point::new(9.0, 5.0));
            assert_eq!(intersections.len(), 0);
            // let first_intersection = unsafe {triangulation.qeds.edge_a_ref(intersections[0])};
            // // The edges found have a specified direction.
            // assert_eq!((first_intersection.edge().point.point(),first_intersection.sym().edge().point.point()),(p2,p3));
            // let second_intersection = unsafe {triangulation.qeds.edge_a_ref(intersections[1])};
            // assert_eq!((second_intersection.edge().point.point(),second_intersection.sym().edge().point.point()),(p2,p4));
        }
        {
            println!("Testing from a vertex (without leaving initial triangle)");
            let intersections =
                triangulation.find_intersections(Point::new(0.0, 0.0), Point::new(1.0, 1.0));
            assert_eq!(intersections.len(), 0);
        }
        {
            // through a vertex (but starting within a triangle).
            println!("Testing from a vertex (with 1 intersection)");
            let intersections =
                triangulation.find_intersections(Point::new(0.0, 0.0), Point::new(8.0, 10.0));
            assert_eq!(intersections.len(), 1);

            let first_intersection = unsafe { triangulation.qeds.edge_a_ref(intersections[0]) };
            assert_eq!(
                (
                    first_intersection.edge().point.point(),
                    first_intersection.sym().edge().point.point()
                ),
                (p5, p3)
            );
        }
    }

    #[test]
    fn simple_constraints() {
        let mut triangulation = ConstrainedTriangulation::new();
        let p1 = Point::new(-10.0, 0.0);
        let p2 = Point::new(10.0, 0.0);
        let p3 = Point::new(5.0, 10.0);
        triangulation.add_constraint(p1, p2);
        assert!(has_constraint(&triangulation, p1, p2));
        // triangulation.add_point(p2);
        // triangulation.add_point(p3);
        // valid_triangulation(&triangulation);
        // assert_eq!(
            // triangulation.find_intersections(Point::new(2.0, 2.0), Point::new(4.0, 3.0)),
            // vec![]
        // );
    }
}

/// Automatically clamped to the segment.
// fn intersection_point(s1: (Point,Point), line)

#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub enum IntersectionResult {
    Parallel,
    LineIntersection(Point),
}

fn segment_intersection(s1: (Point, Point), s2: (Point, Point)) -> IntersectionResult {
    let p = s1.0;
    let rAbs = s1.1;
    let q = s2.0;
    let sAbs = s2.1;
    // line1 t = p + t*r
    // line2 r = q + u*s
    let r = rAbs - p;
    let s = sAbs - q;

    // Sometimes we get numerical errors when the points coincide, so we
    // check for that first
    let t = if p == q {
        0.0
    } else if p == sAbs {
        0.0
    } else if rAbs == q {
        1.0
    } else if rAbs == sAbs {
        1.0
    } else {
        cross(q - p, s) / cross(r, s)
    };

    let u = if q == p {
        0.0
    } else if q == rAbs {
        0.0
    } else if sAbs == p {
        1.0
    } else if sAbs == rAbs {
        1.0
    } else {
        cross(q - p, r) / cross(r, s)
    };

    // If r `cross` s is 0 then the lines are parallel and do not intersect, so
    // return nothing.
    match cross(r, s) {
        x if x == 0.0 => IntersectionResult::Parallel,
        // x if !segments_overlap(s1,s2) && tIn && uIn = todo!(),
        _ => IntersectionResult::LineIntersection(Point::new(t, u)),
    }
}

fn cross(pa: Point, pb: Point) -> f64 {
    pa.x * pb.y - pb.x * pa.y
}
