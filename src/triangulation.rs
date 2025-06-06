use crate::point::*;
use crate::qeds::*;

pub trait HasPoint {
    fn point(&self) -> Point;
}

impl HasPoint for Point {
    fn point(&self) -> Point {
        *self
    }
}

impl<'a, AData: HasPoint, BData> EdgeRefA<'a, AData, BData> {
    pub fn midpoint(&self) -> Point {
        let this_point = self.edge().point.point();
        let next_point = self.sym().edge().point.point();
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
        let pa = self.edge().point.point();
        let pb = self.sym().edge().point.point();
        // We have special treatment of infinite edges.
        match left_or_right(pa, pb, point) {
            Direction::Right => Yes,
            Direction::Straight => On,
            Direction::Left => No,
        }
    }

    pub fn lies_right_strict(&self, point: Point) -> bool {
        self.lies_right(point) == Lies::Yes
    }

    pub fn lies_left(&self, point: Point) -> bool {
        let pa = self.edge().point.point();
        let pb = self.sym().edge().point.point();
        // We have special treatment of infinite edges.
        match left_or_right(pa, pb, point) {
            Direction::Right => false,
            Direction::Straight => false,
            Direction::Left => true,
        }
    }
}

impl<'a, AData: HasPoint, BData> Iterator for FaceVerticesIter<'a, AData, BData> {
    type Item = Point;
    fn next(&mut self) -> Option<Self::Item> {
        if self.next_index >= self.face.edges.len() {
            None
        } else {
            let edge = self.face.edges.get(self.next_index)?;
            self.next_index += 1;
            Some(edge.edge().point.point())
        }
    }
}

impl<AData: HasPoint + Clone, BData: Clone> Face<'_, AData, BData> {
    pub fn centroid(&self) -> Option<Point> {
        let signed_area = self.signed_area();
        if signed_area <= 0.0 {
            return None;
        }
        let mut xsum = 0.0;
        let mut ysum = 0.0;
        for (pa, pb) in self.vertices_pairs() {
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
            sum += pa.x * pb.y - pb.x * pa.y;
        }
        0.5 * sum
    }
    pub fn vertices(&self) -> FaceVerticesIter<'_, AData, BData> {
        FaceVerticesIter::new(self)
    }
    pub fn vertices_pairs(&self) -> VertexPairIter<'_, AData, BData> {
        self.vertices().zip(self.vertices().cycle().skip(1))
    }
}

#[derive(Copy, Clone, Debug, Ord, PartialOrd, Eq, PartialEq, Hash)]
pub enum Lies {
    Yes,
    No,
    On,
}

/// Determine whether a set of 4 points satisfies the Delaunay criterion. This
/// assumes that the pointes are sorted in a CCW fashion.
pub fn del_test_ccw(a: Point, b: Point, c: Point, d: Point) -> bool {
    robust::incircle(
        robust::Coord { x: a.x, y: a.y },
        robust::Coord { x: b.x, y: b.y },
        robust::Coord { x: c.x, y: c.y },
        robust::Coord { x: d.x, y: d.y },
    ) > 0.0
}

pub fn cocircular(a: Point, b: Point, c: Point, d: Point) -> bool {
    robust::incircle(
        robust::Coord { x: a.x, y: a.y },
        robust::Coord { x: b.x, y: b.y },
        robust::Coord { x: c.x, y: c.y },
        robust::Coord { x: d.x, y: d.y },
    ) == 0.0
}

pub fn is_ccw(a: Point, b: Point, c: Point) -> bool {
    robust::orient2d(a.into(), b.into(), c.into()) > 0.0
}

pub fn is_ccw_certain(a: Point, b: Point, c: Point) -> bool {
    robust::orient2d(a.into(), b.into(), c.into()) > 0.1
}

pub fn direction(a: Point, b: Point, c: Point) -> Direction {
    let d = robust::orient2d(a.into(), b.into(), c.into());
    if d > 0.0 {
        Direction::Left
    } else if d < 0.0 {
        Direction::Right
    } else {
        Direction::Straight
    }
}

#[derive(Clone, Debug)]
/// A Qeds data structure specialised to a 2d triangulation.
pub struct Triangulation {
    /// The quad-edge data structure we use as the basis for the triangulation.
    pub qeds: Qeds<Point, ()>,
    pub boundary_edge: EdgeTarget,
    pub bounds: Option<(Point, Point)>,
}

impl Triangulation {
    pub fn new() -> Self {
        let mut qeds = Qeds::new();
        let sw = Point::new(SafeFloat::MIN.0, SafeFloat::MIN.0);
        let se = Point::new(SafeFloat::MAX.0, SafeFloat::MIN.0);
        let ne = Point::new(SafeFloat::MAX.0, SafeFloat::MAX.0);
        let nw = Point::new(SafeFloat::MIN.0, SafeFloat::MAX.0);
        let south = qeds.make_edge_with_a(sw, se).target();
        let east = qeds.make_edge_with_a(se, ne).target();
        {
            qeds.splice(east, south.sym());
            let centre = qeds.connect(east, south).target();

            let north = qeds.make_edge_with_a(ne, nw).target();
            qeds.splice(north, east.sym());
            qeds.connect(north, centre.sym()).target();

            Triangulation {
                qeds,
                boundary_edge: south,
                bounds: None,
            }
        }
    }

    pub fn qeds(&self) -> Option<&Qeds<Point, ()>> {
        Some(&self.qeds)
    }

    pub fn some_edge_a(&self) -> Option<EdgeRefA<'_, Point, ()>> {
        let (i, _) = self.qeds.quads.iter().next()?;
        Some(self.qeds.edge_a_ref(EdgeTarget::new(i, 0, 0)))
    }

    // TODO: we need to be able to locate on an edge etc.
    // TODO: This has not yet been proved to be stable. It may also loop inifintely
    pub fn locate(&self, point: Point) -> Option<EdgeRefA<'_, Point, ()>> {
        use rand::Rng;
        let mut e = self.some_edge_a().unwrap();
        let mut rng = rand::rng();
        let mut current_iterations = 0;
        loop {
            current_iterations += 1;
            if current_iterations > 2000 {
                panic!("locating failed for: {}", point);
            }
            // A random variable to determine if onext is tested first. If not
            // it is tested second.
            let onext_first: bool = rng.random();
            if point == e.edge().point || point == e.sym().edge().point {
                return Some(e);
            } else if e.lies_right_strict(point) {
                e = e.sym();
            } else {
                #[allow(clippy::collapsible_else_if)]
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

    pub fn add_point(&mut self, mut point: Point) {
        point.x *= 1.0e6;
        point.x = point.x.round();
        point.x *= 1.0e-6;

        point.y *= 1.0e6;
        point.y = point.y.round();
        point.y *= 1.0e-6;

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
        let edge_target = self.locate(point).unwrap().target();
        self.add_to_l_face(edge_target, point);
    }

    fn add_to_l_face(&mut self, mut edge_target: EdgeTarget, point: Point) {
        {
            // println!("Edge Target: {} - {}", self.qeds.edge_a_ref(edge_target).edge().point,self.qeds.edge_a_ref(edge_target).sym().edge().point);
            if self.qeds.edge_a_ref(edge_target).edge().point == point
                || self.qeds.edge_a_ref(edge_target).sym().edge().point == point
            {
                return;
            } else if self.qeds.edge_a_ref(edge_target).lies_right(point) == Lies::On {
                // println!("Lies on edge: {} - {}", self.qeds.edge_a_ref(edge_target).edge().point,self.qeds.edge_a_ref(edge_target).sym().edge().point);
                let oprev = self.qeds.edge_a_ref(edge_target).oprev().target();
                self.qeds.delete(edge_target);
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
            // The suspect edges are e(.Onext.Lprev)^k for k=0,1,2,3...
            let mut e = edge_target;
            // println!("Start Swap Loop: first: {}", first);
            loop {
                let t = self.qeds.edge_a_ref(e).oprev().target();
                let t_dest = self.qeds.edge_a_ref(t).sym().edge().point;
                let e_dest = self.qeds.edge_a_ref(e).sym().edge().point;
                let e_org = self.qeds.edge_a_ref(e).edge().point;
                // print!("Inspecting Edge for swap: {} {} : [{},{},{},{}] : ", e_org, e_dest, e_org, t_dest, e_dest, point);
                if self.qeds.edge_a_ref(e).lies_right_strict(t_dest)
                    && del_test_ccw(e_org, t_dest, e_dest, point)
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
        }
    }

    pub fn add_point_no_swap(&mut self, mut point: Point) {
        point.x *= 1.0e6;
        point.x = point.x.round();
        point.x *= 1.0e-6;

        point.y *= 1.0e6;
        point.y = point.y.round();
        point.y *= 1.0e-6;

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
        let edge_target = self.locate(point).unwrap().target();
        self.add_to_l_face_no_swap(edge_target, point);
    }

    fn add_to_l_face_no_swap(&mut self, mut edge_target: EdgeTarget, point: Point) {
        {
            if self.qeds.edge_a_ref(edge_target).edge().point == point
                || self.qeds.edge_a_ref(edge_target).sym().edge().point == point
            {
                return;
            } else if self.qeds.edge_a_ref(edge_target).lies_right(point) == Lies::On {
                println!(
                    "Lies on edge: {} - {}",
                    self.qeds.edge_a_ref(edge_target).edge().point,
                    self.qeds.edge_a_ref(edge_target).sym().edge().point
                );
                let oprev = self.qeds.edge_a_ref(edge_target).oprev().target();
                self.qeds.delete(edge_target);
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
            // The suspect edges are e(.Onext.Lprev)^k for k=0,1,2,3...
            let mut e = edge_target;
            println!("Start Swap Loop");
            loop {
                let t = self.qeds.edge_a_ref(e).oprev().target();
                let t_dest = self.qeds.edge_a_ref(t).sym().edge().point;
                let e_dest = self.qeds.edge_a_ref(e).sym().edge().point;
                let e_org = self.qeds.edge_a_ref(e).edge().point;
                println!("Inspecting Edge for swap (no effect): {} {}", e_org, e_dest);
                if self.qeds.edge_a_ref(e).lies_right_strict(t_dest)
                    && del_test_ccw(e_org, t_dest, e_dest, point)
                {
                    // self.swap(e);
                    e = t;
                } else if e_org == first {
                    break;
                } else {
                    e = self.qeds.edge_a_ref(e).onext().l_prev().target();
                }
            }
        }
    }

    /// Assumes the EdgeTargets have already been put in the correct order.
    pub fn add_external_point_ordered(&mut self, boundary_edges: Vec<EdgeTarget>, point: Point) {
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

    pub fn add_external_point_unordered(
        &mut self,
        mut boundary_edges: Vec<EdgeTarget>,
        point: Point,
    ) {
        // Reorder the edges so that their base points are CCW around the point.
        boundary_edges.sort_by(|a, b| {
            let pa = self.qeds.edge_a_ref(*a).edge().point;
            let pb = self.qeds.edge_a_ref(*b).edge().point;
            match left_or_right(point, pb, pa) {
                Direction::Left => std::cmp::Ordering::Less,
                Direction::Straight => std::cmp::Ordering::Equal,
                Direction::Right => std::cmp::Ordering::Greater,
            }
        });
        self.add_external_point_ordered(boundary_edges, point);
    }

    pub fn boundary(&self) -> BoundaryIter<'_, Point, ()> {
        BoundaryIter::new(&self.qeds, self.boundary_edge)
    }

    pub fn add_external_point(
        &mut self,
        boundary_edge: EdgeTarget,
        p: Point,
    ) -> (EdgeRefA<'_, Point, ()>, EdgeRefA<'_, Point, ()>) {
        let dest = self.qeds.edge_a_ref(boundary_edge).sym().edge().point;
        let e = self.qeds.make_edge_with_a(p, dest).target();
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
        let a = edge.oprev().edge().point;
        let b = edge.r_prev().sym().edge().point;
        let c = edge.l_next().edge().point;
        let d = edge.onext().sym().edge().point;
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
        #[allow(clippy::needless_collect)]
        let edge_targets: Vec<EdgeTarget> =
            self.qeds.base_edges().map(|edge| edge.target()).collect();
        for e in edge_targets.into_iter() {
            unsafe {
                if self.del_test(e) {
                    swaps += 1;
                    self.swap(e);
                }
            }
        }
        swaps
    }

    pub fn swap(&mut self, e: EdgeTarget) {
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

impl Default for Triangulation {
    fn default() -> Self {
        Self::new()
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

pub fn has_edge(triangulation: &Triangulation, pa: Point, pb: Point) -> bool {
    get_edge(triangulation, pa, pb).is_some()
}

pub fn get_edge(
    triangulation: &Triangulation,
    pa: Point,
    pb: Point,
) -> Option<EdgeRefA<'_, Point, ()>> {
    for (i, quad) in triangulation.qeds().unwrap().quads.iter() {
        let edge1 = quad.edges_a[0].point;
        let edge2 = quad.edges_a[1].point;
        let edge_target = if edge1 == pa && edge2 == pb {
            EdgeTarget::new(i, 0, 0)
        } else if edge1 == pb && edge2 == pa {
            EdgeTarget::new(i, 2, 0)
        } else {
            continue;
        };
        return Some(triangulation.qeds().unwrap().edge_a_ref(edge_target));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_point_triangulation_only() {
        let mut triangulation = Triangulation::new();
        let p1 = Point::new(0.0, 0.0);
        triangulation.add_point(p1);
    }

    #[test]
    fn one_point_triangulation_location() {
        let mut triangulation = Triangulation::new();
        // assert_eq!(triangulation.qeds().unwrap().quads.len(), 5);
        let south = { triangulation.qeds.edge_a_ref(EdgeTarget::new(0, 0, 0)) };
        assert_eq!(
            south.edge().point,
            Point::new(SafeFloat::MIN.0, SafeFloat::MIN.0)
        );
        assert_eq!(
            south.sym().edge().point,
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

        let east = { triangulation.qeds.edge_a_ref(EdgeTarget::new(1, 0, 0)) };
        assert_eq!(
            east.edge().point,
            Point::new(SafeFloat::MAX.0, SafeFloat::MIN.0)
        );
        assert_eq!(
            east.sym().edge().point,
            Point::new(SafeFloat::MAX.0, SafeFloat::MAX.0)
        );

        let centre = { triangulation.qeds.edge_a_ref(EdgeTarget::new(2, 0, 0)) };
        assert_eq!(
            centre.edge().point,
            Point::new(SafeFloat::MAX.0, SafeFloat::MAX.0)
        );
        assert_eq!(
            centre.sym().edge().point,
            Point::new(SafeFloat::MIN.0, SafeFloat::MIN.0)
        );

        let north = { triangulation.qeds.edge_a_ref(EdgeTarget::new(3, 0, 0)) };
        assert_eq!(
            north.edge().point,
            Point::new(SafeFloat::MAX.0, SafeFloat::MAX.0)
        );
        assert_eq!(
            north.sym().edge().point,
            Point::new(SafeFloat::MIN.0, SafeFloat::MAX.0)
        );

        let west = { triangulation.qeds.edge_a_ref(EdgeTarget::new(4, 0, 0)) };
        assert_eq!(
            west.edge().point,
            Point::new(SafeFloat::MIN.0, SafeFloat::MAX.0)
        );
        assert_eq!(
            west.sym().edge().point,
            Point::new(SafeFloat::MIN.0, SafeFloat::MIN.0)
        );

        let p1 = Point::new(1.0, 0.0);
        triangulation.add_point(p1);
        assert_eq!(triangulation.qeds().unwrap().quads.len(), 8);
        assert_eq!(triangulation.retriangulate_all(), 0);
    }

    #[test]
    fn two_point_triangulation() {
        let mut triangulation = Triangulation::new();
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(1.0, 1.0);
        triangulation.add_point(p1);
        triangulation.add_point(p2);
        for (i, edge) in triangulation.qeds().unwrap().base_edges().enumerate() {
            println!(
                "Edge[{}]: {}-{}",
                i,
                edge.edge().point,
                edge.sym().edge().point
            );
        }
        assert_eq!(triangulation.qeds().unwrap().quads.len(), 11);
        assert!(has_edge(&triangulation, p1, p2), "missing edge");
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
        for (i, _quad) in triangulation.qeds().unwrap().quads.iter() {
            let target1 = EdgeTarget::new(i, 0, 0);
            let target2 = EdgeTarget::new(i, 2, 0);
            {
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
                        .point,
                    triangulation
                        .qeds()
                        .unwrap()
                        .edge_a_ref(target1)
                        .sym()
                        .edge()
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
                        .point,
                    triangulation
                        .qeds()
                        .unwrap()
                        .edge_a_ref(target2)
                        .sym()
                        .edge()
                        .point
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
        // let p4 = Point::new(11.0, 11.0);
        // let p5 = Point::new(15.0, 2.5);
        // let p6 = Point::new(12.0, 12.5);
        // let p7 = Point::new(15.0, 12.5);
        triangulation.add_point(p1);
        triangulation.add_point(p2);
        triangulation.add_point(p3);
        valid_triangulation(&triangulation);
        // triangulation.add_point(p4);
    }

    #[test]
    fn big_spiral_triangulation() {
        let mut triangulation = Triangulation::new();
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
        let mut triangulation = Triangulation::new();
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
        let mut triangulation = Triangulation::new();
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
        {
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
        {
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
            centre.x /= n as f64;
            centre.y /= n as f64;
            assert_eq!(centre, Point::new(2.5, 5.0 / 3.0))
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

        {
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

    #[test]
    fn quad2() {
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
        {
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
        {
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
            centre.x /= n as f64;
            centre.y /= n as f64;
            assert_eq!(centre, Point::new(2.5, 5.0 / 3.0))
        }
    }

    #[test]
    fn bad_collinear() {
        let b = Point::new(-10000.0, -1.0);
        let c = Point::new(-9000.0, -1.0);
        let d = Point::new(f64::MIN, -1.0 - f64::EPSILON);
        let p = Point::new(10000.0, -1.0);
        // This is a problem as it means we can have a triangle with both sides
        // collinear with another point. This indicates that our collinearity
        // test is insufficient.
        assert_ne!(b, d);
        assert_eq!(left_or_right(b, c, p), Direction::Straight);
        assert_eq!(left_or_right(d, c, p), Direction::Right);
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
        {
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
}
