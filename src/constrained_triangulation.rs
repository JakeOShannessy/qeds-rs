use crate::point::*;
use crate::qeds::*;
use nalgebra::Matrix4;
use std::{
    collections::{HashMap, HashSet, VecDeque},
    num::NonZeroUsize,
};

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
    fn is_constraint_or_point(self, triangulation: &ConstrainedTriangulation) -> bool {
        match self {
            Intersection::Edge(e) => e.edge(triangulation).point.constraint,
            Intersection::Point(_) => true,
        }
    }
}

/// An L3 Node Target is an EdgeTarget that is guranteed to be both a NodeTarget
/// and an L3 Node.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct L3NodeTarget(NodeTarget);

impl L3NodeTarget {
    pub fn new_unchecked(target: EdgeTarget) -> Self {
        L3NodeTarget(NodeTarget(target))
    }

    pub fn connections(self, full_map: &FullMap) -> (L3NodeTarget, L3NodeTarget, L3NodeTarget) {
        match full_map.classification.get(&self.into()).unwrap().level {
            LevelInfo::L3(Some(a), Some(b), Some(c)) => (a, b, c),
            _ => unreachable!(),
        }
    }

    /// Get all the neighbouring L3 nodes excluding itself. This means either 1
    /// or 3 nodes. It cannot be 0, as it would then no longer be an L3 node. It
    /// cannot be 2 as it cannot connect to itself exactly once.
    pub fn l3_neighbours(self, full_map: &FullMap) -> Vec<(u8, L3NodeTarget)> {
        let (a, b, c) = self.connections(full_map);
        let mut nodes = Vec::new();
        if a != self {
            nodes.push((0, a));
        }
        if b != self {
            nodes.push((1, b));
        }
        if c != self {
            nodes.push((2, c));
        }
        nodes
    }

    pub fn corridor<'a>(self, full_map: &'a FullMap, edge_number: u8) -> CorridorIter<'a> {
        let node: NodeTarget = self.into();
        let edge = node.get_edge(&full_map.triangulation, edge_number);
        let next = edge.triangle_across(&full_map.triangulation);
        CorridorIter {
            full_map,
            current: node,
            next: Some(next),
        }
    }
}

impl From<L3NodeTarget> for NodeTarget {
    fn from(nt: L3NodeTarget) -> Self {
        nt.0
    }
}

impl From<L3NodeTarget> for EdgeTarget {
    fn from(nt: L3NodeTarget) -> Self {
        (nt.0).0
    }
}

/// An L3 Node Target is an EdgeTarget that is guranteed to be both a NodeTarget
/// and an L3 Node.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct L2NodeTarget(NodeTarget);

impl L2NodeTarget {
    pub fn new_unchecked(target: EdgeTarget) -> Self {
        L2NodeTarget(NodeTarget(target))
    }

    pub fn as_node(self) -> NodeTarget {
        self.0
    }

    pub fn as_edge(self) -> EdgeTarget {
        (self.0).0
    }

    pub fn connections(self, full_map: &FullMap) -> Option<[L3NodeTarget; 2]> {
        match full_map.classification.get(&self.into()).unwrap().level {
            LevelInfo::L2(Some(a), Some(b)) => Some([a, b]),
            LevelInfo::L2(None, None) => None,
            _ => unreachable!(),
        }
    }

    // /// Get all the neighbouring L3 nodes excluding itself. This means either 1
    // /// or 3 nodes. It cannot be 0, as it would then no longer be an L3 node. It
    // /// cannot be 2 as it cannot connect to itself exactly once.
    // pub fn l3_neighbours(self, full_map: &FullMap) -> Vec<(u8, L3NodeTarget)> {
    //     let (a,b,c) = self.connections(full_map);
    //     let mut nodes = Vec::new();
    //     if a != self {
    //         nodes.push((0, a));
    //     }
    //     if b != self {
    //         nodes.push((1, b));
    //     }
    //     if c != self {
    //         nodes.push((2, c));
    //     }
    //     nodes
    // }

    // pub fn corridor<'a>(self, full_map: &'a FullMap, edge_number: u8) -> CorridorIter<'a> {
    //     let node: NodeTarget = self.into();
    //     let edge = node.get_edge(full_map, edge_number);
    //     let next = edge.triangle_across(&full_map.triangulation);
    //     CorridorIter {
    //         full_map,
    //         prev: node,
    //         next,
    //     }
    // }
}

impl From<L2NodeTarget> for NodeTarget {
    fn from(nt: L2NodeTarget) -> Self {
        nt.0
    }
}

impl From<L2NodeTarget> for EdgeTarget {
    fn from(nt: L2NodeTarget) -> Self {
        (nt.0).0
    }
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

    pub fn get_edge(self, triangulation: &ConstrainedTriangulation, i: u8) -> EdgeTarget {
        let edges = self.edges(triangulation);
        edges[i as usize]
    }

    pub fn edges(self, triangulation: &ConstrainedTriangulation) -> Vec<EdgeTarget> {
        let edge_ref = unsafe { triangulation.qeds.edge_a_ref(self.into()) };
        edge_ref
            .l_face()
            .edges
            .into_iter()
            .map(|e| e.target())
            .collect()
    }

    pub fn adjacent_tris(self, triangulation: &ConstrainedTriangulation) -> Vec<NodeTarget> {
        let edges = self.edges(triangulation);
        let mut adjacents = Vec::new();
        for edge in edges {
            let edge_data = edge.edge(triangulation);
            if !edge_data.point.constraint {
                let tri = edge.triangle_across(triangulation);
                adjacents.push(tri);
            }
        }
        adjacents
    }

    pub fn n_constrained_edges(self, triangulation: &ConstrainedTriangulation) -> usize {
        let edge_ref = unsafe { triangulation.qeds.edge_a_ref(self.into()) };
        edge_ref.n_constrained_edges()
    }

    // How many of the triangles adjacent to this triangle (across an
    // unconstrained edge) are of the give level?
    pub fn num_adjacent_level(
        self,
        triangulation: &ConstrainedTriangulation,
        tri_info: &HashMap<NodeTarget, TriInfo>,
        level: Level,
    ) -> usize {
        let edge_ref = unsafe { triangulation.qeds.edge_a_ref(self.into()) };
        edge_ref.num_adjacent_level(tri_info, level)
    }
}

impl From<NodeTarget> for EdgeTarget {
    fn from(nt: NodeTarget) -> Self {
        nt.0
    }
}

pub struct NodeWalker<'a> {
    full_map: &'a FullMap,
    // current_node: NodeTarget,
    unprocessed_stack: Vec<NodeTarget>,
    processed_stack: HashSet<NodeTarget>,
}

impl<'a> Iterator for NodeWalker<'a> {
    type Item = NodeTarget;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(current_node) = self.unprocessed_stack.pop() {
            let adjacent = current_node.adjacent_tris(&self.full_map.triangulation);
            for tri in adjacent {
                if !self.processed_stack.contains(&tri) {
                    self.unprocessed_stack.push(tri);
                }
            }
            self.processed_stack.insert(current_node);
            Some(current_node)
        } else {
            None
        }
    }
}

/// An [`L3Path`] is a path from one L3 node to another L3 node via other L3
/// nodes. This path is inclusive of the start and end nodes.
#[derive(Clone, Debug, Ord, PartialOrd, Eq, PartialEq, Hash)]
pub struct L3Path {
    start_node: L3NodeTarget,
    path: Vec<(u8, L3NodeTarget)>,
}

impl L3Path {
    /// Create a new path which consists solely of the single start node.
    pub fn new(start_node: L3NodeTarget) -> Self {
        Self {
            start_node,
            path: vec![],
        }
    }

    /// The start node.
    pub fn first(&self) -> L3NodeTarget {
        self.start_node
    }

    /// Get the end node of the path. This always returns a value as a path
    /// must have a start node at a minimum.
    pub fn last(&self) -> L3NodeTarget {
        if let Some((_, node)) = self.path.last() {
            *node
        } else {
            self.start_node
        }
    }

    pub fn contains(&self, node: L3NodeTarget) -> bool {
        if node == self.start_node {
            return true;
        }
        self.path.iter().find(|(_, x)| *x == node).is_some()
    }

    pub fn push(&mut self, path_element: (u8, L3NodeTarget)) {
        self.path.push(path_element);
    }

    /// Connect all the nodes in this [`L3Path`] to create a path of nodes.
    pub fn connect(&self, full_map: &FullMap) -> Vec<NodeTarget> {
        self.node_iter(full_map).collect()
    }

    pub fn node_iter<'a>(&'a self, full_map: &'a FullMap) -> L3PathNodeIter<'a> {
        L3PathNodeIter {
            full_map,
            current_node: self.start_node,
            pending_nodes: Vec::new(),
            path: self.path.as_slice(),
            ended: false,
        }
    }
}

pub struct L3PathNodeIter<'a> {
    full_map: &'a FullMap,
    current_node: L3NodeTarget,
    pending_nodes: Vec<NodeTarget>,
    path: &'a [(u8, L3NodeTarget)],
    ended: bool,
}

impl<'a> Iterator for L3PathNodeIter<'a> {
    type Item = NodeTarget;
    fn next(&mut self) -> Option<Self::Item> {
        // If we already have the next node in a corridor pending in the queue
        // pop and return that.
        if let Some(node) = self.pending_nodes.pop() {
            Some(node)
        } else {
            // Otherwise we must refill the pending nodes. First take the edge
            // number we want to step across from the currrent node.
            if let Some((edge_number, next_node)) = self.path.first() {
                let current_node: NodeTarget = self.current_node.into();
                let corridor = self.current_node.corridor(self.full_map, *edge_number);
                self.current_node = *next_node;
                let mut nodes: Vec<NodeTarget> = corridor.collect();
                // Reverse the nodes so we can use it as a stack.
                nodes.reverse();
                self.pending_nodes = nodes;
                self.path = &self.path[1..];
                Some(current_node)
            } else {
                // If there are no more elements on the path we have finished
                // iterating.
                if self.ended {
                    None
                } else {
                    self.ended = true;
                    Some(self.current_node.into())
                }
            }
        }
    }
}

pub type NodePath = Vec<NodeTarget>;

/// An iterator that walks along a corridor in one direction.
pub struct CorridorIter<'a> {
    full_map: &'a FullMap,
    current: NodeTarget,
    next: Option<NodeTarget>,
}

impl<'a> CorridorIter<'a> {
    pub fn new(
        full_map: &'a FullMap,
        start_node: L2NodeTarget,
    ) -> (CorridorIter<'a>, CorridorIter<'a>) {
        // base_edges is the edges of the start node (also known as the base).
        let base_edges: Vec<EdgeTarget> = start_node
            .as_node()
            .edges(&full_map.triangulation)
            .into_iter()
            .filter(|e| {
                let edge = e.edge(&full_map.triangulation);
                let is_constraint = edge.point.constraint;
                // If it is a constraint we need to short-circuit here as
                // tri_info might not be available.
                if is_constraint {
                    return false;
                }
                let tri_across = e.triangle_across(&full_map.triangulation);
                let tri_across_info = full_map
                    .classification
                    .get(&tri_across)
                    .expect(&format!("could not get tri_info for {:?}", tri_across));
                let is_l1_across = tri_across_info.level.as_level() == Level::L1;
                // We want to ignore any adges that are constraints or that have
                // an L1 node across them.
                !is_l1_across
            })
            .collect();
        if base_edges.len() != 2 {
            panic!("there should only be 2 edges");
        }
        let next_1 = base_edges[0].triangle_across(&full_map.triangulation);
        let next_2 = base_edges[1].triangle_across(&full_map.triangulation);
        // The left value of the tuple should correlate to the left link value
        // in TriInfo.

        (
            CorridorIter {
                full_map,
                current: start_node.into(),
                next: Some(next_1),
            },
            CorridorIter {
                full_map,
                current: start_node.into(),
                next: Some(next_2),
            },
        )
    }
}

impl<'a> Iterator for CorridorIter<'a> {
    type Item = NodeTarget;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(next) = self.next {
            // Current is what we will return.
            let current = self.current;
            self.current = next;
            if self
                .full_map
                .classification
                .get(&next)
                .unwrap()
                .level
                .as_level()
                == Level::L3
            {
                self.next = None;
            } else {
                let edges = next.edges(&self.full_map.triangulation);
                let next_nodes: Vec<NodeTarget> = edges
                    .into_iter()
                    .filter(|e| {
                        let tri_across = e.triangle_across(&self.full_map.triangulation);
                        // The edge is not a constraint.
                        let edge = e.edge(&self.full_map.triangulation);
                        !edge.point.constraint
                        // The triangle across is not the triangle we just visited.
                        && tri_across != current.into()
                        // There is not an L1 (or unknown) node across the edge. (TODO: check that this is performed correctly)
                        && self.full_map.classification
                            .get(&tri_across)
                            .map(|x| x.level.as_level() != Level::L1)
                            .unwrap_or(true)
                    })
                    .map(|e| e.triangle_across(&self.full_map.triangulation))
                    .collect();
                if next_nodes.len() != 1 {
                    panic!("(n1) there should only be 1 edges: {:?}", next_nodes);
                }
                self.next = Some(next_nodes[0]);
            }
            Some(current)
        } else {
            None
        }
    }
}

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
/// Contains trianglulation as well as higher level data.
pub struct FullMap {
    /// The quad-edge data structure we use as the basis for the triangulation.
    pub triangulation: ConstrainedTriangulation,
    pub classification: LinkageMap,
}

impl FullMap {
    pub fn new() -> Self {
        Self {
            triangulation: ConstrainedTriangulation::new(),
            classification: LinkageMap::new(),
        }
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
            .make_edge_with_a(Segment::new_constraint(se), Segment::new_constraint(sw))
            .sym()
            .target();
        let east = qeds
            .make_edge_with_a(Segment::new_constraint(se), Segment::new_constraint(ne))
            .target();
        let north = qeds
            .make_edge_with_a(Segment::new_constraint(ne), Segment::new_constraint(nw))
            .target();
        unsafe {
            qeds.splice(east, south.sym());
            let centre = qeds.connect(east, south).target();
            qeds.edge_a_mut(centre).point.constraint = false;
            qeds.edge_a_mut(centre.sym()).point.constraint = false;

            qeds.splice(north, east.sym());
            let west = qeds.connect(north, centre.sym()).target();
            qeds.edge_a_mut(west).point.constraint = true;
            qeds.edge_a_mut(west.sym()).point.constraint = true;
            ConstrainedTriangulation {
                qeds,
                boundary_edge: south,
                bounds: None,
            }
        }
    }

    pub fn triangles(&self) -> TriangleIter {
        TriangleIter::new(self)
    }
    pub fn nodes(&self) -> NodeIter {
        NodeIter::new(self)
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
        let edge = loop {
            current_iterations += 1;
            if current_iterations > 2000 {
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
    pub fn locate_tri(&self, point: Point) -> Option<EdgeRefA<Segment, ()>> {
        self.locate(point).map(|edge| edge.get_tri_canonical())
    }

    /// The edge this returns should always have the added point at its origin.
    pub fn add_point(&mut self, mut point: Point) -> Option<EdgeTarget> {
        point.snap();
        self.update_bounds(point);
        if let Some(edge_ref) = self.locate(point) {
            let edge_target = edge_ref.target();
            Some(self.add_to_l_face(edge_target, point))
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

    /// TODO: remove the recursion from this function as we don't have tail
    /// recursions. Edges are returned with each edge pointing from right to
    /// left across the intersection line. Second return value is an edge from
    /// the final point. It is necessary that the two points are already in the
    /// trianglulation as it is required for the end condition. If we pass
    /// through a point on a constraint segment, we also report that as an
    /// intersection. TODO: this needs extensive revision.
    fn find_intersections_between_points(&self, a: Point, b: Point) -> IntersectionIter {
        IntersectionIter::new(self, a, b)
    }

    fn find_next_constraint_intersection(&self, a: Point, b: Point) -> Option<Intersection> {
        IntersectionIter::new(self, a, b)
            .filter(|intersection| intersection.is_constraint_or_point(self))
            .next()
    }

    pub fn add_constraint(&mut self, mut pa: Point, mut pb: Point) -> Option<()> {
        let mut pa_edge = unsafe {
            pb = {
                let pb_edge = self.add_point(pb)?;
                self.qeds.edge_a_ref(pb_edge).edge().point.point()
            };
            // by virtue of using add_point, we might invalidate pb_edge.
            let pa_edge = self.add_point(pa)?;

            pa = self.qeds.edge_a_ref(pa_edge).edge().point.point();
            pa_edge
        };

        if pa == pb {
            return Some(());
        }
        // debug_assert_eq!(pa, pa_edge.edge(self).point.point());

        self.update_bounds(pa);
        self.update_bounds(pb);
        // The end points have now been added the and triangulation updated. The
        // next step is to add the edge constraint.
        //
        // Step 1. Wherever the new segment crosses an existing constrained
        // segment add a new point. Also update the edges accordingly. As part
        // of this process, the constraint will be split into a series of
        // smaller constraints.

        // The purpose of this loop is to keep breaking small chunks off the
        // constraint until we are done. Each iteration of this loop is the
        // application of one of these "constraint chunks".

        // TODO: don't consume this all at once, as indices may change.
        // 'main_loop: loop {
        // println!("Intersections from {} to {}", pa, pb);
        let mut iter_count = 0;
        loop {
            if pa == pb {
                break;
            }
            // println!("loop count {}", iter_count);
            if iter_count > 300 {
                panic!("exceeded iter count");
            } else {
                iter_count += 1;
            }

            // println!("constraint intersection found");
            // debug_assert_eq!(pa, pa_edge.edge(self).point.point());
            // In the course of this block, pa_edge gets changed. This makes
            // sense as we are inserting/deleting.
            let (px, px_edge) = {
                match self.find_next_constraint_intersection(pa, pb) {
                    // None => pb, // how do we get pb_edge
                    None => panic!("We require a point at pb"),
                    Some(intersection) => match intersection {
                        // If the intersection is a point, it already exists and we
                        // don't need to do anything.
                        Intersection::Point(p) => {
                            // println!("Intersection at Point: {}", p.edge(self).point.point());
                            (p.edge(self).point.point(), p)
                        }
                        // If the intersection is an edge, we need to add a point to
                        // that edge.
                        Intersection::Edge(edge_target) => {
                            // println!(
                            //     "Intersection at Edge: {}-{}",
                            //     edge_target.edge(self).point.point(),
                            //     edge_target.sym().edge(self).point.point()
                            // );
                            let edge_data = edge_target.edge(self);
                            let edge_sym_data = edge_target.sym().edge(self);
                            let ex = unsafe { self.qeds.edge_a_ref(edge_target) };
                            let ex_sym = unsafe { self.qeds.edge_a_ref(edge_target) }.sym();
                            let e = ex.edge();
                            let e_sym = ex_sym.edge();
                            if !edge_data.point.constraint {
                                // If point is not constrained we don't care about
                                // and continue to loop.
                                continue;
                            } else {
                                if !edge_sym_data.point.constraint {
                                    panic!("Inconsistent edge");
                                }
                                // Otherwise we need to insert this intersection
                                // point.
                                // First find/choose the intersection point.
                                let intersection_point: Point = match segment_intersection(
                                    (e.point.point(), e_sym.point.point()),
                                    (pa, pb),
                                ) {
                                    IntersectionResult::Parallel => panic!("unexpected parallel"),
                                    IntersectionResult::LineIntersection(t, _u) => {
                                        let p = e.point.point();
                                        let r = e_sym.point.point() - p;
                                        Point::new(p.x + t * r.x, p.y + t * r.y)
                                    }
                                };
                                // println!("Intersection Point: {}", intersection_point);
                                // Then insert the point.
                                // assert_eq!(px, unsafe{self.qeds.edge_a_ref(px_edge).edge().point.point()});
                                // inserted_edge_target is an edge leading away from
                                // the inserted target (which will become px).
                                let (inserted_point, inserted_edge_target) = unsafe {
                                    // Major changes occur here, pa_edge could even be deleted.
                                    drop(e);
                                    let target =
                                        self.add_point_to_edge(edge_target, intersection_point);
                                    let inserted_point =
                                        self.qeds.edge_a_ref(target).edge().point.point();
                                    (inserted_point, target)
                                };
                                // TODO: Do we need to update pa_edge here,
                                // possibly? Apparently so, but not always. It's
                                // likely being changed due to the swapping that
                                // occurs in [`add_point_to_edge`]. It's not
                                // immediately clear how to find a new one (at least
                                // in a reliable provable way).
                                // assert_eq!(pa, unsafe{self.qeds.edge_a_ref(pa_edge).edge().point.point()});
                                (inserted_point, inserted_edge_target)
                                // break inserted_point;
                            }
                        }
                    },
                }
            };
            // println!("found or made intersection");
            // debug_assert_eq!(pa, pa_edge.edge(self).point.point());
            // We now know that there are only unconstrained edges between pa
            // and px (possible aligned). There can be at mose one aligned edge,
            // this will be retuened to us as point intersection. This can only
            // not equal px if we have unconstrained points, which we should
            // not.

            // First find all of the intersecting unconstrained edges and put
            // them in a vector. TODO: We should understand this double looping
            // a little more.
            let mut iter_count_swap = 0;
            loop {
                if iter_count_swap > 300 {
                    panic!("Swap iterations exceeded");
                } else {
                    iter_count_swap += 1;
                }
                let mut intersecting_soft: Vec<EdgeTarget> = self
                    .find_intersections_between_points(pa, px)
                    .filter_map(|intersection| match intersection {
                        Intersection::Point(_) => None,
                        Intersection::Edge(e) => {
                            if !e.edge(self).point.constraint {
                                Some(e)
                            } else {
                                None
                            }
                        }
                    })
                    .collect();
                if intersecting_soft.len() == 0 {
                    break;
                }
                let mut i = 0;
                let mut iter_count_swap_inner = 0;
                loop {
                    // println!("intersecting(soft): {}", intersecting_soft.len());
                    if iter_count_swap_inner > 300 {
                        // panic!("Swap iterations (inner) exceeded");
                        // println!("Swap iterations (inner) exceeded");
                        break;
                    } else {
                        iter_count_swap_inner += 1;
                    }
                    if intersecting_soft.len() == 0 {
                        break;
                    }
                    if i >= intersecting_soft.len() {
                        break;
                    }
                    let edge = intersecting_soft[i];
                    if unsafe { self.concave_test(edge) } {
                        unsafe { self.swap(edge) };
                        intersecting_soft.swap_remove(i);
                        i = 0;
                    } else {
                        // println!("could not swap {}, will swap {}", i, i + 1);
                        // debug_assert_ne!(intersecting_soft.len(), 1, "could not swap all");
                        i += 1;
                    }
                }
            }
            debug_assert_eq!(
                self.find_intersections_between_points(pa, px)
                    .filter_map(|intersection| match intersection {
                        Intersection::Point(_) => None,
                        Intersection::Edge(e) => {
                            if !e.edge(self).point.constraint {
                                Some(e)
                            } else {
                                None
                            }
                        }
                    })
                    .count(),
                0
            );
            // At this point everything should be aligned, so we just need to
            // rotate around and add the constraint.
            // println!("edges aligned, adding constraint from {} to {}", pa, px);

            unsafe {
                debug_assert_eq!(px, px_edge.edge(self).point.point());
                let mut edge = px_edge;
                let initial_edge = edge;
                loop {
                    let first_point = self.qeds.edge_a_ref(edge).edge().point.point();
                    let other_point = self.qeds.edge_a_ref(edge).sym().edge().point.point();
                    // println!(
                    //     "pa: {}, first_point: {}, other_point: {}",
                    //     pa, first_point, other_point
                    // );
                    debug_assert_eq!(px, first_point);
                    if other_point == pa {
                        self.qeds.edge_a_mut(edge).point.constraint = true;
                        self.qeds.edge_a_mut(edge.sym()).point.constraint = true;
                        pa = px;
                        // pa_edge = edge.sym();
                        // debug_assert_eq!(pa, pa_edge.edge(self).point.point());
                        break;
                    } else {
                        edge = self.qeds.edge_a_ref(edge).onext().target();
                        if edge == initial_edge {
                            panic!("looped around ring (add_constraint)");
                            // println!("looped around ring (add_constraint)");
                            // return Some(());
                        }
                    }
                }
            }
            // println!("done constraint piece");
            // debug_assert_eq!(pa, pa_edge.edge(self).point.point());
        }
        Some(())
    }

    fn connect(&mut self, a: EdgeTarget, b: EdgeTarget) -> EdgeRefA<Segment, ()> {
        let e = self.qeds.connect(a, b).target();
        unsafe {
            self.qeds.edge_a_mut(e).point.constraint = false;
            self.qeds.edge_a_mut(e.sym()).point.constraint = false;
            self.qeds.edge_a_ref(e)
        }
    }

    /// This function accounts for the point lying on an existing edge or point.
    fn add_to_l_face(&mut self, edge_target: EdgeTarget, point: Point) -> EdgeTarget {
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
                self.add_point_to_edge_unchecked(edge_target, point)
            } else {
                self.add_to_quad_unchecked(edge_target, point)
            }
        }
    }

    /// Add a point to a specified edge. If the point lies on one of the
    /// vertices just add it there.
    fn add_point_to_edge(&mut self, edge_target: EdgeTarget, point: Point) -> EdgeTarget {
        unsafe {
            let point_a = self.qeds.edge_a_ref(edge_target).edge().point.point;
            let point_b = self.qeds.edge_a_ref(edge_target).sym().edge().point.point;

            if point_a == point {
                edge_target
            } else if point_b == point {
                edge_target.sym()
            } else {
                self.add_point_to_edge_unchecked(edge_target, point)
            }
        }
    }

    /// Same as [`add_point_to_edge`] but does not check if the point is on one
    /// of the vertices of the edge. Will also restore constraints where
    /// necessary.
    fn add_point_to_edge_unchecked(
        &mut self,
        mut edge_target: EdgeTarget,
        point: Point,
    ) -> EdgeTarget {
        unsafe {
            let point_a = self.qeds.edge_a_ref(edge_target).edge().point.point;
            let point_b = self.qeds.edge_a_ref(edge_target).sym().edge().point.point;

            let reinstate_as_constraint = self.qeds.edge_a_ref(edge_target).edge().point.constraint;
            {
                // We need to remember if this edge was constrained so
                // that we can reinstate it.
                let oprev_ref = self.qeds.edge_a_ref(edge_target).oprev();
                let oprev = oprev_ref.target();
                drop(oprev_ref);

                self.qeds.delete(edge_target);
                edge_target = oprev;
            }
            let return_value = self.add_to_quad_unchecked(edge_target, point);
            // Now we need to restore the constraint status of the edges we
            // deleted. The [`return_vale`] is one of the edges emenating from
            // this new point, so if we just iterate around that point we should
            // be able to find the edges we want to restore.
            if reinstate_as_constraint {
                let mut edge = return_value;
                loop {
                    let edge_ref = self.qeds.edge_a_ref(edge);
                    // Does the edge we're looking at have the correct points?
                    let org = edge_ref.edge().point.point();
                    let dest = edge_ref.sym().edge().point.point();
                    let has_points = (org == point_a || org == point_b || org == point)
                        && (dest == point_a || dest == point_b || dest == point);
                    if has_points {
                        // TODO: this constraint setting code is not of great design
                        // and sorely needs checking.
                        self.qeds.edge_a_mut(edge).point.constraint = true;
                        self.qeds.edge_a_mut(edge.sym()).point.constraint = true;
                    }
                    edge = self.qeds.edge_a_ref(edge).onext().target();
                    if edge == return_value {
                        break;
                    }
                }
            }
            return_value
        }
    }

    fn add_to_quad_unchecked(&mut self, mut edge_target: EdgeTarget, point: Point) -> EdgeTarget {
        unsafe {
            let first = self.qeds.edge_a_ref(edge_target).edge().point.point;
            let mut base = self
                .qeds
                .make_edge_with_a(Segment::new(first), Segment::new(point))
                .target();
            let return_value = base.sym();
            self.qeds.splice(base, edge_target);
            loop {
                let base_ref = self.connect(edge_target, base.sym());
                edge_target = base_ref.oprev().target();
                base = base_ref.target();
                if self.qeds.edge_a(edge_target.sym()).point.point == first {
                    break;
                }
            }
            edge_target = self.qeds.edge_a_ref(base).oprev().target();
            let e = edge_target;
            self.retriangulate_suspect_edges(e, point, first);
            debug_assert_eq!(
                self.qeds.edge_a_ref(return_value).edge().point.point(),
                point
            );
            return_value
        }
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
            if self.qeds.edge_a_ref(e).edge().point.constraint
                != self.qeds.edge_a_ref(e).edge().point.constraint
            {
                panic!("Error, inconsistent edge");
            }
            // TODO: we need to be cautious of infinite loops now that we're constrained.
            if self.qeds.edge_a_ref(e).lies_right_strict(t_dest)
                && del_test_ccw(e_org, t_dest, e_dest, point)
                && (self.qeds.edge_a_ref(e).edge().point.constraint == false)
                && (self.qeds.edge_a_ref(e).sym().edge().point.constraint == false)
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

    pub fn boundary(&self) -> BoundaryIter<Segment, ()> {
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
                if self.del_test(e)
                    && (self.qeds.edge_a_ref(e).edge().point.constraint == false)
                    && self.concave_test(e)
                {
                    swaps += 1;
                    self.swap(e);
                }
            }
        }
        swaps
    }

    pub unsafe fn swap(&mut self, e: EdgeTarget) {
        let a_is_constrained = self.qeds.edge_a_ref(e).edge().point.constraint;
        let b_is_constrained = self.qeds.edge_a_ref(e).sym().edge().point.constraint;
        let a = self.qeds.edge_a_ref(e).oprev().target();
        let b = self.qeds.edge_a_ref(e).sym().oprev().target();

        self.qeds.splice(e, a);
        self.qeds.splice(e.sym(), b);

        let a_lnext = self.qeds.edge_a_ref(a).l_next().target();
        self.qeds.splice(e, a_lnext);

        let b_lnext = self.qeds.edge_a_ref(b).l_next().target();
        self.qeds.splice(e.sym(), b_lnext);

        let mut a_dest = self.qeds.edge_a_ref(a).sym().edge().point;
        a_dest.constraint = a_is_constrained;
        let mut b_dest = self.qeds.edge_a_ref(b).sym().edge().point;
        b_dest.constraint = b_is_constrained;
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

    pub fn classify_triangles(&self) -> LinkageMap {
        let mut component: NonZeroUsize = NonZeroUsize::new(1).unwrap();
        // Create a map which contains all of our triangle information. The key
        // is the EdgeTarget.
        let mut tri_info: HashMap<NodeTarget, TriInfo> = HashMap::new();
        let mut queue: VecDeque<NodeTarget> = self.abstract_01(&mut tri_info, &mut component);
        // At this point all L0 and L1 nodes have been identified. Next we will
        // identify all L3 nodes. This does not consume the queue
        for triangle in queue.iter() {
            let n = triangle.n_constrained_edges(self);
            let m = triangle.num_adjacent_level(self, &tri_info, Level::L1);
            let this_tri_info = tri_info.get(&triangle);
            if n + m == 0 && this_tri_info.is_none() {
                self.abstract_3(&mut tri_info, *triangle, &mut component);
                component = NonZeroUsize::new(component.get() + 1).unwrap();
            }
        }
        // At this point we have found all the L0, L1, and L3 nodes. By
        // exclusion we have also identified all L2 nodes, but they have not yet
        // been marked. Any unmarked node after this point is L2.
        while let Some(triangle) = queue.pop_front() {
            if tri_info.get(&triangle).is_none() {
                let mut triangle_current = Some(triangle);
                while triangle_current.is_some() {
                    tri_info.insert(
                        triangle_current.unwrap(),
                        TriInfo {
                            level: LevelInfo::L2(None, None),
                            component: Some(component),
                        },
                    );
                    let mut triangle_next = None;
                    // assert!(triangle_current.unwrap().is_tri_canonical());
                    let tri_edges = triangle_current.unwrap().edges(self);
                    for edge in tri_edges {
                        let edge_data = edge.edge(self);
                        let triangle_temp = edge.triangle_across(self);
                        if edge_data.point.constraint
                            || tri_info
                                .get(&triangle_temp)
                                .map(|info| info.level.as_level() == Level::L1)
                                .unwrap_or(false)
                        {
                            if !edge_data.point.constraint {
                                self.collapse_rooted_tree(
                                    &mut tri_info,
                                    triangle_current.unwrap(),
                                    triangle_temp,
                                );
                            }
                        } else {
                            if tri_info.get(&triangle_temp).is_none() {
                                triangle_next = Some(triangle_temp);
                            }
                        }
                    }
                    triangle_current = triangle_next;
                }
                component = NonZeroUsize::new(component.get() + 1).unwrap();
            }
        }
        LinkageMap::from(tri_info)
    }

    // It must be the case that r is L2 and t and all its children are L1
    fn old_collapse_rooted_tree(
        &self,
        tri_info: &mut HashMap<NodeTarget, TriInfo>,
        r: NodeTarget,
        t: NodeTarget,
    ) {
        assert!(tri_info
            .get(&r)
            .map(|info| info.level.as_level() == Level::L2)
            .unwrap_or(false));
        assert!(tri_info
            .get(&t)
            .map(|info| info.level.as_level() == Level::L1)
            .unwrap_or(false));
        let component = tri_info.get(&r).unwrap().component.unwrap();
        let mut s = Vec::new();
        s.push(t);
        let mut a = Vec::new();
        a.push(0);
        while let Some(triangle_current) = s.pop() {
            {
                let tri_info_mut = tri_info.get_mut(&triangle_current).unwrap();
                // Set the component of the current triangle
                tri_info_mut.component = Some(component);
                // Get a mutable reference to Level 1 adjacency info
                if let LevelInfo::L1(ref mut s) = tri_info_mut.level {
                    *s = Some(r);
                } else {
                    panic!("incorrect level, should only be considering L1 triangles")
                }
            }
            for edge in triangle_current.edges(self) {
                let edge_ref = unsafe { self.qeds.edge_a_ref(edge) };
                let edge_data = edge.edge(self);
                let triangle_last = edge.triangle_across(self);
                if edge_data.point.constraint {
                    continue;
                }
                if tri_info
                    .get(&triangle_last)
                    .and_then(|x| x.component)
                    .is_none()
                    && !edge_data.point.constraint
                {
                    if triangle_last == r {
                        // TODO
                    } else {
                        // TODO
                    }
                    // TODO: surely we should only be pushing those with no
                    // constraint as an edge.
                    let e_right = edge_ref.l_next().target();
                    let e_right_edge_data = e_right.edge(self);
                    if !e_right_edge_data.point.constraint
                        && tri_info
                            .get(&e_right.triangle_across(self))
                            .and_then(|x| x.component)
                            .is_none()
                    {
                        s.push(e_right.triangle_across(self));
                        if tri_info.get_mut(&e_right.triangle_across(self)).is_none() {
                            panic!("pushed node (right) with no level to stack")
                        }
                    }
                    // a.push
                    let e_left = edge_ref.l_next().l_next().target();
                    let e_left_edge_data = e_left.edge(self);
                    if !e_left_edge_data.point.constraint
                        && tri_info
                            .get(&e_left.triangle_across(self))
                            .and_then(|x| x.component)
                            .is_none()
                    {
                        s.push(e_left.triangle_across(self));
                        if tri_info.get_mut(&e_left.triangle_across(self)).is_none() {
                            panic!("pushed node (left) with no level to stack")
                        }
                    }
                // a.push
                } else {
                    // {
                    //     let tri_info_mut = tri_info.get_mut(&triangle_current.target()).unwrap();
                    //     // Get a mutable reference to Level 1 adjacency info
                    //     if let LevelInfo::L1(ref mut s) = tri_info_mut.level {
                    //         // Set it to None
                    //         *s = None;
                    //     } else {
                    //         panic!("incorrect level, should only be considering L1 triangles")
                    //     }
                    // }
                    if let Some(_this_tri_info) = tri_info.get(&triangle_last) {
                    } else {
                    }
                    if edge_data.point.constraint {
                        // TODO
                    } else {
                        // TODO
                    }
                }
            }
        }
    }

    // It must be the case that r is L2 and t and all its children are L1
    fn collapse_rooted_tree(
        &self,
        tri_info: &mut HashMap<NodeTarget, TriInfo>,
        r: NodeTarget,
        t: NodeTarget,
    ) {
        debug_assert!(tri_info
            .get(&r)
            .map(|info| info.level.as_level() == Level::L2)
            .unwrap_or(false));
        debug_assert!(tri_info
            .get(&t)
            .map(|info| info.level.as_level() == Level::L1)
            .unwrap_or(false));
        debug_assert!(r.adjacent_tris(self).contains(&t));
        let component = tri_info.get(&r).unwrap().component.unwrap();
        let mut s = Vec::new();
        s.push(t);
        let mut a = Vec::new();
        a.push(0);
        while let Some(triangle_current) = s.pop() {
            {
                // Get a mutable reference to the tri_info for the current
                // triangle (initially the first L1 node).
                let tri_info_mut = tri_info.get_mut(&triangle_current).unwrap();
                // Set the component of the current triangle
                tri_info_mut.component = Some(component);
                // Get a mutable reference to Level 1 adjacency info
                if let LevelInfo::L1(ref mut s) = tri_info_mut.level {
                    *s = Some(r);
                } else {
                    panic!("incorrect level, should only be considering L1 triangles")
                }
            }
            // Iterate through all the edges of the current triangle.
            for edge in triangle_current.edges(self) {
                let edge_data = edge.edge(self);
                if edge_data.point.constraint {
                    continue;
                }
                let triangle_across = edge.triangle_across(self);
                if tri_info
                    .get(&triangle_across)
                    .and_then(|x| x.component)
                    .is_none()
                {
                    s.push(triangle_across)
                }
            }
        }
    }

    fn collapse_unrooted_tree(
        &self,
        tri_info: &mut HashMap<NodeTarget, TriInfo>,
        triangle: NodeTarget,
        component: &mut NonZeroUsize,
    ) {
        let mut s: Vec<NodeTarget> = Vec::new();
        s.push(triangle);
        while let Some(triangle) = s.pop() {
            // Set the component of this triangle to the current component.
            let this_tri_info = if let Some(x) = tri_info.get_mut(&triangle) {
                x
            } else {
                return;
            };
            this_tri_info.component = Some(*component);
            for edge in triangle.edges(self) {
                let edge_data = edge.edge(self);
                if edge_data.point.constraint {
                } else {
                    let next_triangle = {
                        let connected_edge = edge.sym();
                        // Find the canonical edge of this triangle.
                        let canonical_edge = connected_edge.get_tri_canonical(self);
                        canonical_edge
                    };
                    let next_tri_info = tri_info.get(&next_triangle);
                    let next_tri_component = next_tri_info.and_then(|x| x.component);
                    if next_tri_component.is_none() {
                        s.push(next_triangle);
                    }
                }
            }
        }
    }

    fn abstract_01(
        &self,
        tri_info: &mut HashMap<NodeTarget, TriInfo>,
        component: &mut NonZeroUsize,
    ) -> VecDeque<NodeTarget> {
        // q is the first queue, for tris connected to L1 nodes.
        let mut q: VecDeque<NodeTarget> = VecDeque::new();
        // r is the second queue, for all other tris we are yet to process.
        let mut r: VecDeque<NodeTarget> = VecDeque::new();
        for triangle in self.nodes() {
            // How many constrained edges does this triangle have?.
            let n_constraints = triangle.n_constrained_edges(self);
            if n_constraints == 3 {
                // We have found an L0 island. We can therefore give it its own
                // component number and we are finished with this triangle.
                tri_info.insert(
                    triangle,
                    TriInfo {
                        level: LevelInfo::L0,
                        component: Some(*component),
                    },
                );
                *component = NonZeroUsize::new(component.get() + 1).unwrap();
            // TODO: we still need to set the width and choke values.
            } else if n_constraints == 2 {
                // The triangle has two constraints and is therefore guaranteed
                // to be L1. We don't know which component it belongs to, so we
                // will leave that as 0/Null/None.
                tri_info.insert(
                    triangle,
                    TriInfo {
                        level: LevelInfo::L1(None),
                        component: None,
                    },
                );
                // We have considered this triangle and added it to the pile of
                // completed tris. Note that at some point the tri will need to
                // be revisited to add a component number, but that will be
                // handled by another triangle to which it is connected.

                // Find the unconstrained edge and add the triangle that is
                // across this edge to first queue. It is possible that a single
                // triangle is added twice to this queue.
                for edge in triangle.edges(self) {
                    let triangle_across = edge.triangle_across(self);
                    let edge_data = edge.edge(self);
                    // We only add it to q if we haven't already assigned it a level.
                    if !edge_data.point.constraint {
                        let triangle_across_info = tri_info.get(&triangle_across);
                        if triangle_across_info.is_none() {
                            q.push_back(triangle_across);
                            break;
                        } else if triangle_across_info
                            .map(|info| info.level.as_level() == Level::L1)
                            .unwrap_or(false)
                        {
                            self.collapse_unrooted_tree(tri_info, triangle_across, component);
                            *component = NonZeroUsize::new(component.get() + 1).unwrap();
                        }
                    }
                }
            } else if n_constraints == 1 || n_constraints == 0 {
                // Add to the queue of triangles that could be either 1, 2, or 3.
                r.push_back(triangle);
            } else {
                unreachable!(
                    "Each triangle should have a number of constraints in the range [0,3]"
                );
            }
        }
        // All of the triangles are on the pile or in the queues. We now cycle
        // through the first queue. First we pop a triangle off the front of the
        // queue.
        while let Some(triangle) = q.pop_front() {
            // If we don't know the level of this triangle, we must determine
            // it.

            // How any of the adjacent triangles (i.e. across an
            // unconstrained edge) are known to be level 1?
            let n_l1s = triangle.num_adjacent_level(self, tri_info, Level::L1);
            // How many of the edges of this triangle are constrained?
            let n_constraints = triangle.n_constrained_edges(self);
            // It is possible that a triangle appears multiple times in the
            // queue, so we only process it if it's tri_info or component number
            // is None. TODO: we need to consider L1 islands of only 2 elemnts,
            // currently they are not given a component number under this
            // algorithm.
            if tri_info
                .get(&triangle)
                // .and_then(|x| x.component)
                .is_none()
            {
                if (n_constraints + n_l1s) >= 2 {
                    // Set the level of the triangle to 1
                    tri_info.insert(
                        triangle,
                        TriInfo {
                            level: LevelInfo::L1(None),
                            component: None,
                        },
                    );
                    // For each of the edges, get the triangle across and if the
                    // edge is unconstrained and the level of the triangle is
                    // not set then add it to the queue.
                    for edge in triangle.edges(self) {
                        let edge_data = edge.edge(self);
                        let next_triangle = edge.sym().get_tri_canonical(self);
                        let next_triangle_info = tri_info.get(&next_triangle);
                        if !edge_data.point.constraint {}
                        if !edge_data.point.constraint && next_triangle_info.is_none() {
                            q.push_back(next_triangle);
                        }
                    }
                }
                if (n_constraints + n_l1s) == 3 {
                    self.collapse_unrooted_tree(tri_info, triangle, component);
                    *component = NonZeroUsize::new(component.get() + 1).unwrap();
                }
                if (n_constraints + n_l1s) < 2 {}
            }
        }
        r
    }

    fn abstract_3(
        &self,
        tri_info: &mut HashMap<NodeTarget, TriInfo>,
        triangle: NodeTarget,
        component: &mut NonZeroUsize,
    ) {
        let mut q: VecDeque<NodeTarget> = VecDeque::new();
        q.push_back(triangle);
        while let Some(triangle_base) = q.pop_front() {
            // triangle_base is an L3 triangle
            tri_info.insert(
                triangle_base,
                TriInfo {
                    level: LevelInfo::L3(None, None, None),
                    component: Some(*component),
                },
            );
            // edge is one of the edges of the L3 triangle, and i is the index
            // of the edge in the range [0,2], with 0 being the canonical edge.
            for (i, edge) in triangle_base.edges(self).into_iter().enumerate() {
                let edge_data = edge.edge(self);
                if edge_data.point.constraint {
                    panic!("A level 3 tri should never have a constraint as an edge");
                    // continue;
                }
                // For each edge, we walk along the L2 nodes emanating from that
                // edge. [`triangle_current`] is the L2 node we are currently
                // looking at. [`triangle_last`] is the previous L2 node or the
                // L3 node in the case of the first iteration.
                let mut triangle_current: Option<NodeTarget> = Some(edge.triangle_across(self));
                let mut triangle_last = triangle_base.clone();
                loop {
                    let mut triangle_next: Option<NodeTarget> = None;
                    let n = triangle_current.unwrap().n_constrained_edges(self);
                    let m = triangle_current
                        .unwrap()
                        .num_adjacent_level(self, tri_info, Level::L1);
                    if n + m == 0 {
                        // We have reached the next L3 node. We push this L3
                        // node onto the stack ([`q`]) for future processing.
                        if tri_info.get(&triangle_current.unwrap()).is_none() {
                            q.push_back(triangle_current.unwrap());
                        }
                        // This current L3 node [`triangle_current`] must be
                        // connected to our starting L3 node [`triangle_base`].
                        // We therefore mark [`triangle_current`] as one of the
                        // adjacent L3 nodes for [`triangle_base`].
                        {
                            let tri_info_mut = tri_info.get_mut(&triangle_base).unwrap();
                            if let LevelInfo::L3(ref mut s1, ref mut s2, ref mut s3) =
                                tri_info_mut.level
                            {
                                if i == 0 {
                                    *s1 = Some(L3NodeTarget::new_unchecked(
                                        triangle_current.unwrap().into(),
                                    ));
                                } else if i == 1 {
                                    *s2 = Some(L3NodeTarget::new_unchecked(
                                        triangle_current.unwrap().into(),
                                    ));
                                } else if i == 2 {
                                    *s3 = Some(L3NodeTarget::new_unchecked(
                                        triangle_current.unwrap().into(),
                                    ));
                                } else {
                                    panic!("invalid edge index");
                                }
                            } else {
                                panic!("incorrect level, should only be considering L3 triangles")
                            }
                        }
                        break;
                    } else if n + m == 1 {
                        // [`triangle_current`] is an L2 node.

                        // Sometimes we have already labelled the triangle as
                        // L2, but we still need to add some adjacency
                        // information. If, however, we haven't set any
                        // information, set the level to L2 now.
                        if tri_info.get(&triangle_current.unwrap()).is_none() {
                            tri_info.insert(
                                triangle_current.unwrap(),
                                TriInfo {
                                    level: LevelInfo::L2(None, None),
                                    component: Some(*component),
                                },
                            );
                        }
                        // One of the adjacency information slots needs to be
                        // filled with [`triangle_base`] as an adjacent L3 node.
                        // First we need to determine 2 things:
                        //    * Which is the edge index [0,2] of the constrained
                        //      or L1-connected edge.
                        //    * Which is the edge index [0,2] of the edge we
                        //      just stepped across to get to this node.
                        let mut constrained_edge_i: Option<usize> = None;
                        let mut entered_edge_i: Option<usize> = None;
                        let mut next_edge_i: Option<usize> = None;
                        // This L2 node has three edges (as always). One of them
                        // will be a constraint or L1 node, the other 2 will be
                        // connected to either L2 or L3 nodes.
                        for (i, edge) in triangle_current
                            .unwrap()
                            .edges(self)
                            .into_iter()
                            .enumerate()
                        {
                            let triangle_temp = edge.triangle_across(self);
                            let edge_data = edge.edge(self);
                            if triangle_temp == triangle_last {
                                // If the triangle across from this edge is
                                // [`triangle_last`] (i.e. the triangle we came
                                // from), then this is the edge we used to enter
                                // the node.
                                entered_edge_i = Some(i);
                            } else if !edge_data.point.constraint
                                && tri_info
                                    .get(&triangle_temp)
                                    .map(|x| x.level.as_level() != Level::L1)
                                    .unwrap_or(true)
                            {
                                // The triangle across from this edge is the
                                // next L2 or L3 node.
                                next_edge_i = Some(i);
                                triangle_next = Some(triangle_temp);
                            // edge_next = edge;
                            } else if !edge_data.point.constraint
                                && tri_info
                                    .get(&triangle_temp)
                                    .map(|x| x.level.as_level() == Level::L1)
                                    .unwrap_or(false)
                            {
                                // The triangle across from this edge is an L1
                                // tree.
                                constrained_edge_i = Some(i);
                                self.collapse_rooted_tree(
                                    tri_info,
                                    triangle_current.unwrap(),
                                    triangle_temp,
                                );
                            } else if edge_data.point.constraint {
                                // This edge is a constraint
                                constrained_edge_i = Some(i);
                            } else {
                                unreachable!()
                            }
                        }
                        let _constrained_edge_i: usize = constrained_edge_i.unwrap();
                        let entered_edge_i: usize = entered_edge_i.unwrap();
                        let next_edge_i: usize = next_edge_i.unwrap();
                        // Set the adjaceny information. We know the
                        // start L3 which we can attach.
                        {
                            let tri_info_mut =
                                tri_info.get_mut(&triangle_current.unwrap()).unwrap();
                            if let LevelInfo::L2(ref mut s1, ref mut s2) = tri_info_mut.level {
                                if entered_edge_i < next_edge_i {
                                    *s1 = Some(L3NodeTarget::new_unchecked(triangle_base.into()));
                                } else {
                                    *s2 = Some(L3NodeTarget::new_unchecked(triangle_base.into()));
                                }
                            } else {
                                panic!("incorrect level, should only be considering L2 triangles")
                            }
                        }
                    }

                    // TODO: this line is not in the paper, but it makes sense.
                    triangle_last = triangle_current.unwrap();
                    triangle_current = triangle_next;
                }
            }
        }
    }

    /// Return a list of triangles (denoted by canonical [`EdgeTarget`]s)
    /// between two points, including the start and end triangles. Returns
    /// [`None`] if there is no path.
    pub fn find_node_path_points(
        &self,
        tri_info: &LinkageMap,
        pa: Point,
        pb: Point,
    ) -> Option<Vec<EdgeTarget>> {
        // First we need to find the two nodes we are pathfinding between.
        let node_a = self.locate_tri(pa)?;
        let node_b = self.locate_tri(pb)?;
        self.find_node_path_nodes(tri_info, node_a, node_b)
    }

    pub fn find_node_path_nodes(
        &self,
        tri_info: &LinkageMap,
        node_a: EdgeRefA<Segment, ()>,
        node_b: EdgeRefA<Segment, ()>,
    ) -> Option<Vec<EdgeTarget>> {
        // If the two nodes are the same, we return a single path of that node.
        if node_a == node_b {
            return Some(vec![node_a.target()]);
        }
        // If the two nodes lie in different components, there is no path
        // between them.
        let tri_info_a = tri_info.get(&NodeTarget::new_unchecked(node_a.target()))?;
        let tri_info_b = tri_info.get(&NodeTarget::new_unchecked(node_b.target()))?;
        // If the two nodes are in different components there is no path between
        // the two nodes.
        if tri_info_a.component? != tri_info_b.component? {
            return None;
        }
        todo!()
    }
}

#[derive(Clone, Debug)]
pub struct IntersectionIter<'a> {
    triangulation: &'a ConstrainedTriangulation,
    current_intersection: Option<Intersection>,
    a: Point,
    b: Point,
}

impl<'a> IntersectionIter<'a> {
    pub fn new(triangulation: &'a ConstrainedTriangulation, a: Point, b: Point) -> Self {
        Self {
            triangulation,
            current_intersection: None,
            a,
            b,
        }
    }
}

impl<'a> Iterator for IntersectionIter<'a> {
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
                    if intersecting_edge.edge(self.triangulation).point.point() == self.b {
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

#[derive(Clone, Debug)]
pub struct LinkageMap(pub HashMap<NodeTarget, TriInfo>);

impl LinkageMap {
    pub fn new() -> Self {
        LinkageMap(HashMap::new())
    }
    pub fn from(map: HashMap<NodeTarget, TriInfo>) -> Self {
        LinkageMap(map)
    }

    pub fn get(&self, k: &NodeTarget) -> Option<&TriInfo> {
        self.0.get(k)
    }

    pub fn ns(&self) -> (usize, usize, usize, usize) {
        let mut n_l0 = 0;
        let mut n_l1 = 0;
        let mut n_l2 = 0;
        let mut n_l3 = 0;
        for (_, v) in self.0.iter() {
            match v.level {
                LevelInfo::L0 => n_l0 += 1,
                LevelInfo::L1(..) => n_l1 += 1,
                LevelInfo::L2(..) => n_l2 += 1,
                LevelInfo::L3(..) => n_l3 += 1,
            }
        }
        (n_l0, n_l1, n_l2, n_l3)
    }
}

impl<'a> EdgeRefA<'a, Segment, ()> {
    pub fn triangle_across(&self) -> Self {
        let edge_on_next_tri = self.sym();
        let canonical_edge = edge_on_next_tri.get_tri_canonical();
        canonical_edge
    }
    fn n_constrained_edges(&self) -> usize {
        let initial_edge = *self;
        let mut current_edge = initial_edge;
        let mut n_constraints: usize = 0;
        loop {
            if current_edge.edge().point.constraint {
                n_constraints += 1;
            }
            current_edge = current_edge.l_next();
            if current_edge == initial_edge {
                break;
            }
        }
        n_constraints
    }

    // How many of the triangles adjacent to this triangle (across an
    // unconstrained edge) are of the give level?
    fn num_adjacent_level(&self, tri_info: &HashMap<NodeTarget, TriInfo>, level: Level) -> usize {
        let initial_edge = *self;
        let mut current_edge = initial_edge;
        let mut n_tris: usize = 0;
        loop {
            // First check if the edge is a constraint. If it's a constraint we
            // just skip it.
            if !current_edge.edge().point.constraint {
                // The edge is not a constraint so we will test it. First we
                // need to get an edge in the adjacent triangle.
                let adjacent_edge = current_edge.sym();
                // Then we need to get the canonical edge.
                let canonical_edge = adjacent_edge.get_tri_canonical();
                // Then we need to see if we know the level of the triangle.
                let this_level: Option<Level> = tri_info
                    .get(&NodeTarget::new_unchecked(canonical_edge.target()))
                    .map(|x| x.level.into());
                if this_level == Some(level) {
                    n_tris += 1;
                }
            }
            current_edge = current_edge.l_next();
            if current_edge == initial_edge {
                break;
            }
        }
        n_tris
    }
}

#[derive(Copy, Clone, Debug, Ord, PartialOrd, Eq, PartialEq, Hash)]
pub struct TriInfo {
    pub level: LevelInfo,
    pub component: Option<NonZeroUsize>,
}

#[derive(Copy, Clone, Debug, Ord, PartialOrd, Eq, PartialEq, Hash)]
pub enum Level {
    L0,
    L1,
    L2,
    L3,
}

impl From<Level> for usize {
    fn from(level: Level) -> Self {
        match level {
            Level::L0 => 0,
            Level::L1 => 1,
            Level::L2 => 2,
            Level::L3 => 3,
        }
    }
}

impl From<LevelInfo> for Level {
    fn from(level: LevelInfo) -> Self {
        match level {
            LevelInfo::L0 => Level::L0,
            LevelInfo::L1(..) => Level::L1,
            LevelInfo::L2(..) => Level::L2,
            LevelInfo::L3(..) => Level::L3,
        }
    }
}

impl From<&LevelInfo> for Level {
    fn from(level: &LevelInfo) -> Self {
        match level {
            LevelInfo::L0 => Level::L0,
            LevelInfo::L1(..) => Level::L1,
            LevelInfo::L2(..) => Level::L2,
            LevelInfo::L3(..) => Level::L3,
        }
    }
}

/// Each level also has adjacency information which varies depending on the
/// level.
///    * *L0 Nodes:* Has no adjacency information.
///    * *L1 Nodes:* Contains the L2 node at the root of the tree if it has one,
///      or nothing if it's unrooted.
///    * *L2 Nodes:* Contains 2 nodes which are the L3 nodes at either end of
///      the corridor. It is possible that it is the same node. Contains nothing
///      if it is an L2 ring (i.e. the component does not contain any L3 nodes,
///      although it may contain L1 nodes).
///    * *L3 Nodes:* Contain the 3 other L3 nodes it is connected to at the
///      other end of the corridors.
#[derive(Copy, Clone, Debug, Ord, PartialOrd, Eq, PartialEq, Hash)]
pub enum LevelInfo {
    L0,
    L1(Option<NodeTarget>),
    L2(Option<L3NodeTarget>, Option<L3NodeTarget>),
    L3(
        Option<L3NodeTarget>,
        Option<L3NodeTarget>,
        Option<L3NodeTarget>,
    ),
}

impl LevelInfo {
    pub fn as_level(&self) -> Level {
        self.into()
    }
    pub fn is_ring_or_loop(&self) -> bool {
        if let LevelInfo::L2(a, b) = self {
            a == b
        } else {
            false
        }
    }
}

impl From<LevelInfo> for usize {
    fn from(level: LevelInfo) -> Self {
        match level {
            LevelInfo::L0 => 0,
            LevelInfo::L1(..) => 1,
            LevelInfo::L2(..) => 2,
            LevelInfo::L3(..) => 3,
        }
    }
}

#[derive(Clone, Copy)]
pub struct TriangleIter<'a> {
    triangulation: &'a ConstrainedTriangulation,
    next: EdgeTarget,
}

impl<'a> TriangleIter<'a> {
    pub fn new(triangulation: &'a ConstrainedTriangulation) -> Self {
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

impl<'a> Iterator for TriangleIter<'a> {
    type Item = EdgeRefA<'a, Segment, ()>;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // First we check that the edge actually exists for this given
            // EdgeTarget.
            if self.triangulation.qeds.quads.contains(self.next.e) {
                let edge_ref = unsafe { self.triangulation.qeds.edge_a_ref(self.next) };
                self.inc();
                let is_tri_canonical: bool = edge_ref.is_tri_canonical();
                if is_tri_canonical {
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

// TODO: does this actuall return all the nodes?
#[derive(Clone, Copy)]
pub struct NodeIter<'a> {
    triangulation: &'a ConstrainedTriangulation,
    next: EdgeTarget,
}

impl<'a> NodeIter<'a> {
    pub fn new(triangulation: &'a ConstrainedTriangulation) -> Self {
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

impl<'a> Iterator for NodeIter<'a> {
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

impl<'a, BData> EdgeRefA<'a, Segment, BData> {
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
        for edge in self.l_face().edges {
            match edge.target().e {
                0 | 1 | 2 | 4 => return false,
                _ => (),
            }
        }
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

impl EdgeTarget {
    /// An edge is the canonical edge for a tri iff it has the lowest e value
    /// for the tri.
    fn is_tri_canonical(self, triangulation: &ConstrainedTriangulation) -> bool {
        let edge_ref = unsafe { triangulation.qeds.edge_a_ref(self) };
        edge_ref.is_tri_canonical()
    }

    fn is_tri_real(self, triangulation: &ConstrainedTriangulation) -> bool {
        let edge_ref = unsafe { triangulation.qeds.edge_a_ref(self) };
        edge_ref.is_tri_real()
    }

    pub fn get_tri_canonical(self, triangulation: &ConstrainedTriangulation) -> NodeTarget {
        let edge_ref = unsafe { triangulation.qeds.edge_a_ref(self) };
        let canonical_edge = edge_ref.get_tri_canonical();
        NodeTarget(canonical_edge.target())
    }

    pub fn triangle_across(self, triangulation: &ConstrainedTriangulation) -> NodeTarget {
        let edge_ref = unsafe { triangulation.qeds.edge_a_ref(self) };
        let triangle_across = edge_ref.triangle_across();
        NodeTarget(triangle_across.target())
    }

    pub fn edge(self, triangulation: &ConstrainedTriangulation) -> Edge<Segment> {
        let edge_ref = unsafe { triangulation.qeds.edge_a_ref(self) };
        edge_ref.edge().clone()
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
    use quickcheck::quickcheck;

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
    fn empty_triangulation() {
        let triangulation = ConstrainedTriangulation::new();
        let tris: Vec<EdgeRefA<Segment, ()>> = triangulation.triangles().collect();
        assert_eq!(tris.len(), 2);
        let tri_info = triangulation.classify_triangles();
        for (target, info) in tri_info.0 {
            println!("EdgeTarget: {:?} Info: {:?}", target, info);
        }
    }

    #[test]
    fn one_point_triangulation_location() {
        let mut triangulation = ConstrainedTriangulation::new();
        // assert_eq!(triangulation.qeds().unwrap().quads.len(), 5);
        let south = unsafe { triangulation.qeds.edge_a_ref(EdgeTarget::new(0, 2, 0)) };
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
            EdgeTarget::new(2, 0, 0)
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

        let centre = unsafe { triangulation.qeds.edge_a_ref(EdgeTarget::new(3, 0, 0)) };
        assert_eq!(
            centre.edge().point.point,
            Point::new(SafeFloat::MAX.0, SafeFloat::MAX.0)
        );
        assert_eq!(
            centre.sym().edge().point.point,
            Point::new(SafeFloat::MIN.0, SafeFloat::MIN.0)
        );

        let north = unsafe { triangulation.qeds.edge_a_ref(EdgeTarget::new(2, 0, 0)) };
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
        triangulation.add_constraint(p1, p2);
        triangulation.add_constraint(p2, p3);
        triangulation.add_constraint(p3, p1);
        assert_eq!(triangulation.triangles().count(), 8);
        let tri_info = triangulation.classify_triangles();
        assert_eq!(tri_info.ns(), (1, 0, 7, 0));
        for (target, info) in tri_info.0 {
            println!("EdgeTarget: {:?} Info: {:?}", target, info);
        }
    }

    #[test]
    fn quad_triangulation() {
        let mut triangulation = ConstrainedTriangulation::new();
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(5.0, 0.10);
        let p3 = Point::new(2.5, 5.1);
        let p4 = Point::new(5.0, 5.0);
        triangulation.add_constraint(p1, p2);
        triangulation.add_constraint(p3, p1);
        triangulation.add_constraint(p4, p2);
        triangulation.add_constraint(p4, p3);
        // 10 triangles, 2 from the quad and 8 from the bounding box.
        assert_eq!(triangulation.triangles().count(), 10);
        let tri_info = triangulation.classify_triangles();
        assert_eq!(tri_info.ns(), (0, 2, 8, 0));
        for (target, info) in tri_info.0 {
            println!("EdgeTarget: {:?} Info: {:?}", target, info);
        }
    }

    #[test]
    fn star_triangulation() {
        let mut triangulation = ConstrainedTriangulation::new();
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(0.5, 0.5);
        let p3 = Point::new(-0.5, 0.5);
        let p4 = Point::new(-2.0, -1.0);
        let p5 = Point::new(2.0, -1.0);
        let p6 = Point::new(0.0, 2.0);
        triangulation.add_constraint(p4, p1);
        triangulation.add_constraint(p1, p5);
        triangulation.add_constraint(p5, p2);
        triangulation.add_constraint(p2, p6);
        triangulation.add_constraint(p6, p3);
        triangulation.add_constraint(p3, p4);

        // 14 triangles, 4 from the quad and 10 from the bounding box.
        assert_eq!(triangulation.triangles().count(), 14);
        let tri_info = triangulation.classify_triangles();
        assert_eq!(tri_info.ns(), (0, 7, 7, 0));
    }

    #[test]
    fn ring_triangulation() {
        let mut triangulation = ConstrainedTriangulation::new();
        let p1 = Point::new(-1.0, -1.0);
        let p2 = Point::new(1.0, -1.0);
        let p3 = Point::new(0.0, 1.0);
        let p4 = Point::new(-2.0, -1.5);
        let p5 = Point::new(2.0, -1.5);
        let p6 = Point::new(0.0, 2.0);
        triangulation.add_constraint(p1, p2);
        triangulation.add_constraint(p2, p3);
        triangulation.add_constraint(p3, p1);
        triangulation.add_constraint(p4, p5);
        triangulation.add_constraint(p5, p6);
        triangulation.add_constraint(p6, p4);

        // 14 triangles, 7 from the quad and 7 from the bounding box.
        assert_eq!(triangulation.triangles().count(), 14);
        let tri_info = triangulation.classify_triangles();
        assert_eq!(tri_info.ns(), (1, 0, 13, 0));
    }

    #[test]
    fn figure8_triangulation() {
        let mut triangulation = ConstrainedTriangulation::new();
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(5.0, 0.0);
        let p3 = Point::new(5.0, 3.0);
        let p4 = Point::new(0.0, 3.0);

        let q1 = Point::new(1.0, 1.0);
        let q2 = Point::new(2.0, 1.0);
        let q3 = Point::new(2.0, 2.0);
        let q4 = Point::new(1.0, 2.0);

        let r = Point::new(2.0, 0.0);

        let r1 = r + q1;
        let r2 = r + q2;
        let r3 = r + q3;
        let r4 = r + q4;

        triangulation.add_constraint(p1, p2);
        triangulation.add_constraint(p2, p3);
        triangulation.add_constraint(p3, p4);
        triangulation.add_constraint(p4, p1);

        triangulation.add_constraint(q1, q2);
        triangulation.add_constraint(q2, q3);
        triangulation.add_constraint(q3, q4);
        triangulation.add_constraint(q4, q1);

        triangulation.add_constraint(r1, r2);
        triangulation.add_constraint(r2, r3);
        triangulation.add_constraint(r3, r4);
        triangulation.add_constraint(r4, r1);

        // 26 triangles, 18 from the quad and 8 from the bounding box.
        assert_eq!(triangulation.triangles().count(), 26);
        let _tri_info = triangulation.classify_triangles();
        // assert_eq!(tri_info.ns(),(1,0,13,0));
    }

    #[test]
    fn l_shaped_triangulation() {
        let mut triangulation = ConstrainedTriangulation::new();
        // triangulation.add_point(Point::new(-5.5, 5.5));
        // triangulation.add_point(Point::new(-4.5, 5.5));

        triangulation.add_constraint(Point::new(-5.5, 5.5), Point::new(-4.5, 5.5));
        triangulation.add_constraint(Point::new(-4.5, 5.5), Point::new(-4.5, 4.5));

        // 8 triangles, 1 from the l-shape and 7 from the bounding box.
        assert_eq!(triangulation.triangles().count(), 8);
        // let _tri_info = triangulation.classify_triangles();
        // assert_eq!(tri_info.ns(),(1,0,13,0));
    }

    #[test]
    fn u_shaped_triangulation() {
        let mut triangulation = ConstrainedTriangulation::new();
        triangulation.add_constraint(Point::new(-5.5, 5.5), Point::new(-4.5, 5.5));
        triangulation.add_constraint(Point::new(-4.5, 5.5), Point::new(-4.5, 4.5));
        triangulation.add_constraint(Point::new(-5.5, 5.5), Point::new(-5.5, 4.5));

        // 26 triangles, 18 from the quad and 8 from the bounding box.
        // assert_eq!(triangulation.triangles().count(), 26);
        let _tri_info = triangulation.classify_triangles();
        // assert_eq!(tri_info.ns(),(1,0,13,0));
    }

    #[test]
    fn u2_shaped_triangulation() {
        let mut triangulation = ConstrainedTriangulation::new();
        triangulation.add_constraint(Point::new(-5.5, 2.5), Point::new(-4.5, 2.5));
        triangulation.add_constraint(Point::new(-5.5, 1.5), Point::new(-4.5, 1.5));
        triangulation.add_constraint(Point::new(-5.5, 2.5), Point::new(-5.5, 1.5));
        triangulation.add_constraint(Point::new(-3.5, 2.5), Point::new(-3.5, 1.5));

        // 26 triangles, 18 from the quad and 8 from the bounding box.
        // assert_eq!(triangulation.triangles().count(), 26);
        let _tri_info = triangulation.classify_triangles();
        // assert_eq!(tri_info.ns(),(1,0,13,0));
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
        for (i, _quad) in triangulation.qeds().unwrap().quads.iter() {
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
            triangulation
                .find_intersections_between_points(Point::new(2.0, 2.0), Point::new(4.0, 3.0))
                .map(|x| {
                    if let Intersection::Point(p) = x {
                        println!(
                            "IntersectionPoint: {}",
                            p.edge(&triangulation).point.point()
                        );
                    }
                    x
                })
                .collect::<Vec<_>>(),
            vec![]
        );
    }

    #[test]
    fn single_intersection() {
        let mut triangulation = ConstrainedTriangulation::new();
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(5.0, -1.0);
        let p3 = Point::new(10.0, 0.0);
        let p4 = Point::new(5.0, 1.0);
        triangulation.add_constraint(p1, p2);
        triangulation.add_constraint(p2, p3);
        triangulation.add_constraint(p3, p4);
        triangulation.add_constraint(p4, p1);
        triangulation.add_constraint(p1, p3);
    }

    #[test]
    fn double_intersection() {
        let mut triangulation = ConstrainedTriangulation::new();
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(5.0, -1.0);
        let p3 = Point::new(10.0, 0.0);
        let p4 = Point::new(5.5, 1.0);
        let p5 = Point::new(4.5, 1.0);
        triangulation.add_constraint(p1, p2);
        triangulation.add_constraint(p2, p3);
        triangulation.add_constraint(p3, p4);
        triangulation.add_constraint(p4, p5);
        triangulation.add_constraint(p5, p1);
        triangulation.add_constraint(p1, p3);
    }

    #[test]
    fn regress1() {
        let mut triangulation = ConstrainedTriangulation::new();
        let segments = vec![
            (
                Point::new(-75.3645612469715, 10.649471266802962),
                Point::new(11.264242662961536, -74.1427143080474),
            ),
            (
                Point::new(-65.08752072884896, 85.23877185897558),
                Point::new(-36.55529677285707, 25.159802183655742),
            ),
            (
                Point::new(-71.01084733411076, 30.660749902036656),
                Point::new(67.62855075915658, -89.10279376500583),
            ),
        ];
        for (p1, p2) in segments.iter() {
            triangulation.add_constraint(*p1, *p2);
        }
        valid_triangulation(&triangulation);
    }

    #[test]
    fn regress2() {
        let mut triangulation = ConstrainedTriangulation::new();
        let segments = vec![
            (
                Point {
                    x: -40.292118735040106,
                    y: 82.0073097016039,
                },
                Point {
                    x: -42.07946675183436,
                    y: -44.51802168917966,
                },
            ),
            (
                Point {
                    x: -52.00978628295401,
                    y: 70.99340776302836,
                },
                Point {
                    x: 60.41164157780554,
                    y: 77.45771825005286,
                },
            ),
            (
                Point {
                    x: 87.60492632142518,
                    y: 13.83554980163639,
                },
                Point {
                    x: -64.2901267321343,
                    y: 68.80564964035884,
                },
            ),
            (
                Point {
                    x: 40.89557780924511,
                    y: 34.30455081240521,
                },
                Point {
                    x: -96.09647543301776,
                    y: 93.8999414141874,
                },
            ),
            (
                Point {
                    x: -45.63618305392634,
                    y: 73.74453797959046,
                },
                Point {
                    x: -44.04842079916467,
                    y: 19.50193303235106,
                },
            ),
        ];
        for (p1, p2) in segments.iter() {
            triangulation.add_constraint(*p1, *p2);
        }
        valid_triangulation(&triangulation);
    }

    #[test]
    fn regress3() {
        let mut triangulation = ConstrainedTriangulation::new();
        let segments = vec![
            (
                Point {
                    x: -0.4214359290103573,
                    y: 95.82816145300376,
                },
                Point {
                    x: -9.699326930357884,
                    y: 26.125595226387333,
                },
            ),
            (
                Point {
                    x: -13.147969911139384,
                    y: -54.45197788876022,
                },
                Point {
                    x: -44.33529262658724,
                    y: -75.64912335744887,
                },
            ),
            (
                Point {
                    x: -16.705404907959036,
                    y: 35.1977871928658,
                },
                Point {
                    x: 60.46183760891148,
                    y: 10.128488435745282,
                },
            ),
            (
                Point {
                    x: -59.73082630206363,
                    y: -57.035655725149056,
                },
                Point {
                    x: 49.69180679531209,
                    y: 96.08879129160505,
                },
            ),
            (
                Point {
                    x: 12.685519199969605,
                    y: 78.74650844233386,
                },
                Point {
                    x: -0.5177468902998896,
                    y: -41.830391168991895,
                },
            ),
            (
                Point {
                    x: -91.33303213105921,
                    y: 42.425551365690154,
                },
                Point {
                    x: 67.18820631566183,
                    y: 49.651436020714186,
                },
            ),
            (
                Point {
                    x: 19.547046919264545,
                    y: 21.14100619435102,
                },
                Point {
                    x: 43.67012362837994,
                    y: -56.81803213245602,
                },
            ),
        ];
        for (p1, p2) in segments.iter() {
            triangulation.add_constraint(*p1, *p2);
        }
        valid_triangulation(&triangulation);
    }

    #[test]
    fn regress4() {
        let mut triangulation = ConstrainedTriangulation::new();
        let segments = vec![
            (
                Point {
                    x: 49.354794208915905,
                    y: 77.9312265303424,
                },
                Point {
                    x: -11.492740263412088,
                    y: 8.956223279493656,
                },
            ),
            (
                Point {
                    x: 90.35804944791943,
                    y: -64.48450750858385,
                },
                Point {
                    x: 25.29536887506309,
                    y: 8.406416670169662,
                },
            ),
            (
                Point {
                    x: 46.09321063375256,
                    y: 61.22935053925707,
                },
                Point {
                    x: -26.439979432037305,
                    y: 51.94522246412245,
                },
            ),
            (
                Point {
                    x: 96.6694613097255,
                    y: 45.14085139658687,
                },
                Point {
                    x: 58.35546494466675,
                    y: 53.009402448096495,
                },
            ),
            (
                Point {
                    x: 58.99661219941143,
                    y: 65.20711681694809,
                },
                Point {
                    x: -63.66585352841398,
                    y: -59.20567981731186,
                },
            ),
        ];
        for (p1, p2) in segments.iter() {
            triangulation.add_constraint(*p1, *p2);
        }
        valid_triangulation(&triangulation);
    }

    #[test]
    fn regress5() {
        let mut triangulation = ConstrainedTriangulation::new();
        let segments = vec![
            (
                Point {
                    x: -51.400285762967044,
                    y: 78.08416289394077,
                },
                Point {
                    x: 73.94719688650767,
                    y: -9.952271890507092,
                },
            ),
            (
                Point {
                    x: -80.05915835704234,
                    y: 97.55389372322861,
                },
                Point {
                    x: 75.64795651512509,
                    y: -27.169794412596275,
                },
            ),
            (
                Point {
                    x: 81.99617838589236,
                    y: 8.336094155253178,
                },
                Point {
                    x: 99.33459619053124,
                    y: -28.76119989156618,
                },
            ),
            (
                Point {
                    x: 52.06631213621819,
                    y: 91.79642677407745,
                },
                Point {
                    x: 98.04001031185186,
                    y: -81.61857001089774,
                },
            ),
        ];
        for (p1, p2) in segments.iter() {
            triangulation.add_constraint(*p1, *p2);
        }
        valid_triangulation(&triangulation);
    }

    #[test]
    fn regress6() {
        let mut triangulation = ConstrainedTriangulation::new();
        let segments = vec![
            (
                Point {
                    x: -71.48547038818833,
                    y: -74.39322987310905,
                },
                Point {
                    x: -89.94096477395534,
                    y: 19.301100606945766,
                },
            ),
            (
                Point {
                    x: -86.9828479383787,
                    y: 0.05961466115694236,
                },
                Point {
                    x: 43.25946276642907,
                    y: -77.91468968915898,
                },
            ),
            (
                Point {
                    x: -4.142054219662228,
                    y: 70.85270987086068,
                },
                Point {
                    x: -90.96995169850732,
                    y: 36.77453032058807,
                },
            ),
            (
                Point {
                    x: -73.374254157982,
                    y: 87.53245031954737,
                },
                Point {
                    x: -46.959018139029226,
                    y: -48.21755053303174,
                },
            ),
            (
                Point {
                    x: 97.61862767011499,
                    y: 56.159609123931716,
                },
                Point {
                    x: -60.03876659083822,
                    y: -4.969101911120703,
                },
            ),
            (
                Point {
                    x: -99.10800691452306,
                    y: 43.80328990194849,
                },
                Point {
                    x: -58.187890315932435,
                    y: -11.073968979166835,
                },
            ),
        ];
        for (p1, p2) in segments.iter() {
            triangulation.add_constraint(*p1, *p2);
        }
        valid_triangulation(&triangulation);
    }

    #[ignore]
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
            let intersections: Vec<_> = triangulation
                .find_intersections_between_points(Point::new(2.0, 2.0), Point::new(9.0, 8.0))
                .collect();
            println!("intersections: {:?}", intersections);
            assert_eq!(intersections.len(), 1);
            let first_intersection =
                unsafe { triangulation.qeds.edge_a_ref(intersections[0].target()) };
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
            let intersections = triangulation
                .find_intersections_between_points(Point::new(2.0, 2.0), Point::new(18.0, 6.0));
            assert_eq!(intersections.clone().count(), 2);
            let first_intersection = unsafe {
                triangulation
                    .qeds
                    .edge_a_ref(intersections.clone().collect::<Vec<_>>()[0].target())
            };
            // The edges found have a specified direction.
            assert_eq!(
                (
                    first_intersection.edge().point.point(),
                    first_intersection.sym().edge().point.point()
                ),
                (p2, p3)
            );
            let second_intersection = unsafe {
                triangulation
                    .qeds
                    .edge_a_ref(intersections.collect::<Vec<_>>()[1].target())
            };
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
            let intersections = triangulation
                .find_intersections_between_points(Point::new(6.0, 5.0), Point::new(9.0, 5.0));
            assert_eq!(intersections.count(), 0);
            // let first_intersection = unsafe {triangulation.qeds.edge_a_ref(intersections[0])};
            // // The edges found have a specified direction.
            // assert_eq!((first_intersection.edge().point.point(),first_intersection.sym().edge().point.point()),(p2,p3));
            // let second_intersection = unsafe {triangulation.qeds.edge_a_ref(intersections[1])};
            // assert_eq!((second_intersection.edge().point.point(),second_intersection.sym().edge().point.point()),(p2,p4));
        }
        {
            println!("Testing from a vertex (without leaving initial triangle)");
            let intersections = triangulation
                .find_intersections_between_points(Point::new(0.0, 0.0), Point::new(1.0, 1.0));
            assert_eq!(intersections.count(), 0);
        }
        {
            // through a vertex (but starting within a triangle).
            println!("Testing from a vertex (with 1 intersection)");
            let intersections = triangulation
                .find_intersections_between_points(Point::new(0.0, 0.0), Point::new(8.0, 10.0));
            assert_eq!(intersections.clone().count(), 1);

            let first_intersection = unsafe {
                triangulation
                    .qeds
                    .edge_a_ref(intersections.collect::<Vec<_>>()[0].target())
            };
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
        // let p3 = Point::new(5.0, 10.0);
        triangulation.add_constraint(p1, p2);
        assert!(has_constraint(&triangulation, p1, p2));
        // triangulation.add_point(p2);
        // triangulation.add_point(p3);
        // valid_triangulation(&triangulation);
        // assert_eq!(
        // triangulation.find_intersections_between_points(Point::new(2.0, 2.0), Point::new(4.0, 3.0)),
        // vec![]
        // );
    }

    #[test]
    fn cross_constraints() {
        let mut triangulation = ConstrainedTriangulation::new();
        let p1 = Point::new(-10.0, 0.0);
        let p2 = Point::new(10.0, 0.0);
        let p3 = Point::new(0.0, -10.0);
        let p4 = Point::new(0.0, 10.0);
        triangulation.add_constraint(p1, p2);
        triangulation.add_constraint(p3, p4);
        let p0 = Point::new(0.0, 0.0);
        assert!(has_constraint(&triangulation, p0, p1));
        assert!(has_constraint(&triangulation, p0, p2));
        assert!(has_constraint(&triangulation, p0, p3));
        assert!(has_constraint(&triangulation, p0, p4));
        let p5 = Point::new(3.0, 3.0);
        let p6 = Point::new(6.0, 6.0);
        triangulation.add_constraint(p5, p6);
        // assert!(has_constraint(&triangulation, p5, p6));
        // let n_constraints = triangulation
        //     .qeds
        //     .quads
        //     .iter()
        //     .filter(|(_i, q)| q.edges_a[0].point.constraint)
        //     .count();
        // assert_eq!(n_constraints, 9);
    }

    #[test]
    fn path_iter() {
        let mut full_map = FullMap::new();
        {
            let triangulation = &mut full_map.triangulation;
            let p1 = Point::new(0.0, 0.0);
            let p2 = Point::new(5.0, 0.0);
            let p3 = Point::new(5.0, 3.0);
            let p4 = Point::new(0.0, 3.0);

            let q1 = Point::new(1.0, 1.0);
            let q2 = Point::new(2.0, 1.0);
            let q3 = Point::new(2.0, 2.0);
            let q4 = Point::new(1.0, 2.0);

            let r = Point::new(2.0, 0.0);

            let r1 = r + q1;
            let r2 = r + q2;
            let r3 = r + q3;
            let r4 = r + q4;

            triangulation.add_constraint(p1, p2);
            triangulation.add_constraint(p2, p3);
            triangulation.add_constraint(p3, p4);
            triangulation.add_constraint(p4, p1);

            triangulation.add_constraint(q1, q2);
            triangulation.add_constraint(q2, q3);
            triangulation.add_constraint(q3, q4);
            triangulation.add_constraint(q4, q1);

            triangulation.add_constraint(r1, r2);
            triangulation.add_constraint(r2, r3);
            triangulation.add_constraint(r3, r4);
            triangulation.add_constraint(r4, r1);
        }
        full_map.triangulation.classify_triangles();
        let e = L3NodeTarget::new_unchecked(EdgeTarget::new(0, 0, 0));
        let path = L3Path::new(e);
        let nodes: Vec<NodeTarget> = path.node_iter(&full_map).collect();
        assert_eq!(nodes, vec![e.into()]);
    }

    #[test]
    fn unrooted_example() {
        let mut full_map = FullMap::new();
        {
            let triangulation = &mut full_map.triangulation;
            let a = Point::new(0.0, 1.0);
            let b = Point::new(0.0, 2.0);
            let c = Point::new(2.0, 2.0);
            let d = Point::new(2.0, 3.0);
            let e = Point::new(2.0, 4.0);
            let f = Point::new(5.0, 4.0);
            let g = Point::new(5.0, 3.0);
            let h = Point::new(4.0, 3.0);
            let i = Point::new(4.0, 1.0);
            let j = Point::new(6.0, 1.0);
            let k = Point::new(6.0, 0.0);
            let l = Point::new(3.0, 0.0);
            let m = Point::new(2.0, 0.0);
            let n = Point::new(2.0, 1.0);
            let o = Point::new(3.0, 1.0);
            let p = Point::new(3.0, 2.0);
            let q = Point::new(3.0, 3.0);
            let r = Point::new(5.0, 1.0);

            triangulation.add_constraint(a, b);
            triangulation.add_constraint(b, p);
            triangulation.add_constraint(c, e);
            triangulation.add_constraint(d, q);
            triangulation.add_constraint(e, f);
            triangulation.add_constraint(f, g);
            triangulation.add_constraint(g, h);
            triangulation.add_constraint(h, i);
            triangulation.add_constraint(i, j);
            triangulation.add_constraint(j, k);
            triangulation.add_constraint(k, m);
            triangulation.add_constraint(l, o);
            triangulation.add_constraint(m, n);
            triangulation.add_constraint(n, a);
            triangulation.add_constraint(g, r);
        }
        full_map.classification = full_map.triangulation.classify_triangles();
        debug_assert_components_present(&full_map);
    }

    #[test]
    fn coincident_intersection() {
        let mut full_map = FullMap::new();
        {
            let triangulation = &mut full_map.triangulation;
            let a = Point::new(0.0, 0.0);
            let b = Point::new(20.0, 0.0);
            let c = Point::new(20.0, 10.0);
            let d = Point::new(0.0, 10.0);
            let e = Point::new(10.0, 0.0);
            let f = Point::new(10.0, 6.0);
            let g = Point::new(9.0, 6.0);
            let h = Point::new(11.0, 6.0);
            triangulation.add_point(a);
            triangulation.add_point(b);
            triangulation.add_constraint(a, b);
            triangulation.add_constraint(b, c);
            triangulation.add_constraint(c, d);
            triangulation.add_constraint(d, a);
            triangulation.add_constraint(e, f);
            triangulation.add_point(g);
            triangulation.add_point(h);
            let intersections = triangulation.find_intersections_between_points(g, h);
            assert_eq!(intersections.count(), 2);
            triangulation.add_constraint(g, h);
        }
        full_map.classification = full_map.triangulation.classify_triangles();
        debug_assert_components_present(&full_map);
    }

    #[quickcheck]
    fn single_random_segment(segment: (Point, Point)) -> bool {
        let mut full_map = FullMap::new();
        {
            let triangulation = &mut full_map.triangulation;
            triangulation.add_constraint(segment.0, segment.1);
        }
        full_map.classification = full_map.triangulation.classify_triangles();
        debug_assert_components_present(&full_map);
        true
    }

    #[ignore]
    #[quickcheck]
    fn many_random_segment(segments: Vec<(Point, Point)>) -> bool {
        let mut full_map = FullMap::new();
        {
            let triangulation = &mut full_map.triangulation;
            for segment in segments.iter() {
                println!("segment: {}-{}", segment.0, segment.1);
                triangulation.add_constraint(segment.0, segment.1);
            }
        }
        full_map.classification = full_map.triangulation.classify_triangles();
        debug_assert_components_present(&full_map);
        true
    }
}

#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub enum IntersectionResult {
    Parallel,
    LineIntersection(f64, f64),
}

fn segment_intersection(s1: (Point, Point), s2: (Point, Point)) -> IntersectionResult {
    let p = s1.0;
    let r_abs = s1.1;
    let q = s2.0;
    let s_abs = s2.1;
    // line1 t = p + t*r
    // line2 r = q + u*s
    let r = r_abs - p;
    let s = s_abs - q;

    // Sometimes we get numerical errors when the points coincide, so we
    // check for that first
    let t = if p == q {
        0.0
    } else if p == s_abs {
        0.0
    } else if r_abs == q {
        1.0
    } else if r_abs == s_abs {
        1.0
    } else {
        cross(q - p, s) / cross(r, s)
    };

    let u = if q == p {
        0.0
    } else if q == r_abs {
        0.0
    } else if s_abs == p {
        1.0
    } else if s_abs == r_abs {
        1.0
    } else {
        cross(q - p, r) / cross(r, s)
    };

    // If r `cross` s is 0 then the lines are parallel and do not intersect, so
    // return nothing.
    let cross_val = cross(r, s);
    if cross_val == 0.0 {
        IntersectionResult::Parallel
    } else {
        IntersectionResult::LineIntersection(t, u)
    }
}

fn cross(pa: Point, pb: Point) -> f64 {
    pa.x * pb.y - pb.x * pa.y
}

fn show_triangle(triangle: EdgeRefA<Segment, ()>) -> String {
    format!(
        "{}-{}-{}",
        triangle.edge().point.point(),
        triangle.l_next().edge().point.point(),
        triangle.l_next().l_next().edge().point.point()
    )
}

fn show_edge(edge: EdgeRefA<Segment, ()>) -> String {
    format!(
        "{}-{}",
        edge.edge().point.point(),
        edge.sym().edge().point.point(),
    )
}

/// Assert that every node as a component.
pub fn debug_assert_components_present(full_map: &FullMap) {
    // // Check that every node that has a level has a component
    // for (node, tri_info) in full_map.classification.0.iter() {
    //     assert!(tri_info.component.is_some(), "{:?} ({:?}) does not have a component number", node, tri_info);
    // }

    // Check that every node has a complete tri_info
    for node in full_map.triangulation.nodes() {
        let tri_info = full_map.classification.get(&node).unwrap();
        debug_assert!(
            tri_info.component.is_some(),
            "{:?} ({:?}) does not have a component number",
            node,
            tri_info
        );
    }

    // Check that every node has a complete tri_info
    for node in full_map.triangulation.nodes() {
        let tri_info = full_map.classification.get(&node).unwrap();
        debug_assert!(
            tri_info.component.is_some(),
            "{:?} ({:?}) does not have a component number",
            node,
            tri_info
        );
        let component_number = tri_info.component.unwrap();
        // Check that every node that is connected has the correct component number.
        for node in node.adjacent_tris(&full_map.triangulation).into_iter() {
            let tri_info = full_map.classification.get(&node).unwrap();
            debug_assert_eq!(tri_info.component, Some(component_number));
        }
    }
}
