//! A surface triangulation is a triangulation with boundaries.
use crate::point::*;
use crate::qeds::*;
use crate::triangulation::cocircular;
use crate::triangulation::del_test_ccw;
use crate::triangulation::is_ccw;
use crate::triangulation::HasPoint;
use crate::triangulation::Lies;
use std::marker::PhantomData;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::AtomicUsize;

use serde::{Deserialize, Serialize};

pub fn edge_from_target<T: Clone>(
    target: EdgeTarget,
    triangulation: &SurfaceTriangulation<T>,
) -> Edge<VertexIndex> {
    let edge_ref = triangulation.qeds.edge_a_ref(target);
    edge_ref.edge().clone()
}

pub fn triangle_across_from_target<T>(
    target: EdgeTarget,
    triangulation: &SurfaceTriangulation<T>,
) -> NodeTarget {
    let edge_ref = triangulation.qeds.edge_a_ref(target);
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
        let edge_ref = triangulation.qeds.edge_a_ref(self.into());
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

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
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

/// An offset into the vertices vector.
pub type VertexIndex = usize;

#[derive(Clone, Debug, Serialize, Deserialize)]
/// A Qeds data structure specialised to a 2d triangulation.
pub struct SurfaceTriangulation<T> {
    /// The quad-edge data structure we use as the basis for the triangulation.
    pub qeds: Qeds<VertexIndex, ()>,
    pub vertices: Vec<Segment<T>>,
    pub boundary_edge: EdgeTarget,
    pub bounds: Option<(Point, Point)>,
}

impl<T: Serialize> SurfaceTriangulation<T> {
    pub fn debug_table(&self) -> String {
        use prettytable::*;
        use prettytable::{row, Cell, Row, Table};
        // Create the table
        let mut table = Table::new();
        fn to_letter(i: usize) -> char {
            char::from_u32((b'a' + ((i%26) as u8)).into()).unwrap()
        }

        let mut headers = Row::empty();
        headers.add_cell(Cell::new("-"));
        let mut next = Row::empty();
        next.add_cell(Cell::new("next"));
        let mut rot_next = Row::empty();
        rot_next.add_cell(Cell::new("Rot.next"));
        let mut sym_next = Row::empty();
        sym_next.add_cell(Cell::new("Sym.next"));
        let mut sym_rot_next = Row::empty();
        sym_rot_next.add_cell(Cell::new("SymRot.next"));
        let mut org = Row::empty();
        org.add_cell(Cell::new("org"));
        let mut dest = Row::empty();
        dest.add_cell(Cell::new("dest"));
        for (i, quad) in self.qeds.quads.iter() {
            headers.add_cell(Cell::new(&format!("{}", to_letter(i))));
            {
                let target = quad.edges_a[0].next;
                let mut s = format!("{}", to_letter(target.e));
                if target.r == 1 {
                    s.push_str("Rot")
                } else if target.r == 2 {
                    s.push_str("Sym")
                } else if target.r == 3 {
                    s.push_str("SymRot")
                }
                next.add_cell(Cell::new(&s));
            }
            {
                let target = quad.edges_b[0].next;
                let mut s = format!("{}", to_letter(target.e));
                if target.r == 1 {
                    s.push_str("Rot")
                } else if target.r == 2 {
                    s.push_str("Sym")
                } else if target.r == 3 {
                    s.push_str("SymRot")
                }
                rot_next.add_cell(Cell::new(&s));
            }
            {
                let target = quad.edges_a[1].next;
                let mut s = format!("{}", to_letter(target.e));
                if target.r == 1 {
                    s.push_str("Rot")
                } else if target.r == 2 {
                    s.push_str("Sym")
                } else if target.r == 3 {
                    s.push_str("SymRot")
                }
                sym_next.add_cell(Cell::new(&s));
            }
            {
                let target = quad.edges_b[1].next;
                let mut s = format!("{}", to_letter(target.e));
                if target.r == 1 {
                    s.push_str("Rot")
                } else if target.r == 2 {
                    s.push_str("Sym")
                } else if target.r == 3 {
                    s.push_str("SymRot")
                }
                sym_rot_next.add_cell(Cell::new(&s));
            }
            org.add_cell(Cell::new(&format!("P{}",quad.edges_a[0].point)));
            dest.add_cell(Cell::new(&format!("P{}",quad.edges_a[1].point)));
        }
        table.add_row(headers);
        table.add_row(next);
        table.add_row(rot_next);
        table.add_row(sym_next);
        table.add_row(sym_rot_next);
        table.add_row(org);
        table.add_row(dest);
        let mut edge_table = table.to_string();
        let points_table = self.debug_points_table();
        edge_table.push_str("\n");
        edge_table.push_str(&points_table);
        edge_table
    }
    pub fn debug_points_table(&self) -> String {
        use prettytable::*;
        use prettytable::{row, Cell, Row, Table};
        // Create the table
        let mut table = Table::new();
        fn to_letter(i: usize) -> char {
            char::from_u32((b'a' + ((i%26) as u8)).into()).unwrap()
        }

        let mut headers = Row::empty();
        headers.add_cell(Cell::new("Name"));
        headers.add_cell(Cell::new("Point"));
        table.add_row(headers);
        for (i, segment) in self.vertices.iter().enumerate() {
            let mut point_row = Row::empty();
            point_row.add_cell(Cell::new(&format!("P{}", i)));
            point_row.add_cell(Cell::new(&format!("{}", segment.point)));
            table.add_row(point_row);
        }
        table.to_string()
    }
    pub fn debug_dump(&self, msg: Option<&str>) {
        static N: AtomicUsize = AtomicUsize::new(0);
        let old_n = N.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        // Only do this the first 20 times, and don't do it during tests
        #[cfg(not(test))]
        if old_n < 20 {
            let js = serde_json::to_string_pretty(self).unwrap();
            let table = self.debug_table();
            // let filename = if old_n % 2 == 0 { "a.json" } else { "b.json" };
            let filename = if let Some(msg) = msg {
                format!("{} - {}", old_n, msg)
            } else {
                format!("{}", old_n)
            };
            eprintln!("outputting: {}", filename);
            std::fs::write(format!("{}.json", filename), js).unwrap();
            std::fs::write(format!("{}.txt", filename), table).unwrap();
        }
    }
}
impl<T> SurfaceTriangulation<T> {
    fn is_tri_real(&self, edge_ref: EdgeRefA<'_, VertexIndex, ()>) -> bool {
        let first = edge_ref;
        let mut current = first;
        let mut i = 0;
        while i < 3 {
            // TODO: this is very inefficient, fix BData
            if !self.curves_left(current) {
                return false;
            }
            current = current.l_next();
            if current == first {
                break;
            }
            i += 1;
        }
        true
    }
    fn lies_right(&self, edge_ref: EdgeRefA<'_, VertexIndex, ()>, point: Point) -> Lies {
        use Lies::*;
        let pa = self.vertices.get(edge_ref.edge().point).unwrap().point();
        let pb = self
            .vertices
            .get(edge_ref.sym().edge().point)
            .unwrap()
            .point();
        // We have special treatment of infinite edges.
        match left_or_right(pa, pb, point) {
            Direction::Right => Yes,
            Direction::Straight => On,
            Direction::Left => No,
        }
    }

    pub fn lies_right_strict(&self, edge_ref: EdgeRefA<'_, VertexIndex, ()>, point: Point) -> bool {
        self.lies_right(edge_ref, point) == Lies::Yes
    }

    pub fn is_outward_boundary(&self, edge_ref: EdgeRefA<'_, VertexIndex, ()>) -> bool {
        // edge_ref.rot().edge().point == Space::Out || edge_ref.rot().sym().edge().point == Space::Out
        !self.curves_left(edge_ref)
    }

    pub fn is_boundary(&self, edge_ref: EdgeRefA<'_, VertexIndex, ()>) -> bool {
        self.is_outward_boundary(edge_ref) || self.is_outward_boundary(edge_ref.sym())
    }

    pub fn curves_left(&self, e: EdgeRefA<'_, VertexIndex, ()>) -> bool {
        let org = self.vertices.get(e.edge().point).unwrap().point;
        let dest = self.vertices.get(e.sym().edge().point).unwrap().point;
        let next = self
            .vertices
            .get(e.l_next().sym().edge().point)
            .unwrap()
            .point;
        is_ccw(org, dest, next)
    }

    pub fn lies_left_strict(&self, edge_ref: EdgeRefA<'_, VertexIndex, ()>, point: Point) -> bool {
        self.lies_left(edge_ref, point) == Lies::Yes
    }

    pub fn lies_right_or_on(&self, edge_ref: EdgeRefA<'_, VertexIndex, ()>, point: Point) -> bool {
        self.lies_right(edge_ref, point) != Lies::No
    }

    pub fn lies_left(&self, edge_ref: EdgeRefA<'_, VertexIndex, ()>, point: Point) -> Lies {
        let pa = self.vertices.get(edge_ref.edge().point).unwrap().point();
        let pb = self
            .vertices
            .get(edge_ref.sym().edge().point)
            .unwrap()
            .point();
        // We have special treatment of infinite edges.
        match left_or_right(pa, pb, point) {
            Direction::Right => Lies::No,
            Direction::Straight => Lies::On,
            Direction::Left => Lies::Yes,
        }
    }

    pub fn lies_on(&self, edge_ref: EdgeRefA<'_, VertexIndex, ()>, point: Point) -> bool {
        let pa = self.vertices.get(edge_ref.edge().point).unwrap().point();
        let pb = self
            .vertices
            .get(edge_ref.sym().edge().point)
            .unwrap()
            .point();
        // We have special treatment of infinite edges.
        let res_a = match left_or_right(pa, pb, point) {
            Direction::Right => false,
            Direction::Straight => true,
            Direction::Left => false,
        };
        let res_b = match left_or_right(pb, pa, point) {
            Direction::Right => false,
            Direction::Straight => true,
            Direction::Left => false,
        };
        // TODO: we shouldn't have to do this
        res_a || res_b
    }

    pub fn in_left_face(&self, edge_ref: EdgeRefA<'_, VertexIndex, ()>, point: Point) -> bool {
        for edge in edge_ref.l_face().edges.iter() {
            let is_left = match self.lies_left(*edge, point) {
                Lies::No => false,
                Lies::On => false,
                Lies::Yes => true,
            };
            if !is_left {
                return false;
            }
        }
        true
    }

    pub fn on_left_face(&self, edge_ref: EdgeRefA<'_, VertexIndex, ()>, point: Point) -> bool {
        for edge in edge_ref.l_face().edges.iter() {
            let is_left = match self.lies_left(*edge, point) {
                Lies::No => false,
                Lies::On => true,
                Lies::Yes => true,
            };
            if !is_left {
                return false;
            }
        }
        true
    }
}

impl<T: Default + Clone + Serialize> SurfaceTriangulation<T> {
    pub fn add_to_quad_unchecked(
        &mut self,
        mut edge_target: EdgeTarget,
        point: Point,
        data: T,
    ) -> EdgeTarget {
        let first_index = self.qeds.edge_a_ref(edge_target).edge().point;
        let new_index = self.vertices.len();
        self.vertices.push(Segment::new(point, data));
        let mut base = self.qeds.make_edge_with_a(first_index, new_index).target();
        self.qeds.splice(base, edge_target);
        loop {
            let base_ref = self.connect(edge_target, base.sym());
            edge_target = base_ref.oprev().target();
            base = base_ref.target();
            if self.qeds.edge_a(edge_target.sym()).point == first_index {
                break;
            }
        }
        let e = self.qeds.edge_a_ref(base).oprev().target();
        debug_assert_spaces(self);
        self.retriangulate_suspect_edges(e, point, first_index);
        debug_assert_eq!(
            self.vertices
                .get(self.qeds.edge_a_ref(base.sym()).edge().point)
                .unwrap()
                .point(),
            point
        );
        // self.retriangulate_all();
        // for fail in self.fail_del_test() {
        //     eprintln!("{:?} failed del test", fail);
        // }
        // debug_assert_eq!(0, self.n_fail_del_test());
        // eprintln!("inserted: {}", point);
        debug_assert_spaces(self);
        base.sym()
    }

    fn connect(&'_ mut self, a: EdgeTarget, b: EdgeTarget) -> EdgeRefA<'_, VertexIndex, ()> {
        let e = self.qeds.connect(a, b).target();
        self.qeds.edge_a_ref(e)
    }

    pub fn add_point_with_default(&mut self, point: Point) -> Option<EdgeTarget> {
        self.add_point(point, Default::default())
    }

    /// The edge this returns should always have the added point at its origin.
    /// It should not result in non-CCW triangles.
    pub fn add_point(&mut self, mut point: Point, data: T) -> Option<EdgeTarget> {
        self.add_point_raw(point, data, false)
    }

    pub fn add_point_force(&mut self, mut point: Point, data: T) -> Option<EdgeTarget> {
        self.add_point_raw(point, data, true)
    }
    fn add_point_raw(&mut self, mut point: Point, data: T, force: bool) -> Option<EdgeTarget> {
        point.snap();
        self.update_bounds(point);
        debug_assert_spaces(self);
        // assert_eq!(0, self.retriangulate_all());
        match self.locate_raw(point, force)? {
            Location::OnPoint(edge) => Some(edge.target()),
            Location::OnEdge(edge) => {
                let target = edge.target();
                Some(self.add_point_to_edge_unchecked(target, point, data))
            }
            Location::OnFace(edge) => {
                let (pa, pb) = self.get_segment(edge);
                let ccw = is_ccw(pa, pb, point);
                if !ccw {
                    eprintln!(
                        "left or right {} {} {}: {:?}",
                        pa,
                        pb,
                        point,
                        left_or_right(pa, pb, point)
                    );
                    assert!(ccw);
                }
                let target = edge.target();
                Some(self.add_to_quad_unchecked(target, point, data))
            }
        }
    }

    /// Same as [`add_point_to_edge`] but does not check if the point is on one
    /// of the vertices of the edge.
    fn add_point_to_edge_unchecked(
        &mut self,
        mut edge_target: EdgeTarget,
        point: Point,
        data: T,
    ) -> EdgeTarget {
        // eprintln!("adding to edge");
        let is_boundary = self.is_boundary(self.qeds.edge_a_ref(edge_target));
        // let is_boundary = false;
        // eprintln!("is boundary: {}", is_boundary);
        // if is_boundary {
        //     edge_target = edge_target.sym();
        // }
        // TODO: the insertion algorithm doesn't properly handle edge cases.
        if is_boundary {
            if self.is_outward_boundary(self.qeds.edge_a_ref(edge_target)) {
                edge_target = edge_target.sym();
            }
            self.add_to_boundary_unchecked(edge_target, point, data)
        } else {
            let oprev = self.qeds.edge_a_ref(edge_target).oprev().target();
            self.qeds.delete(edge_target);
            edge_target = oprev;
            self.add_to_quad_unchecked(edge_target, point, data)
        }
    }

    fn add_to_boundary_unchecked(
        &mut self,
        mut edge_target: EdgeTarget,
        point: Point,
        data: T,
    ) -> EdgeTarget {
        debug_assert_spaces(self);
        let x_onext = self.qeds.edge_a_ref(edge_target).onext().target();
        let x_dprev = self.qeds.edge_a_ref(edge_target).d_prev().target();
        debug_assert_ne!(x_onext, edge_target);
        debug_assert_ne!(x_dprev, edge_target);
        // let dia_a_indices = self.get_segment_indices(self.qeds.edge_a_ref(edge_target).d_prev().sym());
        self.debug_dump(Some("Before delete"));
        self.qeds.delete(edge_target);
        self.debug_dump(Some("Delete[c]"));
        // let dia_c = self.qeds.edge_a_ref(edge_target);
        // let dia_a = dia_c.sym().onext().target();
        // let dia_e = dia_c.l_next().target();
        // assert_eq!(dia_e,dia_a.sym());
        // edge_target = onext;

        let first_index = self.qeds.edge_a_ref(x_onext).edge().point;
        eprintln!("first index: {}", first_index);
        let new_index = self.vertices.len();
        self.vertices.push(Segment::new(point, data));
        let mut base = self.qeds.make_edge_with_a(first_index, new_index).target();
        // assert_eq!()
        self.debug_dump(Some("MakeEdge"));
        self.qeds.splice(base, x_onext);
        self.debug_dump(Some("Splice[d,c.Dprev]"));
        {
            let base_ref = self.connect(base,x_onext.sym());
            let base_sym = base_ref.target().sym();
            self.debug_dump(Some("Connect[c.Dprev.Sym,d]"));
            let base_ref = self.connect(base_sym,x_dprev.sym());
            self.debug_dump(Some("Connect[e.Oprev.Sym,e.Sym]"));
            // edge_target = base_ref.onext().target();
        }
        // debug_assert_eq!(
        //     self.vertices
        //         .get(self.qeds.edge_a_ref(base.sym()).edge().point)
        //         .unwrap()
        //         .point(),
        //     point
        // );
        // let e = self.qeds.edge_a_ref(base).oprev().target();
        debug_assert_spaces(self);
        // self.retriangulate_suspect_edges(edge_target, point, first_index);
        // TODO: remove this.
        self.retriangulate_all();
        // assert_eq!(0, self.retriangulate_all());
        // self.retriangulate_all();
        // for fail in self.fail_del_test() {
        //     eprintln!("{:?} failed del test", fail);
        // }
        // debug_assert_eq!(0, self.n_fail_del_test());
        // eprintln!("inserted: {}", point);
        // debug_assert_spaces(self);
        base.sym()
    }

    /// Add a point to a specified edge. If the point lies on one of the
    /// vertices just add it there.
    pub fn add_point_to_edge(
        &mut self,
        edge_target: EdgeTarget,
        point: Point,
        data: T,
    ) -> EdgeTarget {
        {
            let point_a = self
                .vertices
                .get(self.qeds.edge_a_ref(edge_target).edge().point)
                .unwrap()
                .point;
            let point_b = self
                .vertices
                .get(self.qeds.edge_a_ref(edge_target).sym().edge().point)
                .unwrap()
                .point;

            if point_a == point {
                edge_target
            } else if point_b == point {
                edge_target.sym()
            } else {
                self.add_point_to_edge_unchecked(edge_target, point, data)
            }
        }
    }
    pub fn get_matching_edge_indices(
        &self,
        i1: VertexIndex,
        i2: VertexIndex,
    ) -> Option<EdgeRefA<'_, VertexIndex, ()>> {
        for edge in self.base_targets() {
            let edge_ref = self.qeds.edge_a_ref(edge);
            let edge_ref_sym = self.qeds.edge_a_ref(edge.sym());
            let x1 = edge_ref.edge().point;
            let x2 = edge_ref_sym.edge().point;
            if [x1, x2] == [i1, i2] {
                return Some(edge_ref);
            }
            if [x1, x2] == [i2, i1] {
                return Some(edge_ref_sym);
            }
        }
        None
    }
    pub fn get_matching_edge(&self, p1: Point, p2: Point) -> Option<EdgeRefA<'_, VertexIndex, ()>> {
        for edge in self.base_targets() {
            let edge_ref = self.qeds.edge_a_ref(edge);
            let segment = self.get_segment(edge_ref);
            if segment == (p1, p2) {
                return Some(edge_ref);
            }
            if segment == (p2, p1) {
                return Some(edge_ref.sym());
            }
        }
        None
    }
}

impl<T: Default> SurfaceTriangulation<T> {
    pub fn new_with_default(a: Point, b: Point, c: Point) -> Self {
        assert!(is_ccw(a, b, c));
        let seg_a = Segment::new(a, Default::default());
        let seg_b = Segment::new(b, Default::default());
        let seg_c = Segment::new(c, Default::default());
        let a_index = 0;
        let b_index = 1;
        let c_index = 2;
        let mut qeds = Qeds::new();
        let edge_a = qeds.make_edge_with_ab(a_index, b_index, (), ()).target();
        let edge_b = qeds.make_edge_with_ab(b_index, c_index, (), ()).target();
        let result = {
            qeds.splice(edge_a.sym(), edge_b);
            qeds.connect(edge_b, edge_a);
            // qeds.quads[2].edges_b[0].point = Space::Out;
            SurfaceTriangulation {
                qeds,
                vertices: vec![seg_a, seg_b, seg_c],
                boundary_edge: edge_a,
                bounds: None,
            }
        };
        result
    }
}

impl<T: Clone> SurfaceTriangulation<T> {
    pub fn new(a: (Point, T), b: (Point, T), c: (Point, T)) -> Self {
        assert!(is_ccw(a.0, b.0, c.0));
        let seg_a = Segment::new(a.0, a.1);
        let seg_b = Segment::new(b.0, b.1);
        let seg_c = Segment::new(c.0, c.1);
        let a_index = 0;
        let b_index = 1;
        let c_index = 2;
        let mut qeds = Qeds::new();
        let edge_a = qeds.make_edge_with_ab(a_index, b_index, (), ()).target();
        let edge_b = qeds.make_edge_with_ab(b_index, c_index, (), ()).target();
        let result = {
            qeds.splice(edge_a.sym(), edge_b);
            qeds.connect(edge_b, edge_a);
            // qeds.quads[2].edges_b[0].point = Space::Out;
            SurfaceTriangulation {
                qeds,
                vertices: vec![seg_a, seg_b, seg_c],
                boundary_edge: edge_a,
                bounds: None,
            }
        };
        result
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
    fn retriangulate_suspect_edges(&mut self, mut e: EdgeTarget, point: Point, first_i: usize) {
        // The suspect edges are e(.Onext.Lprev)^k for k=0,1,2,3...
        loop {
            let t = self.qeds.edge_a_ref(e).oprev().target();
            let t_dest = self
                .vertices
                .get(self.qeds.edge_a_ref(t).sym().edge().point)
                .unwrap()
                .point;
            let e_dest = self
                .vertices
                .get(self.qeds.edge_a_ref(e).sym().edge().point)
                .unwrap()
                .point;
            let e_org_i = self.qeds.edge_a_ref(e).edge().point;
            let e_org = self.vertices.get(e_org_i).unwrap().point;
            if self.lies_right_strict(self.qeds.edge_a_ref(e), t_dest)
                && del_test_ccw(e_org, t_dest, e_dest, point)
            // && !self.is_boundary(self.qeds.edge_a_ref(e))
            {
                self.swap(e);
                // assert!(self.del_test(e)||self.cocircular(e));
                e = self.qeds.edge_a_ref(t).target();
            } else if e_org_i == first_i {
                break;
            } else {
                e = self.qeds.edge_a_ref(e).onext().l_prev().target();
            }
        }
    }

    fn n_fail_del_test(&self) -> usize {
        let mut n_fails = 0;
        let edge_targets = self.qeds.base_edges().map(|edge| edge.target());
        for e in edge_targets {
            {
                if self.del_test(e) && self.concave_test(e) && !self.cocircular(e) {
                    n_fails += 1;
                }
            }
        }
        n_fails
    }

    fn fail_del_test(&self) -> Vec<EdgeTarget> {
        let mut fails = vec![];
        let edge_targets = self.qeds.base_edges().map(|edge| edge.target());
        for e in edge_targets {
            {
                if self.del_test(e) && self.concave_test(e) {
                    fails.push(e);
                }
            }
        }
        fails
    }

    /// Warning: this is very inefficient and just for testing.
    fn retriangulate_all_single_pass(&mut self) -> usize {
        let mut swaps = 0;
        #[allow(clippy::needless_collect)]
        let edge_targets: Vec<EdgeTarget> =
            self.qeds.base_edges().map(|edge| edge.target()).collect();
        for e in edge_targets.into_iter() {
            if self.del_test(e) && self.concave_test(e) {
                swaps += 1;
                self.swap(e);
            }
        }
        swaps
    }

    /// Perform Delaunay swapping on the entire triangulation until complete.
    /// Should not be necessary, mainly included for testing.
    pub fn retriangulate_all(&mut self) -> usize {
        let mut iterations = 0;
        let mut total_swaps = 0;
        loop {
            if iterations > 100 {
                // panic!("too many triangulation iterations");
                break total_swaps;
            }
            let swaps = self.retriangulate_all_single_pass();
            total_swaps += swaps;
            if swaps == 0 {
                break total_swaps;
            }
            iterations += 1;
        }
    }

    // fn get_outside(&self) -> Option<EdgeTarget> {
    //     // Find at least one edge with an Out property.
    //     for (i, quad) in self.qeds.quads.iter() {
    //         if quad.edges_b[0].point == Space::Out {
    //             return Some(EdgeTarget::new(i, 1, 0));
    //         }
    //         if quad.edges_b[1].point == Space::Out {
    //             return Some(EdgeTarget::new(i, 3, 0));
    //         }
    //     }
    //     None
    // }

    pub fn neighbouring_points<'a>(
        &'a self,
        edge: EdgeRefA<'a, VertexIndex, ()>,
    ) -> Vec<EdgeRefA<'a, VertexIndex, ()>> {
        let mut neighbouring_refs = vec![edge.sym()];
        let mut current_edge = edge;
        loop {
            current_edge = current_edge.onext();
            if current_edge == edge {
                break;
            }
            neighbouring_refs.push(current_edge.sym());
        }
        neighbouring_refs
    }

    pub fn locate(&self, mut point: Point) -> Option<Location<'_>> {
        self.locate_raw(point, false)
    }

    pub fn locate_force(&self, mut point: Point) -> Option<Location<'_>> {
        self.locate_raw(point, true)
    }

    // TODO: This algorithm is known to fail in non-Delaunay triangulations.
    fn locate_raw(&self, mut point: Point, force: bool) -> Option<Location<'_>> {
        point.snap();
        use rand::Rng;
        let mut e = self.some_edge_a().unwrap();
        let mut rng = rand::thread_rng();
        let mut current_iterations = 0;
        let location: Location<'_> = loop {
            current_iterations += 1;
            if current_iterations > 200 {
                // for (a, b, c) in self.triangles() {
                //     println!("Triangle: {:?}", (a.point, b.point, c.point));
                // }
                for vertex in self.vertices.iter() {
                    eprintln!("{}", vertex.point);
                }
                // #[cfg(debug_assertions)]
                // panic!("locating failed for: {}", point);
                return None;
            }
            if point == self.vertices.get(e.edge().point).unwrap().point {
                break Location::OnPoint(e);
            } else if point == self.vertices.get(e.sym().edge().point).unwrap().point {
                // If we lie on the eDest, return eSym(). This is a departure
                // from the paper which just returns e.
                break Location::OnPoint(e.sym());
            }
            // A random variable to determine if onext is tested first. If not
            // it is tested second.
            let onext_first: bool = rng.gen();
            // let segment = self.get_segment(e);
            // eprintln!(
            //     "checking: {}-{} Left/Right: {:?}",
            //     segment.0,
            //     segment.1,
            //     left_or_right(segment.0, segment.1, point)
            // );
            let has_left_face = self.curves_left(e);
            // If it liest not_right to a boundary edge it must be either out of bounds or on this triangular edge
            if !has_left_face && !self.lies_right_strict(e, point) {
                if self.lies_right_strict(e.sym(), point) {
                    // Is out of bounds
                    return None;
                } else {
                    // Find the right edge
                    let mut boundary_edge = e;
                    loop {
                        let (p1, p2) = self.get_segment(boundary_edge);
                        let edge_length = p1.distance(p2);
                        let p1d = p1.distance(point);
                        let p2d = p2.distance(point);
                        if (p1d + p2d - edge_length).abs() < edge_length * 0.001 {
                            return Some(Location::OnEdge(boundary_edge.sym()));
                        }
                        boundary_edge = boundary_edge.l_next();
                        if boundary_edge == e {
                            // If we've looped all the way around, we can't
                            // find a match and return None.
                            return None;
                        }
                    }
                }
            }
            if self.lies_right_strict(e, point) {
                // eprintln!(
                //     "point {} lies strictly right of {}-{}",
                //     point, segment.0, segment.1
                // );
                e = e.sym();
            } else if has_left_face {
                // eprintln!("e HasLeftFace: {:?}", has_left_face);
                #[allow(clippy::collapsible_else_if)]
                if onext_first {
                    if !self.lies_right_strict(e.onext(), point) {
                        e = e.onext();
                    } else if !self.lies_right_strict(e.d_prev(), point) {
                        e = e.d_prev();
                    } else {
                        if self.lies_left_strict(e, point) {
                            break Location::OnFace(e);
                        } else {
                            break Location::OnEdge(e);
                        }
                    }
                } else {
                    if !self.lies_right_strict(e.d_prev(), point) {
                        e = e.d_prev();
                    } else if !self.lies_right_strict(e.onext(), point) {
                        e = e.onext();
                    } else {
                        if self.lies_left_strict(e, point) {
                            break Location::OnFace(e);
                        } else {
                            break Location::OnEdge(e);
                        }
                    }
                }
            } else if self.lies_left_strict(e, point) {
                // Lies out of bounds
                return None;
            } else {
                break Location::OnEdge(e.sym());
            }
        };
        Some(location)
    }
    /// Return the canonical tri within which this point is located.
    pub fn locate_tri(&self, point: Point) -> Option<EdgeRefA<'_, VertexIndex, ()>> {
        Some(self.locate(point)?.edge().get_tri_canonical())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Location<'a> {
    OnPoint(EdgeRefA<'a, VertexIndex, ()>),
    OnEdge(EdgeRefA<'a, VertexIndex, ()>),
    OnFace(EdgeRefA<'a, VertexIndex, ()>),
}

impl<'a> Location<'a> {
    pub fn edge(&self) -> EdgeRefA<'a, VertexIndex, ()> {
        match *self {
            Location::OnPoint(edge) => edge,
            Location::OnEdge(edge) => edge,
            Location::OnFace(edge) => edge,
        }
    }
}

pub struct SegmentsIter<'a, T> {
    slab_iter: slab::Iter<'a, Quad<VertexIndex, ()>>,
    current_quad: Option<Quad<VertexIndex, ()>>,
    _data: PhantomData<T>,
}

impl<'a, T> SegmentsIter<'a, T> {
    pub fn new(triangulation: &'a SurfaceTriangulation<T>) -> Self {
        Self {
            slab_iter: triangulation.qeds.quads.iter(),
            current_quad: None,
            _data: PhantomData,
        }
    }
}

impl<'a, T: Clone> Iterator for SegmentsIter<'a, T> {
    type Item = Edge<VertexIndex>;
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(current_quad) = &self.current_quad {
            let edge = current_quad.edges_a[1].clone();
            self.current_quad = None;
            Some(edge)
        } else if let Some((_, quad)) = self.slab_iter.next() {
            let edge = quad.edges_a[0].clone();
            self.current_quad = Some(quad.clone());
            Some(edge)
        } else {
            None
        }
    }
}

pub struct SegmentsIterMut<'a, T> {
    vertex_iter: std::slice::IterMut<'a, Segment<T>>,
}

impl<'a, T> SegmentsIterMut<'a, T> {
    pub fn new(triangulation: &'a mut SurfaceTriangulation<T>) -> Self {
        Self {
            vertex_iter: triangulation.vertices.iter_mut(),
        }
    }
}

impl<'a, T: Clone> Iterator for SegmentsIterMut<'a, T> {
    type Item = &'a mut Segment<T>;
    fn next(&mut self) -> Option<Self::Item> {
        self.vertex_iter.next()
    }
}

impl<T> SurfaceTriangulation<T> {
    pub fn segments(&self) -> SegmentsIter<'_, T> {
        SegmentsIter::new(self)
    }
    pub fn base_targets(&self) -> impl Iterator<Item = EdgeTarget> + '_ {
        self.qeds
            .quads
            .iter()
            .map(|(i, _)| EdgeTarget::new(i, 0, 0))
    }
    pub fn segments_mut(&mut self) -> SegmentsIterMut<'_, T> {
        SegmentsIterMut::new(self)
    }
    pub fn raw_triangles(&self) -> RawTriangleIter<'_, T> {
        RawTriangleIter::new(self)
    }
    pub fn triangles(&self) -> TriangleIter<'_, T> {
        TriangleIter::new(self)
    }
    pub fn triangle_edges(&self) -> TriangleEdgeIter<'_, T> {
        TriangleEdgeIter::new(self)
    }
    // pub fn triangles_mut(&mut self) -> TriangleIterMut<T> {
    //     TriangleIterMut::new(self)
    // }
    pub fn nodes(&self) -> NodeIter<'_, T> {
        NodeIter::new(self)
    }

    pub fn qeds(&self) -> Option<&Qeds<VertexIndex, ()>> {
        Some(&self.qeds)
    }

    pub fn some_edge_a(&self) -> Option<EdgeRefA<'_, VertexIndex, ()>> {
        let (i, _) = self.qeds.quads.iter().next()?;
        Some(self.qeds.edge_a_ref(EdgeTarget::new(i, 0, 0)))
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

    pub fn boundary(&self) -> BoundaryIter<'_, VertexIndex, ()> {
        BoundaryIter::new(&self.qeds, self.boundary_edge)
    }

    /// Determine whether an Edge (i.e. two NavTris) satisfies the Delaunay
    /// criterion. An edge with only a single NavTri will return True.
    fn del_test(&self, e: EdgeTarget) -> bool {
        // Get the edge.
        let edge = self.qeds.edge_a_ref(e);
        // Get all of the vertices around this egdge in a CCW order.
        let a = self.vertices.get(edge.oprev().edge().point).unwrap().point;
        let b = self
            .vertices
            .get(edge.r_prev().sym().edge().point)
            .unwrap()
            .point;
        let c = self.vertices.get(edge.l_next().edge().point).unwrap().point;
        let d = self
            .vertices
            .get(edge.onext().sym().edge().point)
            .unwrap()
            .point;
        // eprintln!("del_test edge: {:?}", self.get_segment(edge));
        // eprintln!("del_test points: {:?}", [a, b, c, d]);
        del_test_ccw(a, b, c, d)
    }
    fn cocircular(&self, e: EdgeTarget) -> bool {
        // Get the edge.
        let edge = self.qeds.edge_a_ref(e);
        // Get all of the vertices around this egdge in a CCW order.
        let a = self.vertices.get(edge.oprev().edge().point).unwrap().point;
        let b = self
            .vertices
            .get(edge.r_prev().sym().edge().point)
            .unwrap()
            .point;
        let c = self.vertices.get(edge.l_next().edge().point).unwrap().point;
        let d = self
            .vertices
            .get(edge.onext().sym().edge().point)
            .unwrap()
            .point;
        cocircular(a, b, c, d)
    }
    pub fn get_segment(&self, edge: EdgeRefA<'_, VertexIndex, ()>) -> (Point, Point) {
        let i1 = edge.edge().point;
        let i2 = edge.sym().edge().point;
        debug_assert_ne!(i1, i2);
        let p1 = self.vertices.get(i1).unwrap().point;
        let p2 = self.vertices.get(i2).unwrap().point;
        // debug_assert_ne!(p1, p2);
        (p1, p2)
    }
    /// Test whether an edge is the diagonal of a concave quad.
    fn concave_test(&self, e: EdgeTarget) -> bool {
        // Get the edge.
        let edge = self.qeds.edge_a_ref(e);
        // Get all of the vertices around this egdge in a CCW order.
        let a = self.vertices.get(edge.oprev().edge().point).unwrap().point;
        let b = self
            .vertices
            .get(edge.r_prev().sym().edge().point)
            .unwrap()
            .point;
        let c = self.vertices.get(edge.l_next().edge().point).unwrap().point;
        let d = self
            .vertices
            .get(edge.onext().sym().edge().point)
            .unwrap()
            .point;
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

impl<'a> EdgeRefA<'a, VertexIndex, ()> {
    pub fn triangle_across(&self) -> Self {
        let edge_on_next_tri = self.sym();
        edge_on_next_tri.get_tri_canonical()
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
    type Item = EdgeRefA<'a, VertexIndex, ()>;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // First we check that the edge actually exists for this given
            // EdgeTarget.
            if self.triangulation.qeds.quads.contains(self.next.e) {
                let edge_ref = { self.triangulation.qeds.edge_a_ref(self.next) };
                self.inc();
                if edge_ref.is_tri_canonical() && self.triangulation.is_tri_real(edge_ref) {
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
            self.raw_triangles
                .triangulation
                .vertices
                .get(tri_a.edge().point)
                .unwrap()
                .point,
            self.raw_triangles
                .triangulation
                .vertices
                .get(tri_b.edge().point)
                .unwrap()
                .point,
            self.raw_triangles
                .triangulation
                .vertices
                .get(tri_c.edge().point)
                .unwrap()
                .point,
        );
        let a = self
            .raw_triangles
            .triangulation
            .vertices
            .get(tri_a.edge().point)
            .unwrap()
            .clone();
        let b = self
            .raw_triangles
            .triangulation
            .vertices
            .get(tri_b.edge().point)
            .unwrap()
            .clone();
        let c = self
            .raw_triangles
            .triangulation
            .vertices
            .get(tri_c.edge().point)
            .unwrap()
            .clone();
        if !ccw {
            println!("a: {}", a.point);
            println!("b: {}", b.point);
            println!("c: {}", c.point);
        }
        // assert!(ccw);
        Some((a, b, c))
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
    type Item = [(EdgeTarget, Segment<T>); 3];
    fn next(&mut self) -> Option<Self::Item> {
        let tri_a = self.raw_triangles.next()?;
        let tri_b = tri_a.l_next();
        let tri_c = tri_b.l_next();
        assert_eq!(tri_c.l_next().target, tri_a.target);
        Some([
            (
                tri_a.target,
                self.raw_triangles
                    .triangulation
                    .vertices
                    .get(tri_a.edge().point)
                    .unwrap()
                    .clone(),
            ),
            (
                tri_b.target,
                self.raw_triangles
                    .triangulation
                    .vertices
                    .get(tri_b.edge().point)
                    .unwrap()
                    .clone(),
            ),
            (
                tri_c.target,
                self.raw_triangles
                    .triangulation
                    .vertices
                    .get(tri_c.edge().point)
                    .unwrap()
                    .clone(),
            ),
        ])
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
                let edge_ref = { self.triangulation.qeds.edge_a_ref(self.next) };
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

impl<'a> EdgeRefA<'a, VertexIndex, ()> {
    /// An edge is the canonical edge for a tri iff it has the lowest e value
    /// for the tri.
    fn is_tri_canonical(&self) -> bool {
        let current_e = self.target().e;
        let next_edge = self.sym_rot().onext();
        let next_e = next_edge.target().e;
        let second_e = next_edge.onext().target().e;
        (current_e < next_e) && (current_e < second_e)
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

pub fn has_edge<T>(triangulation: &SurfaceTriangulation<T>, pa: Point, pb: Point) -> bool {
    get_edge(triangulation, pa, pb).is_some()
}

pub fn get_edge<T>(
    triangulation: &SurfaceTriangulation<T>,
    pa: Point,
    pb: Point,
) -> Option<EdgeRefA<'_, VertexIndex, ()>> {
    for (i, quad) in triangulation.qeds().unwrap().quads.iter() {
        let edge1 = triangulation
            .vertices
            .get(quad.edges_a[0].point)
            .unwrap()
            .point;
        let edge2 = triangulation
            .vertices
            .get(quad.edges_a[1].point)
            .unwrap()
            .point;
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
                        triangulation
                            .vertices
                            .get(quad.edges_a[0].point)
                            .unwrap()
                            .point,
                        triangulation
                            .vertices
                            .get(quad.edges_a[1].point)
                            .unwrap()
                            .point,
                        quad.edges_b[0].point,
                        quad.edges_b[1].point
                    );
                }
                triangulation.add_point(centroid, (0.0, 0.0));
                for (i, quad) in triangulation.qeds.quads.iter() {
                    println!(
                        "AfterQuad[{}]: {}-{} {:?}-{:?}",
                        i,
                        triangulation
                            .vertices
                            .get(quad.edges_a[0].point)
                            .unwrap()
                            .point,
                        triangulation
                            .vertices
                            .get(quad.edges_a[1].point)
                            .unwrap()
                            .point,
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
            for (a, _b, _c) in
                triangulation
                    .triangle_edges()
                    .map(|[(ta, edge_a), (tb, edge_b), (tc, edge_c)]| {
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
            {
                let oprev = triangulation.qeds.edge_a_ref(edge_target).oprev().target();

                triangulation.qeds.delete(edge_target);
                let mut edge_target = oprev;
                let first_index = triangulation.qeds.edge_a_ref(edge_target).edge().point;
                let new_index = triangulation.vertices.len();
                triangulation.vertices.push(Segment::new(point, (0.0, 0.0)));
                let mut base = triangulation
                    .qeds
                    .make_edge_with_a(first_index, new_index)
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
                    if triangulation.qeds.edge_a(edge_target.sym()).point == first_index {
                        break;
                    }
                }
                debug_assert_eq!(
                    triangulation
                        .vertices
                        .get(triangulation.qeds.edge_a_ref(return_value).edge().point)
                        .unwrap()
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
        debug_assert_spaces(&triangulation);
        assert_eq!(triangulation.triangles().count(), 2);
    }

    #[test]
    fn surface_neighbours() {
        let v1 = Point::new(0.0, 0.0);
        let v2 = Point::new(5.0, 0.0);
        let v3 = Point::new(5.0, 5.0);
        let triangulation = SurfaceTriangulation::new((v1, 0.0), (v2, 0.0), (v3, 0.0));
        assert_eq!(triangulation.triangles().count(), 1);
        debug_assert_spaces(&triangulation);
        let start_edge = triangulation.some_edge_a().unwrap();
        let neighbours = triangulation.neighbouring_points(start_edge);
        assert_eq!(neighbours.len(), 2);
    }

    #[test]
    fn edge_location() {
        let v1 = Point::new(0.0, 0.0);
        let v2 = Point::new(5.0, 0.0);
        let v3 = Point::new(5.0, 5.0);
        let v4 = Point::new(2.5, 0.0);
        let v5 = Point::new(3.0, 0.0);
        let mut triangulation = SurfaceTriangulation::new((v1, 0), (v2, 0), (v3, 0));
        triangulation.locate(v4);
        triangulation.add_point(v4, 0);
        eprintln!("Locationg #5");
        let location = triangulation.locate(v5);
        eprintln!("edge_location: {:?}", location);
    }

    #[test]
    fn troublesome_insertion() {
        let vs: Vec<_> = vec![(0.0, 0.0), (2.5, -4.33), (5.0, 0.0)]
            .into_iter()
            .map(|(x, y)| Point::new(x, y))
            .collect();
        let mut triangulation = SurfaceTriangulation::new((vs[0], 0), (vs[1], 0), (vs[2], 0));
        for point in vs {
            triangulation.add_point(point, 0);
        }
        debug_assert_spaces(&triangulation);
        eprintln!("n_triangles: {}", triangulation.triangles().count());
        triangulation.add_point(Point::new(3.4375, -2.70625), 0);
        debug_assert_spaces(&triangulation);
        eprintln!("n_triangles: {}", triangulation.triangles().count());
    }

    #[test]
    fn troublesome_insertion_2() {
        let vs: Vec<_> = vec![(0.0, 0.0), (2.5, -4.33), (5.0, 0.0)]
            .into_iter()
            .map(|(x, y)| Point::new(x, y))
            .collect();
        let mut triangulation = SurfaceTriangulation::new((vs[0], 0), (vs[1], 0), (vs[2], 0));
        for point in vs {
            triangulation.add_point(point, 0);
        }
        debug_assert_spaces(&triangulation);
        eprintln!("n_triangles: {}", triangulation.triangles().count());
        triangulation.add_point(Point::new(1.640625, -0.135313), 0);
        debug_assert_spaces(&triangulation);
        eprintln!("n_triangles: {}", triangulation.triangles().count());
    }

    #[test]
    fn troublesome_insertion_3() {
        let vs: Vec<_> = vec![
            (0.0, 0.0),
            (2.5, -4.33),
            (5.0, 0.0),
            (1.25, -2.165),
            (1.875, -3.2475),
            (2.1875, -3.78875),
        ]
        .into_iter()
        .map(|(x, y)| Point::new(x, y))
        .collect();
        let mut triangulation = SurfaceTriangulation::new((vs[0], 0), (vs[1], 0), (vs[2], 0));
        for point in vs {
            triangulation.add_point(point, 0);
        }
        debug_assert_spaces(&triangulation);
        eprintln!("n_triangles: {}", triangulation.triangles().count());
    }

    #[test]
    fn troublesome_insertion_4() {
        let v1 = Point::new(0.0, 0.0);
        let v2 = Point::new(2.5, -4.33);
        let v3 = Point::new(5.0, 0.0);
        let vs: Vec<_> = vec![].into_iter().map(|(x, y)| Point::new(x, y)).collect();
        let mut triangulation = SurfaceTriangulation::new((v1, ()), (v2, ()), (v3, ()));
        triangulation.add_point((v1 + v2) * 0.5, ());
        triangulation.add_point((v2 + v3) * 0.5, ());
        triangulation.add_point((v3 + v1) * 0.5, ());
        debug_assert_spaces(&triangulation);
        eprintln!("n_triangles: {}", triangulation.triangles().count());
    }

    #[test]
    fn out_of_bounds() {
        // TODO: this test may need fixing
        let vs: Vec<_> = vec![(0.0, 0.0), (2.5, -4.33), (5.0, 0.0)]
            .into_iter()
            .map(|(x, y)| Point::new(x, y))
            .collect();
        let triangulation = SurfaceTriangulation::new((vs[0], 0), (vs[1], 0), (vs[2], 0));
        debug_assert_spaces(&triangulation);
        let location = triangulation.locate(Point::new(3.4375, -2.70625));
        eprintln!(
            "{:?}",
            location.map(|x| triangulation.get_segment(x.edge()))
        );
        debug_assert_spaces(&triangulation);
    }

    #[test]
    fn out_of_bounds_2() {
        let vs: Vec<_> = vec![(0.0, 0.0), (2.5, -4.33), (5.0, 0.0)]
            .into_iter()
            .map(|(x, y)| Point::new(x, y))
            .collect();
        let mut triangulation = SurfaceTriangulation::new((vs[0], 0), (vs[1], 0), (vs[2], 0));
        triangulation.add_point(Point::new(1.25, -2.165), 0);
        debug_assert_spaces(&triangulation);
        let p = Point::new(1.875, -3.2475);
        let location = triangulation.locate(p);
        assert!(location.is_some());
        triangulation.add_point(p, 0);
        debug_assert_spaces(&triangulation);
    }

    #[test]
    fn split_edges() {
        let vs: Vec<_> = vec![(0.0, 0.0), (5.0, 0.0), (5.0, 5.0)]
            .into_iter()
            .map(|(x, y)| Point::new(x, y))
            .collect();
        let mut triangulation = SurfaceTriangulation::<()>::new_with_default(vs[0], vs[1], vs[2]);
        debug_assert_spaces(&triangulation);
        let mut next = None;
        'the_loop: loop {
            if let Some((edge, midpoint)) = next {
                let (p1, p2) = triangulation.get_segment(triangulation.qeds.edge_a_ref(edge));
                eprintln!("adding {} to {}-{}", midpoint, p1, p2);
                triangulation.add_point_to_edge(edge, midpoint, ());
            }
            for edge in triangulation.base_targets() {
                let e = triangulation.qeds.edge_a_ref(edge);
                let (p1, p2) = triangulation.get_segment(e);
                let d = p1.distance(p2);
                if d > 1.0 && !triangulation.is_boundary(e) {
                    let midpoint = (p1 + p2) * 0.5;
                    next = Some((edge, midpoint));
                    continue 'the_loop;
                }
            }
            break;
        }
        // triangulation.add_point(Point::new(1.25, -2.165), 0);
        // let p = Point::new(1.875, -3.2475);
        // let location = triangulation.locate(p);
        // assert!(location.is_some());
        // triangulation.add_point(p, 0);
        // debug_assert_spaces(&triangulation);
    }

    #[ignore]
    #[test]
    fn add_to_boundary() {
        let p1 = Point::new(0.0, 0.0);
        let p2 = Point::new(2.5, -4.33);
        let p3 = Point::new(5.0, 0.0);
        let vs: Vec<_> = vec![p1, p2, p3];
        let p4 = (p1 + p2) * 0.5;
        let p5 = (p2 + p3) * 0.5;
        let p6 = (p3 + p1) * 0.5;

        let mut triangulation = SurfaceTriangulation::<()>::new_with_default(vs[0], vs[1], vs[2]);
        {
            let e1 = triangulation.get_matching_edge(p1, p2).unwrap();
            let e2 = triangulation.get_matching_edge(p2, p3).unwrap();
            let e3 = triangulation.get_matching_edge(p3, p1).unwrap();

            assert_eq!(e1.onext(), e3.sym());
            assert_eq!(e1.oprev(), e3.sym());
            assert_eq!(e1.rot(), e2.rot().onext());

            assert_eq!(e2.onext(), e1.sym());
            assert_eq!(e2.oprev(), e1.sym());

            assert_eq!(e3.onext(), e2.sym());
            assert_eq!(e3.oprev(), e2.sym());
        }
        let location = triangulation.locate(p4).unwrap();
        match location {
            Location::OnEdge(edge) => {
                eprintln!("checking is_outward_boundary");
                assert!(!triangulation.is_outward_boundary(edge));
                let segment = triangulation.get_segment(edge);
                assert_eq!(segment, (vs[0], vs[1]));
            }
            _ => panic!("wrong location"),
        }
        assert_eq!(1, triangulation.triangles().count());
        triangulation.add_point_with_default(p4);
        triangulation.add_point_with_default(p5);
        triangulation.add_point_with_default(p6);
        debug_assert_spaces(&triangulation);
        assert_eq!(4, triangulation.triangles().count());
        let p7 = (p2 + p5) * 0.5;
        triangulation.add_point_with_default(p7);
        debug_assert_spaces(&triangulation);
        {
            let e1 = triangulation.get_matching_edge(p1, p4).unwrap();
            let e2 = triangulation.get_matching_edge(p4, p2).unwrap();
            let e3 = triangulation.get_matching_edge(p2, p7).unwrap();
            let e4 = triangulation.get_matching_edge(p7, p5).unwrap();
            let e5 = triangulation.get_matching_edge(p5, p3).unwrap();
            let e6 = triangulation.get_matching_edge(p3, p6).unwrap();
            let e7 = triangulation.get_matching_edge(p6, p1).unwrap();
            let e8 = triangulation.get_matching_edge(p6, p4).unwrap();
            let e9 = triangulation.get_matching_edge(p4, p5).unwrap();
            let e10 = triangulation.get_matching_edge(p5, p6).unwrap();
            let e11 = triangulation.get_matching_edge(p4, p7).unwrap();

            assert_eq!(e1.onext(), e7.sym());
            assert_eq!(e1.oprev(), e7.sym());
            assert_eq!(e1.l_next(), e8.sym());
            assert_eq!(e1.d_prev(), e8);

            assert_eq!(e1.sym().l_next(), e7.sym());
            assert_eq!(e1.sym().l_next().l_next(), e6.sym());
            assert_eq!(e1.sym().l_next().l_next().l_next(), e5.sym());
            assert_eq!(e1.sym().l_next().l_next().l_next().l_next(), e4.sym());
            assert_eq!(
                e1.sym().l_next().l_next().l_next().l_next().l_next(),
                e3.sym()
            );
            assert_eq!(
                e1.sym()
                    .l_next()
                    .l_next()
                    .l_next()
                    .l_next()
                    .l_next()
                    .l_next(),
                e2.sym()
            );
            assert_eq!(
                e1.sym()
                    .l_next()
                    .l_next()
                    .l_next()
                    .l_next()
                    .l_next()
                    .l_next()
                    .l_next(),
                e1.sym()
            );

            assert_eq!(e2.onext(), e11);
            assert_eq!(e2.oprev(), e1.sym());
            assert_eq!(e2.l_next(), e3);
            assert_eq!(e2.l_next().l_next(), e11.sym());

            assert_eq!(e7.l_next(), e1);
            assert_eq!(e7.l_next().l_next(), e8.sym());
        }
        let p8 = (p3 + p5) * 0.5;
        eprintln!("adding point: {}", p8);
        let location = triangulation.locate_force(p8).unwrap();
        {
            let e = triangulation.get_matching_edge(p3, p5).unwrap();
            eprintln!("e: {:?}", e);
            eprintln!("e.sym(): {:?}", e.sym(),);
            eprintln!("location.edge(): {:?}", location.edge());
            eprintln!(
                "specified: {:?} found: {:?}",
                p8,
                triangulation.get_segment(e)
            );
            let e = if triangulation.curves_left(e) {
                e
            } else {
                e.sym()
            };
            assert_eq!(e, location.edge(), "p8 located on wrong edge");
        }
        triangulation.add_point_with_default(p8);
        debug_assert_spaces(&triangulation);
    }

    #[test]
    fn edge_location_outside() {
        let v1 = Point {
            x: -12.499999999998696,
            y: 0.000001527352651415037,
        };
        let v2 = Point {
            x: -16.249999999998355,
            y: 0.000002148013987930656,
        };
        let v3 = Point {
            x: -14.37499999999848,
            y: -0.9374981413664982,
        };
        let v4 = Point::new(-6.875, 0.000001);
        let triangulation = SurfaceTriangulation::new((v1, 0), (v2, 0), (v3, 0));
        assert!(triangulation.locate(v4).is_none());
        // triangulation.add_point(v4, 0, true);
    }

    #[test]
    fn t10() {
        // use qeds::point::Point;
        let points = vec![
            Point {
                x: -17.870000839231405,
                y: 2.7686743613015225e-6,
            },
            Point {
                x: -4.218110561368454,
                y: -15.000001861413416,
            },
            Point {
                x: -4.218110561370102,
                y: -2.1966163308234274e-6,
            },
            Point {
                x: -11.04405570029993,
                y: -7.499999546369527,
            },
            Point {
                x: -7.631083130834192,
                y: -11.250000703891471,
            },
        ];
        let mut triangulation =
            SurfaceTriangulation::<()>::new_with_default(points[0], points[1], points[2]);
        debug_assert_spaces(&triangulation);
        {
            for edge in triangulation.base_targets() {
                let (p1, p2) = triangulation.get_segment(triangulation.qeds.edge_a_ref(edge));
                eprintln!("edge_distance: {}", p1.distance(p2));
            }
            // panic!("end");
        }
        {
            let edge = triangulation.get_matching_edge_indices(0, 1).unwrap();
            let edge_target = edge.target();
            let (p1, p2) = triangulation.get_segment(edge);
            assert!(p1.distance(p2) > 11.0);
            let midpoint = p1.midpoint(p2);
            triangulation.add_point_to_edge(edge_target, midpoint, ());
            debug_assert_spaces(&triangulation);
        }
        {
            let edge = triangulation.get_matching_edge_indices(1, 2).unwrap();
            let edge_target = edge.target();
            let (p1, p2) = triangulation.get_segment(edge);
            assert!(p1.distance(p2) > 11.0);
            let midpoint = p1.midpoint(p2);
            triangulation.add_point_to_edge(edge_target, midpoint, ());
            debug_assert_spaces(&triangulation);
        }
        {
            let edge = triangulation.get_matching_edge_indices(2, 0).unwrap();
            let edge_target = edge.target();
            let (p1, p2) = triangulation.get_segment(edge);
            assert!(p1.distance(p2) > 11.0);
            let midpoint = p1.midpoint(p2);
            triangulation.add_point_to_edge(edge_target, midpoint, ());
            debug_assert_spaces(&triangulation);
        }
        // triangulation.add_point(points[4], ());
        debug_assert_spaces(&triangulation);
        // assert_eq!(4, triangulation.triangles().count());
    }
}

/// Assert that every node as a component.
pub fn debug_assert_spaces<T: Clone + Serialize>(triangulation: &SurfaceTriangulation<T>) {
    #[cfg(debug_assertions)]
    {
        triangulation.debug_dump(None);
        {
            // Debug spokes
            for edge in triangulation.qeds.base_edges() {
                let vertex_index = edge.edge().point;
                let mut current = edge;
                let mut i = 0;
                loop {
                    assert_eq!(vertex_index, current.edge().point);
                    current = current.onext();
                    if current == edge {
                        break;
                    }
                    i += 1;
                    if i > 2000 {
                        panic!("hit iteration limit");
                    }
                }
            }
        }
        {
            for (i, quad) in triangulation.qeds.quads.iter() {
                assert_ne!(quad.edges_a[0].point,quad.edges_a[1].point);
            }
        }
        // Get all B edges.
        let mut b_targets = Vec::new();
        for (i, _quad) in triangulation.qeds.quads.iter() {
            b_targets.push(EdgeTarget::new(i, 1, 0));
            b_targets.push(EdgeTarget::new(i, 3, 0));
        }
        for target in triangulation.base_targets() {
            let edge_ref = triangulation.qeds.edge_a_ref(target);
            let segment = triangulation.get_segment(edge_ref);
            // eprintln!("edge: {}-{}", segment.0, segment.1);
            assert_ne!(edge_ref, edge_ref.d_prev(), "d_prev is self");
            assert_ne!(edge_ref, edge_ref.onext(), "onext is self");
        }
        // eprintln!("Checking triangles");
        for triangle in triangulation.triangles() {
            let ccw =
                crate::triangulation::is_ccw(triangle.0.point, triangle.1.point, triangle.2.point);
            if !ccw {
                eprintln!(
                    "{}-{}-{}",
                    triangle.0.point, triangle.1.point, triangle.2.point
                );
            }
            assert!(ccw);
        }
    }
}
