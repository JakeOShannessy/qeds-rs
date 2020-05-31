use qeds::*;
fn main() {
    println!(
        "Edge Size: {:?}(0x{:x?})",
        std::mem::size_of::<Edge>(),
        std::mem::size_of::<Edge>()
    );
    println!(
        "Quad Size: {:?}(0x{:x?})",
        std::mem::size_of::<Quad>(),
        std::mem::size_of::<Quad>()
    );
    println!(
        "Quad Alignment: {:?}(0x{:x?})",
        std::mem::align_of::<Quad>(),
        std::mem::align_of::<Quad>()
    );
    // let quad_layout = std::alloc::Layout::from_size_align(std::mem::size_of::<Quad>(), 64).unwrap();
    // let quad_layout_standard = std::alloc::Layout::from_size_align(std::mem::size_of::<Quad>(), 64).unwrap();
    // println!("Quad Layout: {:?}", quad_layout);
    // println!("Quad Layout Std: {:?}",quad_layout_standard);

    // Step 1. Create a Qeds data structure.
    let mut qeds = Qeds::new();
    // Step 2. Add some edges to it.
    let q1 = qeds.make_edge(Point::new(0.0,0.0));
    let q2 = qeds.make_edge(Point::new(5.0,0.0));
    println!("QEDS1: {:?}", qeds);
    for quad in qeds.quads.iter() {
        let quad: &Quad = unsafe { &**quad };
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
    println!("--------");
    for quad in qeds.quads.iter() {
        let quad: &Quad = unsafe { &**quad };
        let e0 = &quad.edges[0];
        println!("Edge: {:?} -> {:?}", e0, e0.sym().onext());
    }
    unsafe {
        let e1 = &mut (*q1).edges[2];
        let e2 = &mut (*q2).edges[0];
        qeds.splice(e1, e2)
    }
    println!("--------");
    for quad in qeds.quads.iter() {
        let quad: &Quad = unsafe { &**quad };
        let e0 = &quad.edges[0];
        println!("Edge: {:?} -> {:?}", e0, e0.sym().onext());
    }
    println!("--------");
    println!("QEDS2: {:?}", qeds);
}
