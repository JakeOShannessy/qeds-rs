#![allow(dead_code)]
#[cfg(test)]
extern crate quickcheck;
#[cfg(test)]
#[macro_use(quickcheck)]
extern crate quickcheck_macros;

pub mod constrained_triangulation;
pub mod point;
pub mod qeds;
pub mod triangulation;
