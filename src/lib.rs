// #![warn(rust_2018_idioms)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::too_many_arguments)]
#![allow(dead_code)]
#[cfg(test)]
extern crate quickcheck;
#[cfg(test)]
#[macro_use(quickcheck)]
extern crate quickcheck_macros;

pub mod constrained_triangulation;
pub mod point;
pub mod qeds;
pub mod surface_triangulation;
pub mod triangulation;
