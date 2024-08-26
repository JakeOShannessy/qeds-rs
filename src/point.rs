#[cfg(feature = "serialize")]
use serde::{Deserialize, Serialize};
use std::{
    fmt::{Display, Formatter},
    ops::{Add, Mul, Sub},
};

// TODO: create new wrapper type around point (and possibly the floats within)
// to limit to safe operations. This wrapper will enforce the minimum signficant
// unit value (i.e. the smallest epsilon) and force us to check for overflows
// (unless we opt-out).
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
#[cfg_attr(feature = "serialize", derive(Serialize, Deserialize))]
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

    pub fn magnitude(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2)).sqrt()
    }

    pub fn snap(&mut self) {
        self.x *= 1.0e6;
        self.x = self.x.round();
        self.x *= 1.0e-6;

        self.y *= 1.0e6;
        self.y = self.y.round();
        self.y *= 1.0e-6;
    }
}

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3 {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    pub fn midpoint(self, other: Self) -> Self {
        Self {
            x: (self.x + other.x) / 2.0,
            y: (self.y + other.y) / 2.0,
            z: (self.z + other.z) / 2.0,
        }
    }

    pub fn is_finite(&self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    pub fn distance(&self, other: Self) -> f64 {
        ((other.x - self.x).powi(2) + (other.y - self.y).powi(2) + (other.z - self.z).powi(2))
            .sqrt()
    }

    pub fn snap(&mut self) {
        self.x *= 1.0e6;
        self.x = self.x.round();
        self.x *= 1.0e-6;

        self.y *= 1.0e6;
        self.y = self.y.round();
        self.y *= 1.0e-6;

        self.y *= 1.0e6;
        self.y = self.z.round();
        self.y *= 1.0e-6;
    }
}

#[cfg(test)]
impl quickcheck::Arbitrary for Point {
    fn arbitrary(g: &mut quickcheck::Gen) -> Point {
        Point {
            x: SafeFloat::arbitrary(g).0,
            y: SafeFloat::arbitrary(g).0,
        }
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

impl From<Point> for robust::Coord<f64> {
    fn from(s: Point) -> robust::Coord<f64> {
        robust::Coord { x: s.x, y: s.y }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct SafeFloat(pub f64);

impl SafeFloat {
    // pub const MAX: SafeFloat = SafeFloat( 0xffffffff_i64 as f64);
    // pub const MIN: SafeFloat = SafeFloat( -0xffffffff_i64 as f64);
    // pub const EPSILON: SafeFloat = SafeFloat(1.0/(0xfffff_i64 as f64) as f64);

    // For tests
    pub const MAX: SafeFloat = SafeFloat(4000.0);
    pub const MIN: SafeFloat = SafeFloat(-4000.0);
    pub const EPSILON: SafeFloat = SafeFloat(1.0 / (0xf_i64 as f64));
    pub fn new(f: f64) -> Option<Self> {
        let s = SafeFloat(f);
        if s <= SafeFloat::MAX && s >= SafeFloat::MIN {
            Some(s)
        } else {
            None
        }
    }
}

#[cfg(test)]
impl quickcheck::Arbitrary for SafeFloat {
    fn arbitrary(g: &mut quickcheck::Gen) -> Self {
        loop {
            let frac: f64 = f64::arbitrary(g);
            let f = SafeFloat::new(frac);
            if let Some(f) = f {
                break f;
            }
        }
    }
}

impl Display for SafeFloat {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.0)
    }
}

impl Add for SafeFloat {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        SafeFloat(self.0 + other.0)
    }
}

impl Mul for SafeFloat {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        SafeFloat(self.0 * other.0)
    }
}

impl Sub for SafeFloat {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        SafeFloat(self.0 - other.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn safe_distance() {
        println!("SafeFloat::MAX: {}", SafeFloat::MAX);
        println!("SafeFloat::MIN: {}", SafeFloat::MIN);
        println!("SafeFloat::EPSILON: {}", SafeFloat::EPSILON);
        assert_ne!(SafeFloat::MAX, SafeFloat::MAX - SafeFloat::EPSILON);
    }
}
