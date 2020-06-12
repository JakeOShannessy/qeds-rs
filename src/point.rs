use std::{
    fmt::{Display, Formatter},
    ops::{Add, Mul, Sub},
};

// TODO: create new wrapper type around point (and possibly the floats within)
// to limit to safe operations. This wrapper will enforce the minimum signficant
// unit value (i.e. the smallest epsilon) and force us to check for overflows
// (unless we opt-out).
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

#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct SafeFloat(f64);

impl SafeFloat {
    pub const MAX: SafeFloat = SafeFloat( 0xffffffff_i64 as f64);
    pub const MIN: SafeFloat = SafeFloat(1.0/(0xfffff_i64 as f64) as f64);
}


impl Display for SafeFloat {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.0)
    }
}

impl Add for SafeFloat {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        SafeFloat(self.0+other.0)
    }
}

impl Mul for SafeFloat {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        SafeFloat(self.0*other.0)
    }
}

impl Sub for SafeFloat {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        SafeFloat(self.0-other.0)
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn safe_distance() {
        println!("SafeFloat::MAX: {}", SafeFloat::MAX);
        println!("SafeFloat::MIN: {}", SafeFloat::MIN);
        assert_ne!(SafeFloat::MAX, SafeFloat::MAX - SafeFloat::MIN);
    }
}
