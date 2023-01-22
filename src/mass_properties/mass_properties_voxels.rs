use crate::mass_properties::MassProperties;
use crate::math::{AngularInertia, Point, PrincipalAngularInertia, Real, Vector};
use crate::shape::{Cuboid, VoxelType, Voxels};

impl MassProperties {
    /// Computes the mass properties of a set of voxels.
    pub fn from_voxels(density: Real, voxels: &Voxels) -> Self {
        let mut com = Point::origin();
        let mut num_not_empty = 0;
        let mut angular_inertia = na::zero();
        let block_ref_mprops =
            MassProperties::from_cuboid(density, Vector::repeat(voxels.scale / 2.0));

        for (center, data) in voxels.centers() {
            if !data.is_empty() {
                com += center.coords;
                num_not_empty += 1;
            }
        }

        com.coords /= num_not_empty as Real;

        for (center, data) in voxels.centers() {
            if !data.is_empty() {
                angular_inertia += block_ref_mprops.construct_shifted_inertia_matrix(center - com);
            }
        }

        let mass = block_ref_mprops.mass() * num_not_empty as Real;

        #[cfg(feature = "dim2")]
        return Self::new(com, mass, angular_inertia);
        #[cfg(feature = "dim3")]
        return Self::with_inertia_matrix(com, mass, angular_inertia);
    }
}
