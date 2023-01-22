use crate::math::{Real, Vector};
use crate::query::{Ray, RayCast, RayIntersection};
use crate::shape::{Cuboid, VoxelType, Voxels};

impl RayCast for Voxels {
    #[inline]
    fn cast_local_ray(&self, ray: &Ray, max_toi: Real, solid: bool) -> Option<Real> {
        // TODO: optimize this very naive implementation.
        let base_cuboid = Cuboid::new(Vector::repeat(self.voxel_size() / 2.0));
        let mut result: Option<Real> = None;
        for (center, data) in self.centers() {
            if data.voxel_type() != VoxelType::Empty {
                let shifted_ray = ray.translate_by(-center.coords);
                if let Some(candidate) = base_cuboid.cast_local_ray(&shifted_ray, max_toi, solid) {
                    if let Some(result) = &mut result {
                        *result = result.min(candidate);
                    } else {
                        result = Some(candidate);
                    }
                }
            }
        }

        result
    }

    #[inline]
    fn cast_local_ray_and_get_normal(
        &self,
        ray: &Ray,
        max_toi: Real,
        solid: bool,
    ) -> Option<RayIntersection> {
        // TODO: optimize this very naive implementation.
        let base_cuboid = Cuboid::new(Vector::repeat(self.voxel_size() / 2.0));
        let mut result: Option<RayIntersection> = None;
        for (center, data) in self.centers() {
            if data.voxel_type() != VoxelType::Empty {
                let shifted_ray = ray.translate_by(-center.coords);
                if let Some(candidate) =
                    base_cuboid.cast_local_ray_and_get_normal(&shifted_ray, max_toi, solid)
                {
                    if let Some(result) = &mut result {
                        if result.toi > candidate.toi {
                            *result = candidate;
                        }
                    } else {
                        result = Some(candidate);
                    }
                }
            }
        }

        result
    }
}
