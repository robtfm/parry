use crate::bounding_volume::{Aabb, BoundingVolume};
use crate::math::{Isometry, Real, Vector};
use crate::query::ContactManifold;
use crate::shape::{Shape, VoxelType, Voxels};

/// Computes the contact manifold between a convex shape and a ball, both represented as a `Shape` trait-object.
pub fn contact_manifolds_voxels_voxels_shapes<ManifoldData, ContactData>(
    pos12: &Isometry<Real>,
    shape1: &dyn Shape,
    shape2: &dyn Shape,
    prediction: Real,
    manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
) where
    ManifoldData: Default,
    ContactData: Default + Copy,
{
    if let (Some(voxels1), Some(voxels2)) = (shape1.as_voxels(), shape2.as_voxels()) {
        contact_manifolds_voxels_voxels(pos12, voxels1, voxels2, prediction, manifolds);
    }
}

/// Computes the contact manifold between a convex shape and a ball.
pub fn contact_manifolds_voxels_voxels<'a, ManifoldData, ContactData>(
    pos12: &Isometry<Real>,
    voxels1: &'a Voxels,
    voxels2: &'a Voxels,
    prediction: Real,
    manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
) where
    ManifoldData: Default,
    ContactData: Default + Copy,
{
    // TODO: don’t generate one manifold per voxel.
    manifolds.clear();

    let radius1 = voxels1.scale / 2.0;
    let radius2 = voxels2.scale / 2.0;

    // FIXME: optimize this.
    let aabb1 = voxels1.local_aabb().loosened(prediction / 2.0);
    let aabb2 = voxels2.local_aabb().loosened(prediction / 2.0);

    if let Some((intersection_aabb1, intersection_aabb2)) =
        aabb1.aligned_intersections(pos12, &aabb2)
    {
        let pos21 = pos12.inverse();

        for (center1, data1) in voxels1.voxels_intersecting_local_aabb(&intersection_aabb1) {
            let voxel1 = data1.voxel_type();
            match voxel1 {
                VoxelType::Vertex | VoxelType::Edge => { /* Ok */ }
                _ => continue,
            }

            let centered_aabb1_2 =
                Aabb::from_half_extents(pos21 * center1, Vector::repeat(radius1 + prediction));

            for (center2, data2) in voxels2.voxels_intersecting_local_aabb(&centered_aabb1_2) {
                let voxel2 = data2.voxel_type();
                match (voxel1, voxel2) {
                    (VoxelType::Vertex, VoxelType::Vertex)
                    | (VoxelType::Vertex, VoxelType::Edge)
                    | (VoxelType::Vertex, VoxelType::Face) => {
                        super::detect_hit_voxel_ball(
                            pos21, center2, radius2, data2, center1, radius1, prediction, true,
                            manifolds,
                        );
                    }
                    (VoxelType::Edge, VoxelType::Vertex)
                    | (VoxelType::Face, VoxelType::Vertex)
                    | (VoxelType::Edge, VoxelType::Edge) => {
                        super::detect_hit_voxel_ball(
                            *pos12, center1, radius1, data1, center2, radius2, prediction, false,
                            manifolds,
                        );
                    }
                    _ => continue, /* Ignore */
                }
            }
        }

        for (center2, data2) in voxels2.voxels_intersecting_local_aabb(&intersection_aabb2) {
            if data2.voxel_type() != VoxelType::Vertex {
                continue;
            }

            let centered_aabb2_1 =
                Aabb::from_half_extents(pos12 * center2, Vector::repeat(radius2 + prediction));

            for (center1, data1) in voxels1.voxels_intersecting_local_aabb(&centered_aabb2_1) {
                if data1.voxel_type() != VoxelType::Face {
                    continue;
                }

                super::detect_hit_voxel_ball(
                    *pos12, center1, radius1, data1, center2, radius2, prediction, false, manifolds,
                );
            }
        }
    }
}
