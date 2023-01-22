use crate::bounding_volume::BoundingVolume;
use crate::math::{Isometry, Real};
use crate::query::{ContactManifold, TrackedContact};
use crate::shape::{FeatureId, PackedFeatureId, Shape, VoxelType, Voxels};
use na::Unit;

/// Computes the contact manifold between a convex shape and a ball, both represented as a `Shape` trait-object.
pub fn contact_manifolds_voxels_shape_shapes<ManifoldData, ContactData>(
    pos12: &Isometry<Real>,
    shape1: &dyn Shape,
    shape2: &dyn Shape,
    prediction: Real,
    manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
) where
    ManifoldData: Default,
    ContactData: Default + Copy,
{
    if let Some(voxels1) = shape1.as_voxels() {
        contact_manifolds_voxels_shape(
            &pos12.inverse(),
            shape2,
            voxels1,
            prediction,
            manifolds,
            true,
        );
    } else if let Some(voxels2) = shape2.as_voxels() {
        contact_manifolds_voxels_shape(pos12, shape1, voxels2, prediction, manifolds, false);
    }
}

/// Computes the contact manifold between a convex shape and a ball.
pub fn contact_manifolds_voxels_shape<'a, ManifoldData, ContactData, S1>(
    pos12: &Isometry<Real>,
    shape1: &'a S1,
    voxels2: &'a Voxels,
    prediction: Real,
    manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
    flipped: bool,
) where
    S1: ?Sized + Shape,
    ManifoldData: Default,
    ContactData: Default + Copy,
{
    // TODO: don’t generate one manifold per voxel.
    manifolds.clear();

    let radius = voxels2.scale / 2.0;

    // FIXME: optimize this.
    let aabb1 = shape1.compute_local_aabb().loosened(prediction);
    let aabb2 = voxels2.local_aabb();

    // if let Some((_, intersection_aabb2)) = aabb1.aligned_intersections(pos12, &aabb2) {
    for (center2, data) in voxels2.centers() {
        // .voxels_intersecting_local_aabb(&intersection_aabb2) {
        let voxel = data.voxel_type();
        if voxel == VoxelType::Empty || voxel == VoxelType::Interior {
            continue;
        }

        let center2_1 = pos12 * center2;
        let (proj, feat) = shape1.project_local_point_and_get_feature(&center2_1);

        match (voxel, feat) {
            (VoxelType::Vertex, _) | (_, FeatureId::Vertex(_)) => { /* OK */ }
            #[cfg(feature = "dim3")]
            (VoxelType::Edge, FeatureId::Edge(_)) => { /* OK */ }
            _ => continue,
        }

        let dpos = center2_1 - proj.point;
        let projection = Unit::try_new_and_get(dpos, 0.0);

        if let Some((mut local_n1, mut dist)) = projection {
            if proj.is_inside {
                local_n1 = -local_n1;
                dist = -dist;
            }

            if dist <= radius + prediction {
                let local_n2 = pos12.inverse_transform_vector(&-*local_n1);
                let local_p2 = center2 + local_n2 * radius;
                let contact_point = TrackedContact::<ContactData>::flipped(
                    proj.point,
                    local_p2,
                    PackedFeatureId::UNKNOWN,
                    PackedFeatureId::UNKNOWN,
                    dist - radius,
                    flipped,
                );

                let mut manifold = ContactManifold::<ManifoldData, ContactData>::new();
                manifold.points.push(contact_point);

                if flipped {
                    manifold.local_n1 = local_n2;
                    manifold.local_n2 = *local_n1;
                } else {
                    manifold.local_n1 = *local_n1;
                    manifold.local_n2 = local_n2;
                }

                manifolds.push(manifold);
            }
        } else {
            println!("Failed");
        }
    }
    // }
}
