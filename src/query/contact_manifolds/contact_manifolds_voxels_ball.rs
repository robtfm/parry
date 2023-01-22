use crate::bounding_volume::BoundingVolume;
use crate::math::{Isometry, Point, Real, Vector};
use crate::query::{ContactManifold, TrackedContact};
use crate::shape::{Ball, OctantPattern, PackedFeatureId, Shape, VoxelData, VoxelType, Voxels};

/// Computes the contact manifold between a convex shape and a ball, both represented as a `Shape` trait-object.
pub fn contact_manifolds_voxels_ball_shapes<ManifoldData, ContactData>(
    pos12: &Isometry<Real>,
    shape1: &dyn Shape,
    shape2: &dyn Shape,
    prediction: Real,
    manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
) where
    ManifoldData: Default,
    ContactData: Default + Copy,
{
    if let (Some(voxels1), Some(ball2)) = (shape1.as_voxels(), shape2.as_ball()) {
        contact_manifolds_voxels_ball(pos12, voxels1, ball2, prediction, manifolds, false);
    } else if let (Some(ball1), Some(voxels2)) = (shape1.as_ball(), shape2.as_voxels()) {
        contact_manifolds_voxels_ball(
            &pos12.inverse(),
            voxels2,
            ball1,
            prediction,
            manifolds,
            true,
        );
    }
}

/// Computes the contact manifold between a convex shape and a ball.
pub fn contact_manifolds_voxels_ball<'a, ManifoldData, ContactData>(
    pos12: &Isometry<Real>,
    voxels1: &'a Voxels,
    ball2: &'a Ball,
    prediction: Real,
    manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
    flipped: bool,
) where
    ManifoldData: Default,
    ContactData: Default + Copy,
{
    // TODO: don’t generate one manifold per voxel.
    manifolds.clear();

    let radius1 = voxels1.scale / 2.0;
    let radius2 = ball2.radius;
    let center2 = Point::origin(); // The ball’s center.

    // FIXME: optimize this.
    let aabb1 = voxels1.local_aabb().loosened(prediction / 2.0);
    let aabb2 = ball2.aabb(pos12).loosened(prediction / 2.0);
    if let Some(aabb_intersection) = aabb1.intersection(&aabb2) {
        for (center1, data1) in voxels1.voxels_intersecting_local_aabb(&aabb_intersection) {
            let voxel1 = data1.voxel_type();
            match voxel1 {
                VoxelType::Vertex | VoxelType::Edge | VoxelType::Face => { /* Ok */ }
                _ => continue,
            }

            detect_hit_voxel_ball(
                *pos12, center1, radius1, data1, center2, radius2, prediction, flipped, manifolds,
            );
        }
    }
}

pub(crate) fn detect_hit_voxel_ball<ManifoldData, ContactData>(
    pos12: Isometry<Real>,
    center1: Point<Real>,
    radius1: Real,
    data1: VoxelData,
    center2: Point<Real>,
    radius2: Real,
    prediction: Real,
    flipped: bool,
    manifolds: &mut Vec<ContactManifold<ManifoldData, ContactData>>,
) where
    ManifoldData: Default,
    ContactData: Default + Copy,
{
    let center2_1 = pos12 * center2;

    if let Some((local_n1, center_dist)) =
        project_point_in_pseudo_ball(center1, radius1, data1.octant_mask(), center2_1)
    {
        let dist = center_dist - radius1 - radius2;

        if dist <= prediction {
            let local_n2 = pos12.inverse_transform_vector(&-local_n1);
            let local_p1 = center2_1 - local_n1 * (center_dist - radius1);
            let local_p2 = center2 + local_n2 * radius2;
            let contact_point = TrackedContact::<ContactData>::flipped(
                local_p1,
                local_p2,
                PackedFeatureId::UNKNOWN,
                PackedFeatureId::UNKNOWN,
                dist,
                flipped,
            );

            let mut manifold = ContactManifold::<ManifoldData, ContactData>::new();
            manifold.points.push(contact_point);

            if flipped {
                manifold.local_n1 = local_n2;
                manifold.local_n2 = local_n1;
            } else {
                manifold.local_n1 = local_n1;
                manifold.local_n2 = local_n2;
            }

            manifolds.push(manifold);
        }
    }
}

pub fn project_point_in_pseudo_ball(
    ball_center: Point<Real>,
    ball_radius: Real,
    ball_mask: u32,
    point: Point<Real>,
) -> Option<(Vector<Real>, Real)> {
    let dpos = point - ball_center;

    // This maps our easy-to map octant_key, to the vertex key as
    // used by the Aabb::FACES_VERTEX_IDS and Aabb::EDGES_VERTEX_IDS
    // arrays. Could we find a way to avoid this map by computing the
    // correct key right away?
    const AABB_OCTANT_KEYS: [u32; 8] = [0, 1, 3, 2, 4, 5, 7, 6];
    #[cfg(feature = "dim2")]
    let octant_key = (((dpos.x >= 0.0) as usize) << 0) | (((dpos.y >= 0.0) as usize) << 1);
    #[cfg(feature = "dim3")]
    let octant_key = (((dpos.x >= 0.0) as usize) << 0)
        | (((dpos.y >= 0.0) as usize) << 1)
        | (((dpos.z >= 0.0) as usize) << 2);
    let aabb_octant_key = AABB_OCTANT_KEYS[octant_key];
    let dpos_signs = dpos.map(|x| x.signum());
    let unit_dpos = dpos.component_mul(&dpos_signs) / ball_radius; // Project the point in "local unit octant space".
    let pattern = (ball_mask >> (aabb_octant_key * 3)) & 0b0111;

    let unit_result = match pattern {
        OctantPattern::INTERIOR => None,
        OctantPattern::VERTEX => {
            let mut normal = unit_dpos;
            let dist = normal.try_normalize_mut(1.0e-8)?;
            Some((normal, dist))
        }
        OctantPattern::EDGE_X | OctantPattern::EDGE_Y | OctantPattern::EDGE_Z => {
            let i = pattern as usize - OctantPattern::EDGE_X as usize;

            if unit_dpos[i] > 1.0 || unit_dpos[i] < 0.0 {
                None
            } else {
                let mut normal = unit_dpos;
                normal[i] = 0.0;
                let dist = normal.try_normalize_mut(1.0e-8)?;
                Some((normal, dist))
            }
        }
        OctantPattern::FACE_X | OctantPattern::FACE_Y | OctantPattern::FACE_Z => {
            let i1 = pattern as usize - OctantPattern::FACE_X as usize;
            let i2 = (i1 + 1) % 3;
            let i3 = (i1 + 2) % 3;

            if unit_dpos[i2] > 1.0
                || unit_dpos[i2] < 0.0
                || unit_dpos[i3] > 1.0
                || unit_dpos[i3] < 0.0
            {
                None
            } else {
                let dist = unit_dpos[i1];
                Some((Vector::ith(i1, 1.0), dist))
            }
        }
        _ => unreachable!(),
    };

    unit_result.map(|(n, d)| (n.component_mul(&dpos_signs), d * ball_radius))
}
