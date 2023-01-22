use crate::bounding_volume::Aabb;
use crate::math::{Point, Real, Vector, DIM};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum VoxelType {
    Empty,
    Vertex,
    Edge,
    Face,
    Interior,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct OctantPattern;

impl OctantPattern {
    pub const INTERIOR: u32 = 0;
    pub const VERTEX: u32 = 1;
    pub const EDGE_X: u32 = 2;
    pub const EDGE_Y: u32 = 3;
    pub const EDGE_Z: u32 = 4;
    pub const FACE_X: u32 = 5;
    pub const FACE_Y: u32 = 6;
    pub const FACE_Z: u32 = 7;
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
pub struct VoxelData(u8);

impl VoxelData {
    pub const fn empty() -> Self {
        Self(EMPTY_FACE_MASK)
    }

    pub const fn interior() -> Self {
        Self(INTERIOR_FACE_MASK)
    }

    pub const fn is_empty(self) -> bool {
        self.0 == EMPTY_FACE_MASK
    }

    pub const fn free_faces(self) -> u8 {
        if self.0 == INTERIOR_FACE_MASK || self.0 == EMPTY_FACE_MASK {
            0
        } else {
            (!self.0) & 0b0011_1111
        }
    }
    pub const fn voxel_type(self) -> VoxelType {
        FACES_TO_VOXEL_TYPES[self.0 as usize]
    }

    // Bitmask indicating what vertices, edges, or faces of the voxel are "free".
    pub const fn feature_mask(self) -> u16 {
        FACES_TO_FEATURE_MASKS[self.0 as usize]
    }
    pub const fn octant_mask(self) -> u32 {
        FACES_TO_OCTANT_MASKS[self.0 as usize]
    }
}

#[derive(Clone)]
#[cfg_attr(feature = "serde-serialize", derive(Serialize, Deserialize))]
pub struct Voxels {
    pub(crate) dimensions: [u32; DIM],
    pub(crate) data: Vec<VoxelData>, // TODO: use a Vob?
    pub scale: Real,
    pub origin: Point<Real>,
}

impl Voxels {
    pub fn from_points(points: &[Point<Real>], scale: Real) -> Self {
        if cfg!(feature = "dim2") {
            panic!("2D voxels not supported yet: need to generate the const tables.");
        }
        let mut aabb = Aabb::from_points(points);
        aabb.mins -= Vector::repeat(scale / 2.0);
        aabb.maxs += Vector::repeat(scale / 2.0); // NOTE: is this + necessary?

        // Discretize the points on the grid.
        let dimensions = ((aabb.maxs - aabb.mins) / scale).map(|x| x.ceil() as u32);
        let len = dimensions.product();
        let mut result = Self {
            dimensions: dimensions.into(),
            data: vec![VoxelData::empty(); len as usize],
            scale,
            origin: aabb.mins,
        };

        for pt in points {
            let coords = (pt - aabb.mins).map(|x| (x / scale).floor() as u32);
            let index = result.linear_index(coords.into());
            result.data[index as usize] = VoxelData::interior();
        }

        result.recompute_voxels_data();
        result
    }

    pub fn extents(&self) -> Vector<Real> {
        Vector::from(self.dimensions).cast::<Real>() * self.scale
    }

    pub fn voxel_size(&self) -> Real {
        self.scale
    }

    fn recompute_voxels_data(&mut self) {
        for i in 0..self.data.len() {
            let key = self.to_key(i as u32);

            self.data[i] = self.compute_voxel_data(key);
        }
    }

    pub fn scaled(mut self, scale: &Vector<Real>) -> Option<Self> {
        #[cfg(feature = "dim2")]
        let scale_is_uniform = scale.x == scale.y;
        #[cfg(feature = "dim3")]
        let scale_is_uniform = scale.x == scale.y && scale.y == scale.z;
        if !scale_is_uniform {
            // TODO: what about non-uniform scale?
            None
        } else {
            self.scale *= scale.x;
            self.origin *= scale.x;
            Some(self)
        }
    }

    fn quantify_point(&self, pt: Point<Real>) -> Vector<u32> {
        ((pt - self.origin) / self.scale).map(|x| x.floor().max(0.0) as u32)
    }

    pub fn voxels_intersecting_local_aabb(
        &self,
        aabb: &Aabb,
    ) -> impl Iterator<Item = (Point<Real>, VoxelData)> + '_ {
        let dims = Vector::from(self.dimensions);
        let mins = ((aabb.mins - self.origin) / self.scale)
            .map(|x| x.floor().max(0.0) as u32)
            .inf(&dims);
        let maxs = ((aabb.maxs - self.origin) / self.scale)
            .map(|x| x.ceil().max(0.0) as u32)
            .inf(&dims);

        self.centers_range(mins.into(), maxs.into())
    }

    pub fn centers(&self) -> impl Iterator<Item = (Point<Real>, VoxelData)> + '_ {
        self.centers_range([0; DIM], self.dimensions)
    }

    pub fn split_with_box(&self, aabb: &Aabb) -> (Option<Self>, Option<Self>) {
        // TODO: optimize this?
        let mut in_box = vec![];
        let mut rest = vec![];
        for (center, voxel) in self.centers() {
            if !voxel.is_empty() {
                if aabb.contains_local_point(&center) {
                    in_box.push(center);
                } else {
                    rest.push(center);
                }
            }
        }

        let in_box = if !in_box.is_empty() {
            Some(Voxels::from_points(&in_box, self.scale))
        } else {
            None
        };

        let rest = if !rest.is_empty() {
            Some(Voxels::from_points(&rest, self.scale))
        } else {
            None
        };

        (in_box, rest)
    }

    #[cfg(feature = "dim2")]
    fn centers_range(
        &self,
        mins: [u32; DIM],
        maxs: [u32; DIM],
    ) -> impl Iterator<Item = (Point<Real>, VoxelData)> + '_ {
        (mins[0]..maxs[0]).flat_map(move |ix| {
            (mins[1]..maxs[1]).map(move |iy| {
                let vid = self.linear_index([ix, iy]);
                let center = self.origin + vector![ix as Real + 0.5, iy as Real + 0.5] * self.scale;
                (center, self.data[vid as usize])
            })
        })
    }

    #[cfg(feature = "dim3")]
    fn centers_range(
        &self,
        mins: [u32; DIM],
        maxs: [u32; DIM],
    ) -> impl Iterator<Item = (Point<Real>, VoxelData)> + '_ {
        (mins[0]..maxs[0]).flat_map(move |ix| {
            (mins[1]..maxs[1]).flat_map(move |iy| {
                (mins[2]..maxs[2]).map(move |iz| {
                    let vid = self.linear_index([ix, iy, iz]);
                    let center = self.origin
                        + Vector::new(ix as Real + 0.5, iy as Real + 0.5, iz as Real + 0.5)
                            * self.scale;
                    (center, self.data[vid as usize])
                })
            })
        })
    }

    #[cfg(feature = "dim2")]
    fn linear_index(&self, voxel_id: [u32; DIM]) -> u32 {
        voxel_id[0] + voxel_id[1] * self.dimensions[0]
    }

    #[cfg(feature = "dim3")]
    fn linear_index(&self, voxel_id: [u32; DIM]) -> u32 {
        voxel_id[0]
            + voxel_id[1] * self.dimensions[0]
            + voxel_id[2] * self.dimensions[0] * self.dimensions[1]
    }

    #[cfg(feature = "dim2")]
    fn to_key(&self, linear_index: u32) -> [u32; DIM] {
        let y = linear_index / self.dimensions[0];
        let x = linear_index % self.dimensions[0];
        [x, y]
    }

    #[cfg(feature = "dim3")]
    fn to_key(&self, linear_index: u32) -> [u32; DIM] {
        let d0d1 = self.dimensions[0] * self.dimensions[1];
        let z = linear_index / d0d1;
        let y = (linear_index - z * d0d1) / self.dimensions[0];
        let x = linear_index % self.dimensions[0];
        [x, y, z]
    }

    fn compute_voxel_data(&self, key: [u32; DIM]) -> VoxelData {
        if self.data[self.linear_index(key) as usize].is_empty() {
            return VoxelData::empty();
        }

        let mut occupied_faces = 0;

        for k in 0..3 {
            let (mut prev, mut next) = (key, key);
            prev[k] = prev[k].saturating_sub(1);
            next[k] += 1;

            if key[k] < self.dimensions[k] - 1
                && !self.data[self.linear_index(next) as usize].is_empty()
            {
                occupied_faces |= 1 << (k * 2);
            }
            if key[k] > 0 && !self.data[self.linear_index(prev) as usize].is_empty() {
                occupied_faces |= 1 << (k * 2 + 1);
            }
        }

        VoxelData(occupied_faces)
    }
}

// NOTE: this code is used to generate the constant tables
// FACES_TO_VOXEL_TYPES, FACES_TO_FEATURE_MASKS, FACES_TO_OCTANT_MASKS.
#[allow(dead_code)]
#[cfg(feature = "dim3")]
fn gen_const_tables() {
    let mut faces_adj_to_vtx = [0usize; 8];
    let mut faces_adj_to_edge = [0usize; 12];
    for fid in 0..6 {
        let vids = Aabb::FACES_VERTEX_IDS[fid];
        let key = 1 << fid;
        faces_adj_to_vtx[vids.0] |= key;
        faces_adj_to_vtx[vids.1] |= key;
        faces_adj_to_vtx[vids.2] |= key;
        faces_adj_to_vtx[vids.3] |= key;
    }

    for eid in 0..12 {
        let evids = Aabb::EDGES_VERTEX_IDS[eid];
        for fid in 0..6 {
            let fvids = Aabb::FACES_VERTEX_IDS[fid];
            if (fvids.0 == evids.0
                || fvids.1 == evids.0
                || fvids.2 == evids.0
                || fvids.3 == evids.0)
                && (fvids.0 == evids.1
                    || fvids.1 == evids.1
                    || fvids.2 == evids.1
                    || fvids.3 == evids.1)
            {
                let key = 1 << fid;
                faces_adj_to_edge[eid] |= key;
            }
        }
    }

    /*
     * FACES_TO_VOXEL_TYPES
     */
    println!("const FACES_TO_VOXEL_TYPES: [VoxelType; 64] = [");
    'outer: for i in 0usize..64 {
        for adjs in faces_adj_to_vtx.iter() {
            if (*adjs & i) == 0 {
                println!("VoxelType::Vertex,");
                continue 'outer;
            }
        }

        for adjs in faces_adj_to_edge.iter() {
            if (*adjs & i) == 0 {
                println!("VoxelType::Edge,");
                continue 'outer;
            }
        }

        for fid in 0..6 {
            if ((1 << fid) & i) == 0 {
                println!("VoxelType::Face,");
                continue 'outer;
            }
        }
    }
    println!("VoxelType::Interior,");
    println!("VoxelType::Empty,");
    println!("];");

    /*
     * FACES_TO_FEATURE_MASKS
     */
    println!("const FACES_TO_FEATURE_MASKS: [u16; 64] = [");
    for i in 0usize..64 {
        // First test if we have vertices.
        let mut vtx_key = 0;
        for (vid, adjs) in faces_adj_to_vtx.iter().enumerate() {
            if (*adjs & i) == 0 {
                vtx_key |= 1 << vid;
            }
        }

        if vtx_key != 0 {
            println!("0b{:b},", vtx_key as u16);
            continue;
        }

        let mut edge_key = 0;
        for (vid, adjs) in faces_adj_to_edge.iter().enumerate() {
            if (*adjs & i) == 0 {
                edge_key |= 1 << vid;
            }
        }

        if edge_key != 0 {
            println!("0b{:b},", edge_key as u16);
            continue;
        }

        let mut face_key = 0;
        for fid in 0..6 {
            if ((1 << fid) & i) == 0 {
                face_key |= 1 << fid;
            }
        }

        if face_key != 0 {
            println!("0b{:b},", face_key as u16);
            continue;
        }
    }

    println!("0b{:b},", u16::MAX);
    println!("0,");
    println!("];");

    /*
     * Faces to octant masks.
     */
    println!("const FACES_TO_OCTANT_MASKS: [u32; 65] = [");
    for i in 0usize..64 {
        // First test if we have vertices.
        let mut octant_mask = 0;
        let mut set_mask = |mask, octant| {
            if (octant_mask >> (octant * 3)) & 0b0111 == 0 {
                octant_mask |= mask << (octant * 3);
            }
        };

        for (vid, adjs) in faces_adj_to_vtx.iter().enumerate() {
            if (*adjs & i) == 0 {
                set_mask(1, vid);
            }
        }

        // This is the index of the axis porting the edges given by
        // Aabb::EDGES_VERTEX_IDS.
        const EDGE_AXIS: [u32; 12] = [0, 1, 0, 1, 0, 1, 0, 1, 2, 2, 2, 2];
        for (eid, adjs) in faces_adj_to_edge.iter().enumerate() {
            if (*adjs & i) == 0 {
                let vid = Aabb::EDGES_VERTEX_IDS[eid];
                // We add OctantPattern::EDGE_X so it matches the values in the OctantPattern enum.
                let mask = EDGE_AXIS[eid] + OctantPattern::EDGE_X as u32;

                set_mask(mask, vid.0);
                set_mask(mask, vid.1);
            }
        }

        // This is the index of the normal of the faces given by
        // Aabb::EDGES_VERTEX_IDS.
        const FACE_NORMALS: [u32; 6] = [0, 0, 1, 1, 2, 2];
        for fid in 0..6 {
            if ((1 << fid) & i) == 0 {
                let vid = Aabb::FACES_VERTEX_IDS[fid];
                // We add OctantPattern::FACE_X so it matches the values in the OctantPattern enum.
                let mask = FACE_NORMALS[fid] + OctantPattern::FACE_X as u32;

                set_mask(mask, vid.0);
                set_mask(mask, vid.1);
                set_mask(mask, vid.2);
                set_mask(mask, vid.3);
            }
        }
        println!("0b{:b},", octant_mask);
    }
    println!("0,");
    println!("];");
}

// Index to the item of FACES_TO_VOXEL_TYPES which identifies interior voxels.
const INTERIOR_FACE_MASK: u8 = 63;
// Index to the item of FACES_TO_VOXEL_TYPES which identifies empty voxels.
const EMPTY_FACE_MASK: u8 = 64;
const FACES_TO_VOXEL_TYPES: [VoxelType; 65] = [
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Face,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Face,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Vertex,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Face,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Face,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Face,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Edge,
    VoxelType::Face,
    VoxelType::Face,
    VoxelType::Face,
    VoxelType::Face,
    VoxelType::Interior,
    VoxelType::Empty,
];

// Deduces the voxel type from the face adjascency information.
// In 3D there are 6 neighbor faces => 64 cases.
const FACES_TO_FEATURE_MASKS: [u16; 65] = [
    0b11111111,
    0b10011001,
    0b1100110,
    0b1010101,
    0b110011,
    0b10001,
    0b100010,
    0b10001,
    0b11001100,
    0b10001000,
    0b1000100,
    0b1000100,
    0b10101010,
    0b10001000,
    0b100010,
    0b110000,
    0b1111,
    0b1001,
    0b110,
    0b101,
    0b11,
    0b1,
    0b10,
    0b1,
    0b1100,
    0b1000,
    0b100,
    0b100,
    0b1010,
    0b1000,
    0b10,
    0b100000,
    0b11110000,
    0b10010000,
    0b1100000,
    0b1010000,
    0b110000,
    0b10000,
    0b100000,
    0b10000,
    0b11000000,
    0b10000000,
    0b1000000,
    0b1000000,
    0b10100000,
    0b10000000,
    0b100000,
    0b10000,
    0b111100000000,
    0b100100000000,
    0b11000000000,
    0b1100,
    0b1100000000,
    0b100000000,
    0b1000000000,
    0b1000,
    0b110000000000,
    0b100000000000,
    0b10000000000,
    0b100,
    0b11,
    0b10,
    0b1,
    0b1111111111111111,
    0,
];

const FACES_TO_OCTANT_MASKS: [u32; 65] = [
    0b1001001001001001001001,
    0b1010010001001010010001,
    0b10001001010010001001010,
    0b10010010010010010010010,
    0b11011001001011011001001,
    0b11111010001011111010001,
    0b111011001010111011001010,
    0b111111010010111111010010,
    0b1001011011001001011011,
    0b1010111011001010111011,
    0b10001011111010001011111,
    0b10010111111010010111111,
    0b11011011011011011011011,
    0b11111111011011111111011,
    0b111011011111111011011111,
    0b111111111111111111111111,
    0b100100100100001001001001,
    0b100110110100001010010001,
    0b110100100110010001001010,
    0b110110110110010010010010,
    0b101101100100011011001001,
    0b101000110100011111010001,
    0b101100110111011001010,
    0b110110111111010010,
    0b100100101101001001011011,
    0b100110000101001010111011,
    0b110100101000010001011111,
    0b110110000000010010111111,
    0b101101101101011011011011,
    0b101000000101011111111011,
    0b101101000111011011111,
    0b111111111111,
    0b1001001001100100100100,
    0b1010010001100110110100,
    0b10001001010110100100110,
    0b10010010010110110110110,
    0b11011001001101101100100,
    0b11111010001101000110100,
    0b111011001010000101100110,
    0b111111010010000000110110,
    0b1001011011100100101101,
    0b1010111011100110000101,
    0b10001011111110100101000,
    0b10010111111110110000000,
    0b11011011011101101101101,
    0b11111111011101000000101,
    0b111011011111000101101000,
    0b111111111111000000000000,
    0b100100100100100100100100,
    0b100110110100100110110100,
    0b110100100110110100100110,
    0b110110110110110110110110,
    0b101101100100101101100100,
    0b101000110100101000110100,
    0b101100110000101100110,
    0b110110000000110110,
    0b100100101101100100101101,
    0b100110000101100110000101,
    0b110100101000110100101000,
    0b110110000000110110000000,
    0b101101101101101101101101,
    0b101000000101101000000101,
    0b101101000000101101000,
    0b0,
    0,
];
