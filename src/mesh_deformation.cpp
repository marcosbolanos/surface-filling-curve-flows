#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/direction_fields.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <polyscope/volume_mesh.h>

#include <igl/predicates/delaunay_triangulation.h>
#include <igl/triangle/cdt.h>
#include <igl/edge_topology.h>
#include <igl/PI.h>
#include <igl/copyleft/tetgen/cdt.h>
#include <igl/copyleft/tetgen/mesh_to_tetgenio.h>
#include <igl/writeOBJ.h>
#include <igl/circulation.h>
#include <igl/harmonic.h>
#include <igl/readDMAT.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/AABB.h>
#include <igl/barycentric_coordinates.h>
#include <Eigen/Eigenvalues>

#include <tetgen-src/tetgen.h>
// #include <igl/copyleft/cgal/delaunay_triangulation.h>

#include <queue>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include "args/args.hxx"
#include "imgui.h"

#include "modules/scene_file.h"
#include "modules/curve_evolution.h"
#include "modules/medial_axis_2d.h"
#include "modules/evolution_mode.h"
#include "modules/curve_remesh_2d.h"
#include "modules/circumcenter.h"
#include "modules/surface_point_to_cartesian.h"
#include "modules/read_curve_on_mesh.h"
#include "modules/connect_surface_points.h"
#include "modules/remesh_curve_on_surface.h"
#include "modules/surface_filling_energy_geodesic.h"
#include "modules/surface_path_evolution.h"
#include "modules/get_tangent_basis.h"
#include "modules/surface_point_tangent_basis.h"
#include "modules/write_curve.h"
#include "modules/curve_to_arclength.h"
#include "modules/curve_curvature.h"
#include "modules/curvature_barrier_threshold.h"
#include "modules/cut_mesh_with_curve.h"
#include "modules/gc_igl_convert.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

using namespace modules;

float radius = 30;
float timestep = 1;
float h = igl::PI * radius / 25;
float rmax = 0.5;
// double branchRatio = 1.18;

// Polyscope visualization handle, to quickly add data to the surface
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

modules::SceneObject scene;

std::vector<SurfacePoint> nodes = {}, restNodes = {}, initialNodes = {};
std::vector<bool> isFixedNode = {}, isInitialFixedNode = {};
std::vector<std::array<int, 2>> segments = {}, restSegments = {}, initialSegments = {};

std::vector<std::vector<SurfacePoint>> segmentSurfacePoints = {}, restSegmentSurfacePoints = {}, initialSegmentSurfacePoints = {};
std::vector<double> segmentLengths = {}, restSegmentLengths = {}, initialSegmentLengths = {};

VertexData<Vector2> vField;
VertexData<double> smoothedFunction;

double bc_frac = 1.0;
double bc_dir = -0.03;
Eigen::MatrixXd V,U,V_bc,U_bc;
Eigen::VectorXd Z;
Eigen::MatrixXi F;
Eigen::VectorXi b;

// Example computation function -- this one computes and registers a scalar
// quantity
int iteration = 0;
bool writeCurve = false;
bool writeData = false;
bool runLoop = false;
double adaptiveCurvatureBarrierKmax = -1.0;

std::ofstream data_out;

std::tuple<
  std::vector<Vector3> /* nodes */,
  std::vector<std::array<int, 2>> /* segments */,
  std::vector<Vector3> /* cartesianCoords of the nodes */
> formatGeodesicCurve(
  ManifoldSurfaceMesh &_mesh,
  VertexPositionGeometry &_geometry,
  const std::vector<SurfacePoint> &_nodes,
  const std::vector<std::array<int, 2>> &_segments,
  const std::vector<std::vector<SurfacePoint>> &_segmentSurfacePoints
) {
  std::vector<Vector3> cartesianCoords = modules::surface_point_to_cartesian(_mesh, _geometry, _nodes);

  std::vector<Vector3> allNodes = cartesianCoords;
  std::vector<std::array<int, 2>> allSegments = {};

  for (int i = 0; i < _segments.size(); i++) {
    int v0 = _segments[i][0], v1 = _segments[i][1];

    auto cc = surface_point_to_cartesian(_mesh, _geometry, _segmentSurfacePoints[i]);

    std::vector<int> ccIdx = {};
    ccIdx.emplace_back(v0);
    for (int j = 0; j < cc.size(); j++) {
      allNodes.emplace_back(cc[j]);
      ccIdx.emplace_back(allNodes.size() - 1);
    }
    ccIdx.emplace_back(v1);

    for (int j = 0; j < ccIdx.size() - 1; j++) {
      std::array<int, 2> line = {ccIdx[j], ccIdx[j + 1]};
      allSegments.push_back(line);
    }
  }

  return {
    allNodes,
    allSegments,
    cartesianCoords
  };
}

std::tuple<
  std::vector<Vector3> /* nodes */,
  std::vector<std::array<int, 2>> /* segments */,
  std::vector<double> /* norm of the vector */
> formatVector(
  std::vector<Vector3> &cartesianCoords,
  std::vector<Vector3> &vector
) {
  assert(cartesianCoords.size() == vector.size());

  std::vector<Vector3> dNodes = cartesianCoords;
  std::vector<std::array<int, 2>> dSegments = {};
  std::vector<double> dNorm = {};

  for (int i = 0; i < vector.size(); i++) {
    dNodes.emplace_back(cartesianCoords[i] + vector[i]);
    dNorm.emplace_back(norm(vector[i]));
  }

  for (int i = 0; i < vector.size(); i++) {
    int oi = i + vector.size();
    std::array<int, 2> line = {i, oi};
    dSegments.push_back(line);
  }

  return {
    dNodes,
    dSegments,
    dNorm
  };
}

std::tuple<
  std::vector<Vector3> /* nodes */,
  std::vector<std::array<int, 2>> /* segments */
> formatGeodesicPaths(
  ManifoldSurfaceMesh &_mesh,
  VertexPositionGeometry &_geometry,
  std::vector<std::vector<SurfacePoint>> &paths
) {
  std::vector<Vector3> nodes = {};
  std::vector<std::array<int, 2>> segments = {};

  for (int i = 0; i < paths.size(); i++) {
    auto path = paths[i];

    if (path.size() == 0) {
      continue;
    }

    auto cc = surface_point_to_cartesian(_mesh, _geometry, path);

    for (int j = 0; j < path.size() - 1; j++) {
      int v0 = nodes.size();
      nodes.emplace_back(cc[j]);
      int v1 = nodes.size();
      nodes.emplace_back(cc[j + 1]);

      std::array<int, 2> line = {v0, v1};
      segments.push_back(line);
    }
  }

  return {
    nodes,
    segments
  };
}

std::tuple<
  std::vector<Vector3> /* nodes */,
  std::vector<double> /* radius */
> formatMedialAxis(
  const std::vector<std::vector<Vector3>> &medialAxis,
  const std::vector<Vector3> &cartesianCoords
){
  std::vector<Vector3> maNodes = {};
  std::vector<double> maRadius = {};

  for (int i = 0; i < medialAxis.size(); i++) {
    assert(medialAxis[i].size() == 2);

    for (int j = 0; j < medialAxis[i].size(); j++) {
      maNodes.emplace_back(medialAxis[i][j]);

      double r = norm(medialAxis[i][j] - cartesianCoords[i]);
      maRadius.emplace_back(r);
    }
  }

  return {
    maNodes,
    maRadius
  };
}

SurfacePoint project_to_surface_point(
  ManifoldSurfaceMesh& _mesh,
  VertexPositionGeometry& _geometry,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  igl::AABB<Eigen::MatrixXd, 3>& tree,
  const Vector3& query
) {
  Eigen::MatrixXd P(1, 3);
  P << query.x, query.y, query.z;

  Eigen::VectorXd sqrD;
  Eigen::VectorXi I;
  Eigen::MatrixXd C;
  tree.squared_distance(V, F, P, sqrD, I, C);

  int fid = I(0);
  Eigen::MatrixXd A(1, 3), B(1, 3), D(1, 3);
  A.row(0) = V.row(F(fid, 0));
  B.row(0) = V.row(F(fid, 1));
  D.row(0) = V.row(F(fid, 2));

  Eigen::MatrixXd bary;
  igl::barycentric_coordinates(C, A, B, D, bary);

  Eigen::Vector3d bcoord = bary.row(0).transpose();
  for (int i = 0; i < 3; i++) {
    bcoord(i) = std::max(0.0, bcoord(i));
  }
  double sum = bcoord.sum();
  if (sum < 1e-12) {
    bcoord = Eigen::Vector3d(1.0, 0.0, 0.0);
  } else {
    bcoord /= sum;
  }

  Vector3 gcBary {bcoord(0), bcoord(1), bcoord(2)};
  return SurfacePoint(_mesh.face(fid), gcBary);
}

void initialize_parallel_rings(
  ManifoldSurfaceMesh& _mesh,
  VertexPositionGeometry& _geometry,
  const modules::SceneObject& _scene,
  const double rminWorld,
  std::vector<SurfacePoint>& _nodes,
  std::vector<std::array<int, 2>>& _segments,
  std::vector<bool>& _isFixedNode
) {
  int ringCount = std::max(1, _scene.initTriangleCount);
  int ringVertices = std::max(3, _scene.initRingVertexCount);

  double nmPerUnit = std::max(_scene.nanometersPerUnit, 1e-12);
  double spacing = _scene.initRingSpacing;
  if (_scene.initRingSpacingNm > 0.0) {
    spacing = _scene.initRingSpacingNm / nmPerUnit;
  }
  if (spacing <= 0.0) {
    spacing = std::max(rminWorld / 3.0, 1e-3);
  }

  double outerRadius = _scene.initRingRadius;
  if (_scene.initRingRadiusNm > 0.0) {
    outerRadius = _scene.initRingRadiusNm / nmPerUnit;
  }
  if (outerRadius <= 0.0) {
    outerRadius = std::max(3.0 * rminWorld, (ringCount + 1.0) * spacing);
  }

  Eigen::MatrixXd V(_mesh.nVertices(), 3);
  for (int i = 0; i < static_cast<int>(_mesh.nVertices()); i++) {
    auto v = _mesh.vertex(i);
    auto p = _geometry.vertexPositions[v];
    V.row(i) << p.x, p.y, p.z;
  }

  Eigen::MatrixXi F(_mesh.nFaces(), 3);
  for (int i = 0; i < static_cast<int>(_mesh.nFaces()); i++) {
    auto f = _mesh.face(i);
    int j = 0;
    for (auto v : f.adjacentVertices()) {
      F(i, j++) = v.getIndex();
    }
  }

  igl::AABB<Eigen::MatrixXd, 3> tree;
  tree.init(V, F);

  Vector3 centroid {0.0, 0.0, 0.0};
  int vCount = 0;
  for (auto v : _mesh.vertices()) {
    centroid += _geometry.vertexPositions[v];
    vCount++;
  }
  centroid /= static_cast<double>(std::max(vCount, 1));

  Eigen::Vector3d mean(centroid.x, centroid.y, centroid.z);
  Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
  for (int i = 0; i < V.rows(); i++) {
    Eigen::Vector3d d = V.row(i).transpose() - mean;
    covariance += d * d.transpose();
  }
  covariance /= std::max(1, static_cast<int>(V.rows()));

  std::string axisMode = _scene.initStackAxisMode;
  std::transform(axisMode.begin(), axisMode.end(), axisMode.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

  int axisColumn = 2;
  if (axisMode == "short" || axisMode == "lowest" || axisMode == "min" || axisMode == "smallest") {
    axisColumn = 0;
  } else if (axisMode == "mid" || axisMode == "middle" || axisMode == "median") {
    axisColumn = 1;
  }

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(covariance);
  Eigen::Vector3d axisEig = eig.eigenvectors().col(axisColumn).normalized();
  Vector3 axis {axisEig(0), axisEig(1), axisEig(2)};

  Vector3 ref {_scene.initStackRefX, _scene.initStackRefY, _scene.initStackRefZ};
  if (norm(ref) < 1e-8) {
    ref = Vector3{0.0, 0.0, 1.0};
  }
  ref = unit(ref);
  if (std::abs(dot(axis, ref)) > 0.95) {
    ref = std::abs(dot(axis, Vector3{1.0, 0.0, 0.0})) < 0.95 ? Vector3{1.0, 0.0, 0.0} : Vector3{0.0, 1.0, 0.0};
  }

  Vector3 tangentX = unit(cross(ref, axis));
  Vector3 tangentY = unit(cross(axis, tangentX));

  double tilt = _scene.initStackTiltDeg * igl::PI / 180.0;
  Vector3 stackAxis = unit(std::cos(tilt) * axis + std::sin(tilt) * tangentX);

  Vector3 refAfterTilt = std::abs(dot(stackAxis, ref)) < 0.95 ? ref : tangentY;
  Vector3 basisX = unit(cross(refAfterTilt, stackAxis));
  Vector3 basisY = unit(cross(stackAxis, basisX));

  double rotate = _scene.initStackRotateDeg * igl::PI / 180.0;
  Vector3 rotatedX = std::cos(rotate) * basisX + std::sin(rotate) * basisY;
  Vector3 rotatedY = -std::sin(rotate) * basisX + std::cos(rotate) * basisY;

  double meshRadiusAroundAxis = 0.0;
  for (auto v : _mesh.vertices()) {
    Vector3 p = _geometry.vertexPositions[v];
    Vector3 rel = p - centroid;
    double along = dot(rel, stackAxis);
    Vector3 radial = rel - along * stackAxis;
    meshRadiusAroundAxis = std::max(meshRadiusAroundAxis, norm(radial));
  }

  double projectedRingRadius = outerRadius;
  if (_scene.initRingRadius <= 0.0 && _scene.initRingRadiusNm <= 0.0) {
    projectedRingRadius = std::max(meshRadiusAroundAxis * 1.1, spacing);
  }

  double offsetCenter = 0.5 * static_cast<double>(ringCount - 1);

  for (int ring = 0; ring < ringCount; ring++) {
    double ringOffset = (static_cast<double>(ring) - offsetCenter) * spacing;
    Vector3 ringCenter = centroid + ringOffset * stackAxis;
    int baseNode = _nodes.size();
    double phase = (ring % 2 == 0) ? 0.0 : (igl::PI / static_cast<double>(ringVertices));

    for (int k = 0; k < ringVertices; k++) {
      double angle = phase + 2.0 * igl::PI * static_cast<double>(k) / static_cast<double>(ringVertices);
      Vector3 dir = std::cos(angle) * rotatedX + std::sin(angle) * rotatedY;
      Vector3 query = ringCenter + projectedRingRadius * dir;

      SurfacePoint sp = project_to_surface_point(_mesh, _geometry, V, F, tree, query);
      _nodes.emplace_back(sp);
      _isFixedNode.emplace_back(false);
    }

    for (int k = 0; k < ringVertices; k++) {
      _segments.emplace_back(std::array<int, 2> {baseNode + k, baseNode + (k + 1) % ringVertices});
    }
  }
}

int count_curve_components(
  const size_t nodeCount,
  const std::vector<std::array<int, 2>>& _segments
) {
  if (nodeCount == 0) {
    return 0;
  }

  std::vector<int> parent(nodeCount);
  std::iota(parent.begin(), parent.end(), 0);
  std::vector<bool> active(nodeCount, false);

  std::function<int(int)> findRoot = [&](int x) {
    if (parent[x] != x) {
      parent[x] = findRoot(parent[x]);
    }
    return parent[x];
  };

  auto unite = [&](int a, int b) {
    int ra = findRoot(a);
    int rb = findRoot(b);
    if (ra != rb) {
      parent[rb] = ra;
    }
  };

  for (auto segment : _segments) {
    int a = segment[0];
    int b = segment[1];
    if (a < 0 || b < 0 || a >= static_cast<int>(nodeCount) || b >= static_cast<int>(nodeCount)) {
      continue;
    }
    active[a] = true;
    active[b] = true;
    unite(a, b);
  }

  std::unordered_set<int> roots;
  for (int i = 0; i < static_cast<int>(nodeCount); i++) {
    if (active[i]) {
      roots.insert(findRoot(i));
    }
  }

  return static_cast<int>(roots.size());
}

std::vector<std::vector<int>> curve_components(
  const size_t nodeCount,
  const std::vector<std::array<int, 2>>& _segments
) {
  if (nodeCount == 0) {
    return {};
  }

  std::vector<int> parent(nodeCount);
  std::iota(parent.begin(), parent.end(), 0);
  std::vector<bool> active(nodeCount, false);

  std::function<int(int)> findRoot = [&](int x) {
    if (parent[x] != x) {
      parent[x] = findRoot(parent[x]);
    }
    return parent[x];
  };

  auto unite = [&](int a, int b) {
    int ra = findRoot(a);
    int rb = findRoot(b);
    if (ra != rb) {
      parent[rb] = ra;
    }
  };

  for (auto segment : _segments) {
    int a = segment[0];
    int b = segment[1];
    if (a < 0 || b < 0 || a >= static_cast<int>(nodeCount) || b >= static_cast<int>(nodeCount)) {
      continue;
    }
    active[a] = true;
    active[b] = true;
    unite(a, b);
  }

  std::unordered_map<int, int> rootToComp;
  std::vector<std::vector<int>> components;
  for (int i = 0; i < static_cast<int>(nodeCount); i++) {
    if (!active[i]) {
      continue;
    }
    int root = findRoot(i);
    if (rootToComp.find(root) == rootToComp.end()) {
      rootToComp[root] = components.size();
      components.emplace_back(std::vector<int>{});
    }
    components[rootToComp[root]].emplace_back(i);
  }

  return components;
}

bool insert_dynamic_outer_ring(
  ManifoldSurfaceMesh& _mesh,
  VertexPositionGeometry& _geometry,
  const modules::SceneObject& _scene,
  const int _iteration,
  std::vector<SurfacePoint>& _nodes,
  std::vector<std::array<int, 2>>& _segments,
  std::vector<std::vector<SurfacePoint>>& _segmentSurfacePoints,
  std::vector<double>& _segmentLengths,
  std::vector<bool>& _isFixedNode,
  const double fallbackSpacing
) {
  if (!_scene.ringInsertionEnabled || _scene.ringInsertionMax <= 0 || _scene.ringInsertionEvery <= 0) {
    return false;
  }
  if (_iteration <= _scene.ringInsertionAfterIter) {
    return false;
  }
  if ((_iteration - _scene.ringInsertionAfterIter) % _scene.ringInsertionEvery != 0) {
    return false;
  }

  auto components = curve_components(_nodes.size(), _segments);
  int currentComponents = static_cast<int>(components.size());
  if (currentComponents >= _scene.ringInsertionMax || _nodes.empty()) {
    return false;
  }

  double nmPerUnit = std::max(_scene.nanometersPerUnit, 1e-12);
  double insertionSpacing = _scene.ringInsertionSpacing;
  if (_scene.ringInsertionSpacingNm > 0.0) {
    insertionSpacing = _scene.ringInsertionSpacingNm / nmPerUnit;
  }
  if (insertionSpacing <= 0.0) {
    insertionSpacing = std::max(fallbackSpacing, 1e-3);
  }

  int ringVertices = std::max(3, _scene.initRingVertexCount);

  Eigen::MatrixXd V(_mesh.nVertices(), 3);
  for (int i = 0; i < static_cast<int>(_mesh.nVertices()); i++) {
    auto v = _mesh.vertex(i);
    auto p = _geometry.vertexPositions[v];
    V.row(i) << p.x, p.y, p.z;
  }

  Eigen::MatrixXi F(_mesh.nFaces(), 3);
  for (int i = 0; i < static_cast<int>(_mesh.nFaces()); i++) {
    auto f = _mesh.face(i);
    int j = 0;
    for (auto v : f.adjacentVertices()) {
      F(i, j++) = v.getIndex();
    }
  }

  igl::AABB<Eigen::MatrixXd, 3> tree;
  tree.init(V, F);

  auto cartesianNodes = modules::surface_point_to_cartesian(_mesh, _geometry, _nodes);
  Vector3 centroid {0.0, 0.0, 0.0};
  for (auto p : cartesianNodes) {
    centroid += p;
  }
  centroid /= static_cast<double>(std::max(static_cast<int>(cartesianNodes.size()), 1));

  SurfacePoint centerSp = project_to_surface_point(_mesh, _geometry, V, F, tree, centroid);
  Vector3 center = modules::surface_point_to_cartesian(_mesh, _geometry, {centerSp})[0];
  auto [basisX, basisY] = modules::surface_point_tangent_basis(_geometry, centerSp);

  std::vector<double> radialDistance(_nodes.size(), 0.0);
  double outerRadius = 0.0;
  for (int i = 0; i < static_cast<int>(cartesianNodes.size()); i++) {
    auto p = cartesianNodes[i];
    Vector3 rel = p - center;
    double x = dot(rel, basisX);
    double y = dot(rel, basisY);
    double r = std::sqrt(x * x + y * y);
    radialDistance[i] = r;
    outerRadius = std::max(outerRadius, r);
  }

  std::vector<double> componentRadii;
  componentRadii.reserve(components.size());
  for (const auto& comp : components) {
    double meanRadius = 0.0;
    for (int nodeId : comp) {
      meanRadius += radialDistance[nodeId];
    }
    meanRadius /= static_cast<double>(std::max(static_cast<int>(comp.size()), 1));
    componentRadii.emplace_back(meanRadius);
  }
  std::sort(componentRadii.begin(), componentRadii.end(), std::greater<double>());

  double ringRadius = std::max(outerRadius + insertionSpacing, insertionSpacing);
  if (componentRadii.size() >= 2) {
    double largestGap = -1.0;
    int gapIdx = -1;
    for (int i = 0; i + 1 < static_cast<int>(componentRadii.size()); i++) {
      double gap = componentRadii[i] - componentRadii[i + 1];
      if (gap > largestGap) {
        largestGap = gap;
        gapIdx = i;
      }
    }
    if (gapIdx >= 0 && largestGap > 1e-6) {
      ringRadius = componentRadii[gapIdx + 1] + 0.5 * largestGap;
    }
  }
  int baseNode = _nodes.size();
  double phase = (currentComponents % 2 == 0) ? 0.0 : (igl::PI / static_cast<double>(ringVertices));

  for (int k = 0; k < ringVertices; k++) {
    double angle = phase + 2.0 * igl::PI * static_cast<double>(k) / static_cast<double>(ringVertices);
    Vector3 dir = std::cos(angle) * basisX + std::sin(angle) * basisY;
    Vector3 query = center + ringRadius * dir;
    SurfacePoint sp = project_to_surface_point(_mesh, _geometry, V, F, tree, query);
    _nodes.emplace_back(sp);
    _isFixedNode.emplace_back(false);
  }

  for (int k = 0; k < ringVertices; k++) {
    _segments.emplace_back(std::array<int, 2> {baseNode + k, baseNode + (k + 1) % ringVertices});
  }

  std::tie(_segmentSurfacePoints, _segmentLengths) = modules::connect_surface_points(_mesh, _geometry, _nodes, _segments);

  std::cout << "inserted dynamic outer ring at iteration " << _iteration
            << ": component count " << currentComponents << " -> " << (currentComponents + 1)
            << ", radius " << ringRadius << std::endl;
  return true;
}

std::tuple<
  Eigen::MatrixXd /* V */,
  Eigen::MatrixXi /* T */
> tetrahedralization(
  ManifoldSurfaceMesh &_mesh,
  VertexPositionGeometry &_geometry,
  std::vector<Vector3> &cartesianPaths
) {
  Eigen::MatrixXd V(cartesianPaths.size(), 3);

  for (int i = 0; i < cartesianPaths.size(); i++) {
    V(i, 0) = cartesianPaths[i].x;
    V(i, 1) = cartesianPaths[i].y;
    V(i, 2) = cartesianPaths[i].z;
  }

  auto start = std::chrono::high_resolution_clock::now();

  Eigen::MatrixXi T;

  {
    // Eigen::MatrixXd T_V;
    // Eigen::MatrixXi _, __;
    // std::string flags = "S0cQ";
    // igl::copyleft::tetgen::tetrahedralize(V, _, flags, T_V, T, __);
  }

  {
    tetgenio in, out;

    in.firstnumber = 0;

    in.numberofpoints = V.rows();
    in.pointlist = new REAL[in.numberofpoints * 3];
    // loop over points
    for(int i = 0; i < V.rows(); i++)
    {
      in.pointlist[i*3+0] = V(i, 0);
      in.pointlist[i*3+1] = V(i, 1);
      in.pointlist[i*3+2] = V(i, 2);
    }
    in.numberoffacets = 0;
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];

    char* flags = "S0cQ";

    tetrahedralize(flags, &in, &out);

    assert(out.numberofpoints == V.rows());
    assert(out.numberofcorners == 4);

    T.resize(out.numberoftetrahedra, 4);
    int min_index = 1e7;
    int max_index = -1e7;
    // loop over tetrahedra
    for(int i = 0; i < out.numberoftetrahedra; i++)
    {
      for(int j = 0; j<out.numberofcorners; j++)
      {
        int index = out.tetrahedronlist[i * out.numberofcorners + j];
        T(i, j) = index;
        min_index = (min_index > index ? index : min_index);
        max_index = (max_index < index ? index : max_index);
      }
    }
  }

  auto end = std::chrono::high_resolution_clock::now();

  std::cout << "tetrahedralize: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

  return {
    V,
    T
  };
}

void doWork() {
  iteration++;

  insert_dynamic_outer_ring(
    *mesh,
    *geometry,
    scene,
    iteration,
    nodes,
    segments,
    segmentSurfacePoints,
    segmentLengths,
    isFixedNode,
    std::max(radius / 3.0, 1e-3)
  );

  restNodes = nodes;
  restSegments = segments;
  restSegmentSurfacePoints = segmentSurfacePoints;
  restSegmentLengths = segmentLengths;

  std::cout << "===== iteration: " << iteration << "=====" << std::endl;

  std::cout << "numNodes: " << nodes.size() << std::endl;
  std::cout << "numSegments: " << segments.size() << std::endl;
  std::cout << std::endl;

  auto start = std::chrono::high_resolution_clock::now();

  auto cartesianCoords = modules::surface_point_to_cartesian(*mesh, *geometry, nodes);

  double effectiveCurvatureBarrierKmax = scene.curvatureBarrierThreshold;
  if (scene.curvatureConstraint == "barrier" && scene.w_curvatureBarrier > 0.0) {
    auto threshold = modules::curvature_barrier_threshold(cartesianCoords, segments, segmentLengths, scene.curvatureBarrierMinLength);
    double feasibleKmax = threshold.minFeasibleKmax;
    double margin = std::max(scene.curvatureBarrierEpsilon, 1e-6);
    double targetKmax = scene.curvatureBarrierThreshold;
    double feasibleFloor = feasibleKmax + std::max(10.0 * margin, 1e-4);

    if (adaptiveCurvatureBarrierKmax < 0.0) {
      double initOffset = std::max(0.05 * std::max(feasibleKmax, 1.0), 2.0 * margin);
      adaptiveCurvatureBarrierKmax = std::max(targetKmax, feasibleKmax + initOffset);
    }

    double ramp = std::max(0.02 * adaptiveCurvatureBarrierKmax, 5.0 * margin);
    if (adaptiveCurvatureBarrierKmax > targetKmax) {
      adaptiveCurvatureBarrierKmax = std::max(targetKmax, adaptiveCurvatureBarrierKmax - ramp);
    }

    if (adaptiveCurvatureBarrierKmax < feasibleFloor) {
      adaptiveCurvatureBarrierKmax = feasibleFloor;
    }

    effectiveCurvatureBarrierKmax = adaptiveCurvatureBarrierKmax;

    std::cout << "curvature barrier threshold info: feasible k_max >= " << feasibleKmax
              << ", effective k_max = " << effectiveCurvatureBarrierKmax
              << ", target k_max = " << targetKmax
              << ", valid nodes = " << threshold.validNodeCount << std::endl;
  }

  std::vector<Vector3> d, g;
  double f;
  std::vector<std::vector<Vector3>> medialAxis;

  // 1. get the descent direction
  if (runLoop) {
    modules::SurfaceFillingEnergy::Options options = {
      .p = scene.p,
      .q = scene.q,
      .w_fieldAlignedness = scene.w_fieldAlignedness,
      .vectorField = vField,
      .w_curvatureAlignedness = scene.w_curvatureAlignedness,
      .useCurvatureBarrier = scene.curvatureConstraint == "barrier",
      .w_curvatureBarrier = scene.w_curvatureBarrier,
      .curvatureBarrierThreshold = effectiveCurvatureBarrierKmax,
      .curvatureBarrierEpsilon = scene.curvatureBarrierEpsilon,
      .curvatureBarrierMinLength = scene.curvatureBarrierMinLength,
      .w_bilaplacian = scene.w_bilaplacian,
      .useAnisotropicAlphaOnMesh = scene.varyingAlpha,
      .alphaRatioOnMesh = smoothedFunction,
      .useGeodesicMedialAxis = scene.useGeodesicMedialAxis
    };
    std::tie(d, g, f, medialAxis) = modules::surface_filling_energy_geodesic(*mesh, *geometry, nodes, segments, segmentSurfacePoints, segmentLengths, isFixedNode, cartesianCoords, radius, rmax, options);

    // 2. evolve the curve but be careful not to cause self-intersection
    std::vector<std::vector<SurfacePoint>> retractionPaths = {};
    std::tie(nodes, segments, segmentSurfacePoints, segmentLengths, isFixedNode, retractionPaths) = modules::surface_path_evolution(*mesh, *geometry, nodes, segments, segmentSurfacePoints, segmentLengths, isFixedNode, h, d);
  } else {
    // std::tie(segmentSurfacePoints, segmentLengths) = modules::connect_surface_points(*mesh, *geometry, nodes, segments);
  }

  auto updatedCartesianCoords = modules::surface_point_to_cartesian(*mesh, *geometry, nodes);
  auto curvatureData = modules::curve_curvature(updatedCartesianCoords, segments);
  double maxCurvature = curvatureData.maxCurvature;

  auto end = std::chrono::high_resolution_clock::now();

  int time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  std::cout << "total time: " << time << "ms" << std::endl;
  std::cout << "max curvature: " << maxCurvature << std::endl;
  std::cout << "===== iteration end =====" << std::endl << std::endl;

  // visualize
  {
    auto [prevNodes, prevSegments, prevCartesianCoords] = formatGeodesicCurve(*mesh, *geometry, restNodes, restSegments, restSegmentSurfacePoints);
    polyscope::registerCurveNetwork("prev curve", prevNodes, prevSegments)
      ->setEnabled(false)
    ;

    auto [allNodes, allSegments, newCartesianCoords] = formatGeodesicCurve(*mesh, *geometry, nodes, segments, segmentSurfacePoints);

    auto arclength = modules::curve_to_arclength(allNodes, allSegments);

    polyscope::registerCurveNetwork("curve", allNodes, allSegments)
      ->addNodeScalarQuantity("arclength", arclength)
      // ->setColor({1.0, 1.0, 1.0})
      // ->setRadius(0.003)
    ;

    // polyscope::registerPointCloud("initial nodes", allNodes);

    std::vector<Vector3> fixedNodes = {}, freeNodes = {};
    for (int i = 0; i < nodes.size(); i++) {
      if (isFixedNode[i]) {
        fixedNodes.emplace_back(newCartesianCoords[i]);
      } else {
        freeNodes.emplace_back(newCartesianCoords[i]);
      }
    }

    if (writeCurve) {
      polyscope::registerPointCloud("free nodes", freeNodes)
        // ->setPointColor({1.0, 1.0, 1.0})
      ;
      polyscope::registerPointCloud("fixed nodes", fixedNodes)
        // ->setPointColor({.5, .5, .5})
      ;

      std::string meshfname = "./surfaces/surface_" + std::to_string(iteration) + ".obj";
      std::ofstream meshf(meshfname);
      for (int i = 0; i < mesh->nVertices(); i++) {
        auto v = geometry->vertexPositions[i];
        meshf << "v " << v.x << " " << v.y << " " << v.z << std::endl;
      }
      meshf << std::endl;
      for (int i = 0; i < mesh->nFaces(); i++) {
        auto f = mesh->face(i);
        meshf << "f " << f.halfedge().vertex().getIndex() + 1 << " " << f.halfedge().next().vertex().getIndex() + 1 << " " << f.halfedge().next().next().vertex().getIndex() + 1 << std::endl;
      }
      // meshf.close();

      // auto [dNodes, dSegments, dNorm] = formatVector(cartesianCoords, d);
      // auto [gNodes, gSegments, gNorm] = formatVector(cartesianCoords, g);

      // auto dn = polyscope::registerCurveNetwork("descent direction", dNodes, dSegments);
      // dn->setEnabled(false);
      // dn->addEdgeScalarQuantity("norm", dNorm);

      // auto gn = polyscope::registerCurveNetwork("gradient direction", gNodes, gSegments);
      // gn->setEnabled(false);
      // gn->addEdgeScalarQuantity("norm", gNorm);

      // auto [maNodes, maRadius] = formatMedialAxis(medialAxis, cartesianCoords);
      // auto pma = polyscope::registerPointCloud("medial axis", maNodes);
      // pma->setEnabled(false);
      // pma->addScalarQuantity("radius", maRadius);

      // std::vector<Vector3> normalDirections, tangentX, tangentY;
      // for (auto node : restNodes) {
      //   auto [x, y] = modules::get_tangent_basis(*geometry, node);
      //   tangentX.emplace_back(x / 30);
      //   tangentY.emplace_back(y / 30);
      //   normalDirections.emplace_back(cross(x, y) / 30);
      // }

      // auto [nn, ns, _] = formatVector(cartesianCoords, normalDirections);
      // polyscope::registerCurveNetwork("normal directions", nn, ns)->setEnabled(false);

      // auto [txn, txs, __] = formatVector(cartesianCoords, tangentX);
      // polyscope::registerCurveNetwork("tangent x", txn, txs)->setEnabled(false);

      // auto [tyn, tys, ___] = formatVector(cartesianCoords, tangentY);
      // polyscope::registerCurveNetwork("tangent y", tyn, tys)->setEnabled(false);

      // int ten = segments.size();
      // std::vector<std::array<int, 2>> firstTen(ten);
      // std::copy(segments.begin(), segments.begin() + ten, firstTen.begin());

      std::string fname = "./objs/curve_" + std::to_string(iteration) + ".obj";
      modules::write_curve(fname, allNodes, allSegments);

      // std::string radiusfilename = "./radii/radii_" + std::to_string(iteration - 1) + ".csv";
      // std::ofstream radiusfile(radiusfilename);

      // radiusfile << "radius, weight, alpha" << std::endl;

      // std::vector<double> nodeWeight(cartesianCoords.size(), 0);
      // for (int i = 0; i < restSegments.size(); i++) {
      //   auto [v0, v1] = restSegments[i];
      //   nodeWeight[v0] += restSegmentLengths[i] / 2;
      //   nodeWeight[v1] += restSegmentLengths[i] / 2;
      // }

      // double branchRatio = std::sqrt(std::sqrt(2));
      // double branchRadius = radius * branchRatio;
      // double alpha = 4 / (branchRadius * branchRadius);

      // for (int i = 0; i < cartesianCoords.size(); i++) {
      //   auto v = cartesianCoords[i];
      //   for (int j = 0; j < 2; j++) {
      //     auto c = medialAxis[i][j];
      //     double radius = norm(v - c);
      //     radiusfile << radius << ", " << nodeWeight[i] << ", " << alpha << std::endl;
      //   }
      // }

      // std::string _fname = "./fixed_points/fixed_points_" + std::to_string(iteration) + ".obj";
      // std::ofstream _myfile(_fname);

      // for (auto v : fixedNodes) {
      //   _myfile << "v " << v.x << " " << v.y << " " << v.z << std::endl;
      // }

      // std::string _fname2 = "./medial_axis/medial_axis_" + std::to_string(iteration) + ".obj";
      // std::ofstream _myfile2(_fname2);

      // for (auto v : maNodes) {
      //   _myfile2 << "v " << v.x << " " << v.y << " " << v.z << std::endl;
      // }

      // std::string _fname3 = "./descent/descent_" + std::to_string(iteration) + ".obj";
      // std::ofstream _myfile3(_fname3);

      // for (auto d : dNodes) {
      //   _myfile3 << "v " << d.x << " " << d.y << " " << d.z << std::endl;
      // }
      // for (auto d : dSegments) {
      //   _myfile3 << "l " << d[0]+1 << " " << d[1]+1 << std::endl;
      // }

      // auto [cutVs, cutFs, _V, _F] = modules::cut_mesh_with_curve(*mesh, *geometry, nodes, segments, segmentSurfacePoints);
      // for (int i = 0; i < cutVs.size(); i++) {
      //   std::string _fname = "./cut_mesh/cut_mesh_" + std::to_string(iteration) + "_" + std::to_string(i) + ".obj";
      //   std::ofstream _myfile4(_fname);

      //   auto _V = cutVs[i];
      //   auto _F = cutFs[i];

      //   for (int j = 0; j < _V.rows(); j++) {
      //     _myfile4 << "v " << _V(j, 0) << " " << _V(j, 1) << " " << _V(j, 2) << std::endl;
      //   }

      //   for (int j = 0; j < _F.rows(); j++) {
      //     _myfile4 << "f " << _F(j, 0) + 1 << " " << _F(j, 1) + 1 << " " << _F(j, 2) + 1 << std::endl;
      //   }
      // }
    }

    if (writeData) {
      Eigen::VectorXd _d(d.size() * 3);
      _d.setZero();
      for (int i = 0; i < d.size(); i++) {
        _d(i * 3 + 0) = d[i].x;
        _d(i * 3 + 1) = d[i].y;
        _d(i * 3 + 2) = d[i].z;
      }

      double l1 = _d.lpNorm<1>();
      double l2 = _d.lpNorm<2>();
      double linf = _d.lpNorm<Eigen::Infinity>();

      // iteration, time, numNodes, f, descent norm (L2), descent norm (L1), descent norm (L∞), max curvature
      data_out << iteration << ", " << time << ", " << nodes.size() << ", " << f << ", " << l2 << ", " << l1 << ", " << linf << ", " << maxCurvature << std::endl;
    }
  }
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
  ImGui::Checkbox("write obj", &writeCurve);
  ImGui::Checkbox("write data", &writeData);
  ImGui::Checkbox("run loop", &runLoop);

  ImGui::SliderFloat("radius", &radius, scene.radius / 100, scene.radius * 100, "%.3f", ImGuiSliderFlags_Logarithmic);
  ImGui::SliderFloat("rmax", &rmax, scene.rmax / 100, scene.rmax * 100, "%.3f", ImGuiSliderFlags_Logarithmic);

  {
    // biharmonic deformation
    using namespace Eigen;

    // Determine boundary conditions
    bc_frac += bc_dir;
    bc_dir *= (bc_frac>=1.0 || bc_frac<=0.0?-1.0:1.0);

    const MatrixXd U_bc_anim = V_bc+bc_frac*(U_bc-V_bc);
    MatrixXd D;
    MatrixXd D_bc = U_bc_anim - V_bc;
    igl::harmonic(V,F,b,D_bc,2,D);
    U = V+D;

    geometry = std::make_unique<VertexPositionGeometry>(*mesh, U);
    geometry->requireVertexTangentBasis();
    geometry->requireFaceTangentBasis();
    geometry->requireVertexNormals();
    geometry->requireVertexPrincipalCurvatureDirections();

    polyscope::getSurfaceMesh("mesh")->updateVertexPositions(U);
  }

  doWork();

  if (ImGui::Button("cut mesh with curve (expermental)")) {
    auto [cutVs, cutFs, _V, _F] = modules::cut_mesh_with_curve(*mesh, *geometry, nodes, segments, segmentSurfacePoints);
    polyscope::registerSurfaceMesh("cut mesh", _V, _F);
    for (int i = 0; i < cutVs.size(); i++) {
      std::cout << "cut mesh " << i << ": " << cutVs[i].rows() << ", " << cutFs[i].rows() << std::endl;
      polyscope::registerSurfaceMesh("cut mesh " + std::to_string(i), cutVs[i], cutFs[i]);
    }
  }

  if (iteration % 100 == 0) {
    // runLoop = false;
  }
}

int main(int argc, char **argv) {

  // Configure the argument parser
  args::ArgumentParser parser("geometry-central & Polyscope example project");
  args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help &h) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a mesh name was given
  if (!inputFilename) {
    std::cerr << "Please specify a mesh file as argument" << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize polyscope
  const char* polyscopeBackend = std::getenv("POLYSCOPE_BACKEND");
  if (polyscopeBackend != nullptr && std::string(polyscopeBackend).size() > 0) {
    polyscope::init(polyscopeBackend);
  } else {
    polyscope::init();
  }

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  // Load scene
  scene = modules::read_scene(args::get(inputFilename));

  radius = scene.radius;
  h = scene.h;
  timestep = scene.timestep;
  rmax = scene.rmax == 0 ? radius * 10 : scene.rmax;

  if (scene.rminNm > 0.0 && scene.nanometersPerUnit > 0.0) {
    radius = scene.rminNm / scene.nanometersPerUnit;
    if (scene.rmax == 0) {
      rmax = radius * 10;
    }
    scene.curvatureBarrierThreshold = scene.nanometersPerUnit / scene.rminNm;

    std::cout << "using physical Rmin: " << scene.rminNm << "nm, "
              << scene.nanometersPerUnit << "nm/unit -> radius = " << radius
              << ", target k_max = " << scene.curvatureBarrierThreshold << std::endl;
  }

  if (scene.meshFileName == "") {
    std::cerr << "Please specify a mesh file in the scene file" << std::endl;
    return EXIT_FAILURE;
  }

  std::tie(mesh, geometry) = readManifoldSurfaceMesh(scene.meshFileName);

  auto psMesh = polyscope::registerSurfaceMesh("mesh", geometry->inputVertexPositions, mesh->getFaceVertexList());

  geometry->requireVertexTangentBasis();
  geometry->requireFaceTangentBasis();
  geometry->requireVertexNormals();
  geometry->requireVertexPrincipalCurvatureDirections();

  smoothedFunction = VertexData<double>(*mesh);

  if (scene.scalarFileName != "") {
    std::ifstream inFile;
    inFile.open(scene.scalarFileName);

    std::cout << "reading scalar file: " << scene.scalarFileName << std::endl;

    if (!inFile) {
      std::cerr << "Could not open file " << scene.scalarFileName << std::endl;
      exit(1);
    }

    int i = 0;
    for (std::string line; std::getline(inFile, line ); ) {
      if (line == "") {
        continue;
      }

      if (i >= mesh->nVertices()) {
        std::cerr << "error!: scalar file has more entries than the mesh" << std::endl;
        exit(1);
      }

      double d = std::stod(line);

      auto v = mesh->vertex(i);
      smoothedFunction[v] = d;

      std::cout << i << ", " << d << std::endl;
      i++;
    }

    inFile.close();
  }
  else
  {
    for (Vertex v : mesh->vertices()) {
      auto coord = geometry->vertexPositions[v];
      // smoothedFunction[v] = std::pow(std::abs(coord.x), 1.5) + std::pow(std::abs(coord.y), 1.5);
      smoothedFunction[v] = 3*std::abs(geometry->vertexGaussianCurvature(v));
    }
  }

  {
    std::string scalarFilename = "scalar.txt";
    std::ofstream scalarfile(scalarFilename);

    for (int i = 0; i < smoothedFunction.size(); i++) {
      auto v = smoothedFunction[i];
      scalarfile << v << std::endl;
    }
  }
  

  VertexData<Vector3> vBasisX(*mesh);
  for (Vertex v : mesh->vertices()) {
    vBasisX[v] = geometry->vertexTangentBasis[v][0];
  }
  psMesh->setVertexTangentBasisX(vBasisX);
  
  vField =
      geometrycentral::surface::computeSmoothestVertexDirectionField(*geometry);

  polyscope::getSurfaceMesh("mesh")->addVertexIntrinsicVectorQuantity("vector field", vField);
  polyscope::getSurfaceMesh("mesh")->addVertexIntrinsicVectorQuantity("principal curvature", geometry->vertexPrincipalCurvatureDirections);
  polyscope::getSurfaceMesh("mesh")->addVertexScalarQuantity("smooth function", smoothedFunction);

  {
    std::string vfieldFilename = "vfield.obj";
    std::ofstream vfieldfile(vfieldFilename);
    
    for (int i = 0; i < vField.size(); i++) {
      auto v = vField[i];
      auto x = geometry->vertexTangentBasis[i][0];
      auto y = geometry->vertexTangentBasis[i][1];

      auto v_cartesian = v.x * x + v.y * y;

      vfieldfile << "v " << v_cartesian.x << " " << v_cartesian.y << " " << v_cartesian.z << std::endl;
    }
  }

  if (scene.curveFileName != "") {
    std::tie(nodes, segments, isFixedNode) = modules::readCurveOnMesh(scene.curveFileName, *mesh, *geometry);
  }

  std::string datafilename = "data.csv";
  data_out = std::ofstream(datafilename);
  std::cout << data_out.is_open() << ", " << datafilename << std::endl;
  data_out << "iteration, time, numNodes, f, descent norm (L2), descent norm (L1), descent norm (L∞), max curvature" << std::endl;

  if (nodes.size() == 0) {
    double rminWorld = radius;
    if (scene.rminNm > 0.0 && scene.nanometersPerUnit > 0.0) {
      rminWorld = scene.rminNm / scene.nanometersPerUnit;
    }

    initialize_parallel_rings(
      *mesh,
      *geometry,
      scene,
      rminWorld,
      nodes,
      segments,
      isFixedNode
    );

    std::cout << "initialized " << scene.initTriangleCount << " parallel ring curve(s)" << std::endl;
  }

  // add boundary loops to fixed nodes
  {
    std::vector<std::array<int, 2>> addedSegments = {};
    for (auto l : mesh->boundaryLoops()) {
      std::vector<int> addedNodes = {};
      for (auto v : l.adjacentVertices()) {
        nodes.emplace_back(SurfacePoint(v));
        isFixedNode.emplace_back(true);
        addedNodes.emplace_back(nodes.size() - 1);
      }

      for (int i = 0; i < addedNodes.size(); i++) {
        std::array<int, 2> newsegment = {addedNodes[i], addedNodes[(i + 1) % addedNodes.size()]};
        segments.emplace_back(newsegment);
        addedSegments.emplace_back(newsegment);
      }
    }

    auto ssp = std::vector<std::vector<SurfacePoint>>(addedSegments.size());
    auto [allNodes, allSegments, _] = formatGeodesicCurve(*mesh, *geometry, nodes, addedSegments, ssp);

    modules::write_curve("boundary_loops.obj", allNodes, allSegments);
  }

  std::tie(segmentSurfacePoints, segmentLengths) = modules::connect_surface_points(*mesh, *geometry, nodes, segments);

  std::tie(nodes, segments, segmentSurfacePoints, segmentLengths, isFixedNode) = modules::remesh_curve_on_surface(*mesh, *geometry, nodes, segments, segmentSurfacePoints, segmentLengths, isFixedNode, h);

  initialNodes = nodes;
  initialSegments = segments;
  initialSegmentSurfacePoints = segmentSurfacePoints;
  initialSegmentLengths = segmentLengths;
  isInitialFixedNode = isFixedNode;

  std::filesystem::create_directories("./objs");
  std::filesystem::create_directories("./radii");

  std::tie(V, F) = modules::gc_to_igl_vf(*mesh, *geometry);

  {
    U=V;
    // S(i) = j: j<0 (vertex i not in handle), j >= 0 (vertex i in handle j)
    Eigen::VectorXi S;
    igl::readDMAT(scene.dmatFilename,S);
    igl::colon<int>(0,V.rows()-1,b);
    b.conservativeResize(std::stable_partition( b.data(), b.data()+b.size(),
    [&S](int i)->bool{return S(i)>=0;})-b.data());

    // Boundary conditions directly on deformed positions
    U_bc.resize(b.size(),V.cols());
    V_bc.resize(b.size(),V.cols());
    for(int bi = 0;bi<b.size();bi++)
    {
      V_bc.row(bi) = V.row(b(bi));
      switch(S(b(bi)))
      {
        case 0:
          // Don't move handle 0
          U_bc.row(bi) = V.row(b(bi));
          break;
        case 1:
          // move handle 1 down
          U_bc.row(bi) = V.row(b(bi)) + Eigen::RowVector3d(-50,-75,0);
          break;
        case 2:
        default:
          // move other handles forward
          U_bc.row(bi) = V.row(b(bi)) + Eigen::RowVector3d(75,50,-75);
          break;
      }
    }
  }

  {
    // for (auto segment : initialSegments) {
    //   std::cout << "segment: " << segment[0] << ", " << segment[1] << std::endl;
    // }

    auto [allNodes, allSegments, cartesianCoords] = formatGeodesicCurve(*mesh, *geometry, initialNodes, initialSegments, initialSegmentSurfacePoints);

    polyscope::registerCurveNetwork("initial curve", allNodes, allSegments)
      // ->setColor({1.0, 1.0, 1.0})
      // ->setRadius(0.003)
    ;

    polyscope::registerCurveNetwork("curve", allNodes, allSegments);

    // polyscope::registerPointCloud("initial nodes", allNodes);

    std::vector<Vector3> fixedNodes = {}, freeNodes = {};
    for (int i = 0; i < initialNodes.size(); i++) {
      if (isInitialFixedNode[i]) {
        fixedNodes.emplace_back(cartesianCoords[i]);
      } else {
        freeNodes.emplace_back(cartesianCoords[i]);
      }
    }

    std::vector<Vector3> normalDirections, tangentX, tangentY;
    for (auto node : nodes) {
      auto [x, y] = modules::get_tangent_basis(*geometry, node);
      tangentX.emplace_back(x / 100);
      tangentY.emplace_back(y / 100);
      normalDirections.emplace_back(cross(x, y) / 100);
    }

    auto [nn, ns, _] = formatVector(cartesianCoords, normalDirections);
    polyscope::registerCurveNetwork("normal directions", nn, ns)->setEnabled(false);

    auto [txn, txs, __] = formatVector(cartesianCoords, tangentX);
    polyscope::registerCurveNetwork("tangent x", txn, txs)->setEnabled(false);

    auto [tyn, tys, ___] = formatVector(cartesianCoords, tangentY);
    polyscope::registerCurveNetwork("tangent y", tyn, tys)->setEnabled(false);

    modules::write_curve("./objs/curve_0.obj", allNodes, allSegments);

    // polyscope::registerPointCloud("free nodes", freeNodes)
    //   ->setPointColor({1.0, 1.0, 1.0})
    // ;
    // polyscope::registerPointCloud("fixed nodes", fixedNodes)
    //   ->setPointColor({.5, .5, .5})
    // ;
  }

  // Give control to the polyscope gui
  if (scene.excecuteOnly) {
    // for debug
    writeData = true;
    writeCurve = true;

    for (int i = 0; i < scene.maxIterations; i++) {
      doWork();
    }
  } else {
    polyscope::show();
  }

  auto finalCartesianCoords = modules::surface_point_to_cartesian(*mesh, *geometry, nodes);
  auto finalCurvature = modules::curve_curvature(finalCartesianCoords, segments);
  std::cout << "final curvature stats -- mean: " << finalCurvature.meanCurvature
            << ", max: " << finalCurvature.maxCurvature
            << ", q75: " << finalCurvature.q75Curvature << std::endl;

  return EXIT_SUCCESS;
}
