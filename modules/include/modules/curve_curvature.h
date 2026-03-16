#pragma once

#include <vector>

#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace modules {

struct CurveCurvatureData {
  std::vector<double> nodeCurvature;
  std::vector<double> tangentDot;
  std::vector<double> localLength;
  std::vector<bool> validNode;

  double meanCurvature = 0.0;
  double maxCurvature = 0.0;
  double q75Curvature = 0.0;
};

CurveCurvatureData curve_curvature(
  const std::vector<Vector3>& nodePositions,
  // Curve edges as node-index pairs. These define curve connectivity/order.
  // We keep the name "segments" to match the rest of the codebase.
  const std::vector<std::array<int, 2>>& segments
);

}
