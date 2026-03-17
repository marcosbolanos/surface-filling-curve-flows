#pragma once

#include <array>
#include <vector>

#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace modules {

struct CurvatureBarrierThreshold {
  double minFeasibleKmax = 0.0;
  int validNodeCount = 0;
};

CurvatureBarrierThreshold curvature_barrier_threshold(
  const std::vector<Vector3>& nodePositions,
  const std::vector<std::array<int, 2>>& segments,
  const std::vector<double>& segmentLengths,
  double minLength = 0.0
);

}
