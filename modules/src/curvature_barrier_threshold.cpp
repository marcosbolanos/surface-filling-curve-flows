#include "modules/curvature_barrier_threshold.h"

#include <algorithm>
#include <cmath>

namespace modules {

CurvatureBarrierThreshold curvature_barrier_threshold(
  const std::vector<Vector3>& nodePositions,
  const std::vector<std::array<int, 2>>& segments,
  const std::vector<double>& segmentLengths,
  double minLength
) {
  CurvatureBarrierThreshold out;

  std::vector<std::vector<int>> nodeNeighbors(nodePositions.size());
  std::vector<std::vector<int>> nodeIncidentSegments(nodePositions.size());
  for (int i = 0; i < static_cast<int>(segments.size()); i++) {
    int v0 = segments[i][0];
    int v1 = segments[i][1];

    if (v0 < 0 || v1 < 0 || v0 >= static_cast<int>(nodePositions.size()) || v1 >= static_cast<int>(nodePositions.size())) {
      continue;
    }

    nodeNeighbors[v0].push_back(v1);
    nodeNeighbors[v1].push_back(v0);
    nodeIncidentSegments[v0].push_back(i);
    nodeIncidentSegments[v1].push_back(i);
  }

  for (int i = 0; i < static_cast<int>(nodePositions.size()); i++) {
    if (nodeNeighbors[i].size() != 2) {
      continue;
    }

    int prev = nodeNeighbors[i][0];
    int next = nodeNeighbors[i][1];

    Vector3 ePrev = nodePositions[i] - nodePositions[prev];
    Vector3 eNext = nodePositions[next] - nodePositions[i];

    double lPrev = norm(ePrev);
    double lNext = norm(eNext);
    if (lPrev <= 1e-12 || lNext <= 1e-12) {
      continue;
    }

    Vector3 tPrev = ePrev / lPrev;
    Vector3 tNext = eNext / lNext;

    double d = std::clamp(dot(tPrev, tNext), -1.0, 1.0);
    double theta = std::acos(d);

    double lBar = 0.5 * (lPrev + lNext);
    if (segmentLengths.size() == segments.size() && nodeIncidentSegments[i].size() == 2) {
      int s0 = nodeIncidentSegments[i][0];
      int s1 = nodeIncidentSegments[i][1];
      lBar = 0.5 * (segmentLengths[s0] + segmentLengths[s1]);
    }
    lBar = std::max(lBar, minLength);
    if (lBar <= 0.0) {
      continue;
    }

    out.minFeasibleKmax = std::max(out.minFeasibleKmax, theta / lBar);
    out.validNodeCount++;
  }

  return out;
}

}
