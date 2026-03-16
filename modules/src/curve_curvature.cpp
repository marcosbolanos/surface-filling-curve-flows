#include "modules/curve_curvature.h"

#include <algorithm>
#include <cmath>

namespace modules {

CurveCurvatureData curve_curvature(
  const std::vector<Vector3>& nodePositions,
  const std::vector<std::array<int, 2>>& segments
) {
  CurveCurvatureData out;

  out.nodeCurvature = std::vector<double>(nodePositions.size(), 0.0);
  out.tangentDot = std::vector<double>(nodePositions.size(), 1.0);
  out.localLength = std::vector<double>(nodePositions.size(), 0.0);
  out.validNode = std::vector<bool>(nodePositions.size(), false);

  // Recover local curve connectivity from segment edges.
  // For curvature at node i, we use its two incident neighbors (degree-2 node).
  std::vector<std::vector<int>> nodeNeighbors(nodePositions.size());
  for (int i = 0; i < segments.size(); i++) {
    int v0 = segments[i][0];
    int v1 = segments[i][1];

    if (v0 < 0 || v1 < 0 || v0 >= nodePositions.size() || v1 >= nodePositions.size()) {
      continue;
    }

    nodeNeighbors[v0].push_back(v1);
    nodeNeighbors[v1].push_back(v0);
  }

  std::vector<double> validCurvatures;
  for (int i = 0; i < nodePositions.size(); i++) {
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

    double dotVal = std::clamp(dot(tPrev, tNext), -1.0, 1.0);
    double l = 0.5 * (lPrev + lNext);
    double curvature = norm(tNext - tPrev) / l;

    out.nodeCurvature[i] = curvature;
    out.tangentDot[i] = dotVal;
    out.localLength[i] = l;
    out.validNode[i] = true;
    validCurvatures.push_back(curvature);
  }

  if (validCurvatures.empty()) {
    return out;
  }

  double sum = 0.0;
  for (double k : validCurvatures) {
    sum += k;
  }

  out.meanCurvature = sum / validCurvatures.size();
  out.maxCurvature = *std::max_element(validCurvatures.begin(), validCurvatures.end());

  std::sort(validCurvatures.begin(), validCurvatures.end());
  int q75Id = static_cast<int>(std::ceil(0.75 * validCurvatures.size())) - 1;
  q75Id = std::clamp(q75Id, 0, static_cast<int>(validCurvatures.size()) - 1);
  out.q75Curvature = validCurvatures[q75Id];

  return out;
}

}
