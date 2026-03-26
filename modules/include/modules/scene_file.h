#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <igl/PI.h>

namespace modules {
  enum CurveType {
    Closed,
    Open
  };

  struct SceneObject {
    std::string curveFileName = "";
    std::string fixedCurveFileName = "";
    std::string meshFileName = "";
    std::string dmatFilename = "";
    std::string scalarFileName = "";
    CurveType curveType = CurveType::Closed;
    std::vector<int> startIds;
    std::vector<int> endIds;
    double velocity = 1.0;
    double radius = 0.1;
    double alpha = 1;
    double timestep = 1.;
    double h = igl::PI * radius / 20;
    double p = 2;
    double q = 2;
    double rmax = radius * 5;
    double w_fieldAlignedness = 0;
    double w_curvatureAlignedness = 0;
    std::string curvatureConstraint = "";
    double w_curvatureBarrier = 0;
    double curvatureBarrierThreshold = 30;
    double curvatureBarrierEpsilon = 1e-6;
    double curvatureBarrierMinLength = 0;
    double rminNm = 0;
    double nanometersPerUnit = 1;
    double w_bilaplacian = 0;
    std::string initMode = "rings";
    int initTriangleCount = 1;
    int initStageIterations = 100;
    bool initFreezePreviousStages = true;
    int initRingVertexCount = 24;
    std::string initStackAxisMode = "long";
    double initStackRotateDeg = 0;
    double initStackTiltDeg = 0;
    double initStackRefX = 0;
    double initStackRefY = 0;
    double initStackRefZ = 1;
    double initRingRadius = 0;
    double initRingSpacing = 0;
    double initRingRadiusNm = 0;
    double initRingSpacingNm = 0;
    double initStageClearance = 0;
    double initStageClearanceNm = 0;
    bool ringInsertionEnabled = false;
    int ringInsertionEvery = 10;
    int ringInsertionAfterIter = 10;
    int ringInsertionMax = 0;
    double ringInsertionSpacing = 0;
    double ringInsertionSpacingNm = 0;
    bool varyingAlpha = false;
    bool useGeodesicMedialAxis = false;
    int maxIterations = 100;
    bool excecuteOnly = false;
  };

  void splitString(const std::string &str, std::vector<std::string> &cont, char delim = ' ');

  SceneObject read_scene(std::string filename);
}
