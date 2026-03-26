#include "modules/scene_file.h"

#include <algorithm>

namespace modules {
  using namespace std;

  void splitString(const std::string& str, std::vector<string> &cont, char delim) {
      std::stringstream ss(str);
      std::string token;
      while (std::getline(ss, token, delim)) {
        cont.push_back(token);
      }
  }

  std::string getDirectoryFromPath(std::string str) {
    using namespace std;
    vector<string> parts;
    splitString(str, parts, '/');

    int nParts = parts.size();
    if (nParts == 1) return "./";
    
    string path = "";

    for (int i = 0; i < nParts - 1; i++) {
        path = path + parts[i] + "/";
    }

    return path;
  }

  void processLine(SceneObject &scene, std::string directory, std::vector<std::string> &parts) {
    string key = parts[0];

    if (key == "#" || key == "//") {
      return;
    }

    if (key == "curve") {
      scene.curveFileName = directory + parts[1];
    } else if (key == "mesh") {
      scene.meshFileName = directory + parts[1];
    } else if (key == "curve_mesh") {
      scene.meshFileName = directory + parts[1];
    } else if (key == "dmat") {
      scene.dmatFilename = directory + parts[1];
    } else if (key == "scalar") {
      scene.scalarFileName = directory + parts[1];
    } else if (key == "radius") {
      scene.radius = stod(parts[1]);
      scene.h = scene.radius * igl::PI / 20;
      scene.rmax = scene.radius * 5;
    } else if (key == "timestep") {
      scene.timestep = stod(parts[1]);
    } else if (key == "h") {
      scene.h = stod(parts[1]);
    } else if (key == "p") {
      scene.p = stod(parts[1]);
    } else if (key == "q") {
      scene.q = stod(parts[1]);
    } else if (key == "rmax") {
      scene.rmax = stod(parts[1]);
    } else if (key == "field_aligned") {
      scene.w_fieldAlignedness = stod(parts[1]);
    } else if (key == "curxvature_aligned") {
      scene.w_curvatureAlignedness = stod(parts[1]);
    } else if (key == "curvature_constraint") {
      scene.curvatureConstraint = parts[1];
    } else if (key == "curvature_barrier") {
      scene.w_curvatureBarrier = stod(parts[1]);
    } else if (key == "curvature_barrier_threshold") {
      scene.curvatureBarrierThreshold = stod(parts[1]);
    } else if (key == "curvature_barrier_epsilon") {
      scene.curvatureBarrierEpsilon = stod(parts[1]);
    } else if (key == "curvature_barrier_min_length") {
      scene.curvatureBarrierMinLength = stod(parts[1]);
    } else if (key == "rmin_nm") {
      scene.rminNm = stod(parts[1]);
    } else if (key == "nm_per_unit" || key == "nanometers_per_unit") {
      scene.nanometersPerUnit = stod(parts[1]);
    } else if (key == "bilaplacian") {
      scene.w_bilaplacian = stod(parts[1]);
    } else if (key == "init_mode") {
      scene.initMode = parts[1];
    } else if (key == "init_triangles") {
      scene.initTriangleCount = std::max(1, std::stoi(parts[1]));
    } else if (key == "init_stage_iterations") {
      scene.initStageIterations = std::max(1, std::stoi(parts[1]));
    } else if (key == "init_freeze_previous") {
      scene.initFreezePreviousStages = true;
    } else if (key == "init_keep_previous_active") {
      scene.initFreezePreviousStages = false;
    } else if (key == "init_ring_vertices") {
      scene.initRingVertexCount = std::max(3, std::stoi(parts[1]));
    } else if (key == "init_stack_axis_mode") {
      scene.initStackAxisMode = parts[1];
    } else if (key == "init_stack_axis") {
      scene.initStackAxisMode = parts[1];
      if (scene.initStackAxisMode == "highest") {
        scene.initStackAxisMode = "long";
      } else if (scene.initStackAxisMode == "lowest") {
        scene.initStackAxisMode = "short";
      } else if (scene.initStackAxisMode == "middle" || scene.initStackAxisMode == "median") {
        scene.initStackAxisMode = "mid";
      }
    } else if (key == "init_stack_rotate_deg") {
      scene.initStackRotateDeg = stod(parts[1]);
    } else if (key == "init_stack_tilt_deg") {
      scene.initStackTiltDeg = stod(parts[1]);
    } else if (key == "init_stack_ref") {
      scene.initStackRefX = stod(parts[1]);
      scene.initStackRefY = stod(parts[2]);
      scene.initStackRefZ = stod(parts[3]);
    } else if (key == "init_ring_radius") {
      scene.initRingRadius = stod(parts[1]);
    } else if (key == "init_ring_spacing") {
      scene.initRingSpacing = stod(parts[1]);
    } else if (key == "init_ring_radius_nm") {
      scene.initRingRadiusNm = stod(parts[1]);
    } else if (key == "init_ring_spacing_nm") {
      scene.initRingSpacingNm = stod(parts[1]);
    } else if (key == "init_stage_clearance") {
      scene.initStageClearance = stod(parts[1]);
    } else if (key == "init_stage_clearance_nm") {
      scene.initStageClearanceNm = stod(parts[1]);
    } else if (key == "ring_insertion") {
      scene.ringInsertionEnabled = true;
    } else if (key == "ring_insertion_every") {
      scene.ringInsertionEvery = std::max(1, std::stoi(parts[1]));
    } else if (key == "ring_insertion_after_iter") {
      scene.ringInsertionAfterIter = std::max(0, std::stoi(parts[1]));
    } else if (key == "ring_insertion_max") {
      scene.ringInsertionMax = std::max(0, std::stoi(parts[1]));
    } else if (key == "ring_insertion_spacing") {
      scene.ringInsertionSpacing = stod(parts[1]);
    } else if (key == "ring_insertion_spacing_nm") {
      scene.ringInsertionSpacingNm = stod(parts[1]);
    } else if (key == "max_iterations") {
      scene.maxIterations = std::max(1, std::stoi(parts[1]));
    } else if (key == "varying_alpha") {
      scene.varyingAlpha = true;
    } else if (key == "geodesic_medial_axis") {
      scene.useGeodesicMedialAxis = true;
    } else if (key == "excecute_only") {
      scene.excecuteOnly = true;
    }
  }

  SceneObject read_scene(std::string filename) {
    string directory = getDirectoryFromPath(filename);

    ifstream inFile;
    inFile.open(filename);

    if (!inFile) {
      cerr << "Could not open file " << filename << endl;
      exit(1);
    }

    SceneObject scene;

    std::vector<std::string> parts;
    for (std::string line; std::getline(inFile, line ); ) {
      if (line == "" || line == "\n") continue;
      parts.clear();
      splitString(line, parts, ' ');
      processLine(scene, directory, parts);
    }

    inFile.close();
    return scene;
  }
}
