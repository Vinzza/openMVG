// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP
#define OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP

#include "openMVG/sfm/pipelines/sfm_engine.hpp"

#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"

namespace openMVG{

/// Global SfM Pipeline Reconstruction Engine.
/// - Method: Global Fusion of Relative Motions.
class GlobalSfMReconstructionEngine_RelativeMotions : public ReconstructionEngine
{
public:

  GlobalSfMReconstructionEngine_RelativeMotions(
    const SfM_Data & sfm_data,
    const std::string & soutDirectory,
    const std::string & loggingFile = "");

  ~GlobalSfMReconstructionEngine_RelativeMotions();

  void SetFeaturesProvider(Features_Provider * provider);
  void SetMatchesProvider(Matches_Provider * provider);

  void SetRotationAveragingMethod(globalSfM::ERotationAveragingMethod eRotationAveragingMethod);
  void SetTranslationAveragingMethod(globalSfM::ETranslationAveragingMethod _eTranslationAveragingMethod);

  virtual bool Process();

protected:
  /// Compute from relative rotations the global rotations of the camera poses
  bool Compute_Global_Rotations();

  /// Compute/refine relative translations and compute global translations
  bool Compute_Global_Translations();

  /// Compute the initial structure of the scene
  bool Compute_Initial_Structure();

  // Adjust the scene (& remove outliers)
  bool Adjust();

private:
  /// Compute relative rotations
  void Compute_Relative_Rotations(RelativeInfo_Map & vec_relatives);

  //----
  //-- Data
  //----

  // HTML logger
  std::shared_ptr<htmlDocument::htmlDocumentStream> _htmlDocStream;
  std::string _sLoggingFile;
  bool Log_Display_graph( const std::string graph_name, const graphUtils::indexedGraph & putativeGraph );
  bool Log_Display_graph( const std::string graph_name, const PairWiseMatches & map_matches );
  bool Log_Display_graph( const std::string graph_name, const RelativeInfo_Map & relatives_Rt );

  // Parameter
  globalSfM::ERotationAveragingMethod _eRotationAveragingMethod;
  globalSfM::ETranslationAveragingMethod _eTranslationAveragingMethod;

  //-- Data provider
  Features_Provider  * _features_provider;
  Matches_Provider  * _matches_provider;

  std::shared_ptr<Features_Provider> _normalized_features_provider;

  // Temporary data
  RelativeInfo_Map _relatives_Rt; // The relative [R|t] extracted from essential matrices
  Hash_Map<IndexT, Mat3> _map_globalR; // The global rotations
  matching::PairWiseMatches _tripletWise_matches; // Matches valided by triplet at the translation averaging step

};

} // namespace openMVG

#endif // OPENMVG_SFM_GLOBAL_ENGINE_RELATIVE_MOTIONS_HPP
