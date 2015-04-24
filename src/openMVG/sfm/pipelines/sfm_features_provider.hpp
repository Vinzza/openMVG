
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_FEATURES_PROVIDER_HPP
#define OPENMVG_SFM_FEATURES_PROVIDER_HPP

#include <openMVG/types.hpp>
#include <openMVG/sfm/sfm_data.hpp>
#include <openMVG/features/features.hpp>

namespace openMVG{

typedef std::vector<PointFeature> PointFeatures;

/// Abstract PointFeature provider (read some feature and store them as PointFeature).
/// Allow to load and return the features related to a view
struct Features_Provider
{
  /// PointFeature array per ViewId of the considered SfM_Data container
  Hash_Map<IndexT, PointFeatures> feats_per_view;

  // Load features related to a provide SfM_Data View container
  virtual bool load(const SfM_Data & sfm_data, const std::string & feat_directory)
  {
    // Read for each view the corresponding features and store them as PointFeatures
    for (Views::const_iterator iter = sfm_data.getViews().begin();
      iter != sfm_data.getViews().end(); ++iter)
    {
      const std::string sImageName = stlplus::create_filespec(sfm_data.s_root_path, iter->second.get()->s_Img_path);
      const std::string basename = stlplus::basename_part(sImageName);
      std::vector<SIOPointFeature> vec_feats;
      if (!loadFeatsFromFile( stlplus::create_filespec(feat_directory, basename, ".feat"), vec_feats)) {
        std::cerr << "Invalid feature files for the view: "<< sImageName << std::endl;
        return false;
      }
      // convert loaded SIOPointFeatures to PointFeature
      feats_per_view[iter->second.get()->id_view] = std::vector<PointFeature>(vec_feats.begin(), vec_feats.end());
    }
    return true;
  }

  /// Return the PointFeatures belonging to the View, if the view does not exist
  ///  it return an empty PointFeature array.
  const PointFeatures & getFeatures(const IndexT & id_view) const
  {
    // Have a empty feature set in order to deal with non existing view_id
    static const PointFeatures emptyFeats = PointFeatures();

    Hash_Map<IndexT, PointFeatures>::const_iterator it = feats_per_view.find(id_view);
    if (it != feats_per_view.end())
      return it->second;
    else
      return emptyFeats;
  }
}; // Features_Provider

} // namespace openMVG

#endif // OPENMVG_SFM_FEATURES_PROVIDER_HPP
