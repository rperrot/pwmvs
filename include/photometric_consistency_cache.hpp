#ifndef _PHOTOMETRIC_CONSISTENCY_CACHE_HPP
#define _PHOTOMETRIC_CONSISTENCY_CACHE_HPP

#include "types.hpp"

/**
 * @brief Per thread photometric consistency cache 
 * 
 */
struct PhotometricConistencyCache
{
  std::vector<Normal, Eigen::aligned_allocator<Normal>> _cache_ray;
  std::vector<FloatT, Eigen::aligned_allocator<FloatT>> _cache_ref_color;
};

#endif