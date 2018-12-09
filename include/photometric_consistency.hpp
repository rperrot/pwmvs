#ifndef _PWMVS_PHOTOMETRIC_CONSISTENCY_HPP_
#define _PWMVS_PHOTOMETRIC_CONSISTENCY_HPP_

#include <Eigen/StdVector>

#include "geometry.hpp"
#include "photometric_consistency_cache.hpp"
#include "types.hpp"

#include <openMVG/image/sample.hpp>

struct PhotometricConsistency
{
public: // TODO: for testing purposes
  struct AbstractWindow
  {
    virtual FloatT operator()( int row, int col ) const = 0;
  };

  struct WindowProduct : public AbstractWindow
  {
    WindowProduct( const AbstractWindow& a, const AbstractWindow& b )
        : AbstractWindow(), a( a ), b( b )
    {
      // nothing to do
    }

    virtual FloatT operator()( int window_row, int window_col ) const override
    {
      return a( window_row, window_col ) * b( window_row, window_col );
    }

  private:
    const AbstractWindow& a;
    const AbstractWindow& b;
  };

  struct Window : public AbstractWindow
  {
    Window( const Image<FloatT>& image, const Vector2i& ref_center )
        : image( image ), ref_center( ref_center )
    {
      // nothing to do
    }

    WindowProduct operator*( const Window& other ) const
    {
      return WindowProduct( *this, other );
    }

  protected:
    const Image<FloatT>& image;
    const Vector2i       ref_center;
  };

  struct RefWindow : public Window
  {
    RefWindow( const Image<FloatT>& ref_image, const Vector2i& ref_center )
        : Window( ref_image, ref_center )
    {
      // nothing to do
    }

    virtual FloatT operator()( int window_row, int window_col ) const override
    {
      int row = ref_center.y() + window_row;
      int col = ref_center.x() + window_col;

      if ( !image.Contains( row, col ) )
      {
        row = std::max( image.Height() - 1, std::max( 0, row ) );
        col = std::max( image.Width() - 1, std::max( 0, col ) );
        return image( row, col );
      }

      //return image(row, col);

      return sampler( image, convertToOpenMVG( convertToFloat( row ) ), convertToOpenMVG( convertToFloat( col ) ) );
    }

    openMVG::image::Sampler2d<openMVG::image::SamplerLinear> sampler;
  };

  struct SrcWindow : public Window
  {
    SrcWindow( const Image<FloatT>& src_image, const Vector2i& ref_center, const Matrix3& H )
        : Window( src_image, ref_center ), H( H ), src_center( transfer( 0, 0 ) )
    {
      // nothing to do
    }

    virtual FloatT operator()( int window_row, int window_col ) const override
    {
      Vector2 x = transfer( window_row, window_col );

      if ( !image.Contains( x ) )
      {
        x( 0 ) = std::min( FloatT( image.Height() ), std::max( FloatT( 0 ), x( 0 ) ) );
        x( 1 ) = std::min( FloatT( image.Width() ), std::max( FloatT( 0 ), x( 1 ) ) );
      }

      //return image(convertToInt(x));

      return sampler( image, convertToOpenMVG( x.y() ), convertToOpenMVG( x.x() ) );
    }

  private:
    Vector2 transfer( int window_row, int window_col ) const
    {
      return ( H * convertToFloat( Vector2i( window_col + ref_center.x(), window_row + ref_center.y() ) ).homogeneous() ).hnormalized();
    }

  private:
    const Matrix3& H;
    const Vector2  src_center;

    openMVG::image::Sampler2d<openMVG::image::SamplerLinear> sampler;
  };

public:
  /// NOTE: this implementation is not optimized at all. It serves as reference implementation. Many calculations repeat many times and could be cached.
  PhotometricConsistency( const RefView& ref, const int window_size, const FloatT color_dispersion, FloatT spatial_dispersion )
      : window_size( window_size ), color_dispersion( 1 / ( 2 * color_dispersion * color_dispersion ) ), spatial_dispersion( 1 / ( 2 * spatial_dispersion * spatial_dispersion ) ), ref( ref )
  {
    assert( window_size > 0 );
    assert( color_dispersion > 0 );
    assert( spatial_dispersion > 0 );
  }

public: // TODO: for testing purposes
  FloatT operator()( const SrcView& src, const Vector2i& x )
  {
    return ( *this )( src, x, ref.normal( x ), ref.depth( x ) );
  }

  FloatT operator()( const SrcView& src, const Vector2i& x, const Normal& normal, FloatT depth )
  {
    if ( !src.isVisible( src.project( ref.unproject( convertToFloat( x ), depth ) ) ) )
      return -1;

    const Matrix3 H = ComputeHomography( ref, src, normal, depth, convertToFloat( x ) );

    RefWindow ref_window( ref.image, x );
    SrcWindow src_window( src.image, x, H );

    FloatT cov_src_src = covariance( ref_window, src_window, src_window );
    if ( cov_src_src < 1e-5 )
      return -1;
    FloatT cov_ref_ref = covariance( ref_window, ref_window, ref_window );
    if ( cov_src_src < 1e-5 )
      return -1;

    FloatT cov_ref_src = covariance( ref_window, ref_window, src_window );

    FloatT ncc = cov_ref_src / std::sqrt( cov_ref_ref * cov_src_src );
    ncc        = std::max( static_cast<FloatT>( -1 ), std::min( static_cast<FloatT>( +1 ), ncc ) );
    return ncc;
  }

  FloatT covariance( const RefWindow& ref, const Window& a, const Window& b )
  {
    return weightedAverage( ref, a * b ) - weightedAverage( ref, a ) * weightedAverage( ref, b );
  }

  FloatT weightedAverage( const RefWindow& ref_window, const AbstractWindow& window )
  {
    FloatT ref_center_color = ref_window( 0, 0 );

    FloatT weighted_color_sum = 0;
    FloatT weight_sum         = 0;
    for ( int row = -window_size; row <= +window_size; row++ )
    {
      for ( int col = -window_size; col <= +window_size; col++ )
      {
        FloatT color = window( row, col );
        FloatT w     = weight( Vector2( col, row ), ref_window( row, col ) - ref_center_color );

        weighted_color_sum += w * color;
        weight_sum += w;
      }
    }

    return weighted_color_sum / weight_sum;
  }

  FloatT weight( const Vector2& spatial_diff, const FloatT color_diff )
  {
    return std::exp( -( color_diff * color_diff ) * color_dispersion - spatial_diff.squaredNorm() * spatial_dispersion );
  }

private:
  int    window_size;
  FloatT color_dispersion;
  FloatT spatial_dispersion;

  const RefView& ref;
};

template <class Sampler = openMVG::image::SamplerLinear>
struct PhotometricConsistencyOptimized
{
  PhotometricConsistencyOptimized( const RefView& ref, const int window_size, const FloatT color_dispersion, const FloatT spatial_dispersion )
      : window_size( window_size ),
        color_dispersion( 1 / ( 2 * color_dispersion * color_dispersion ) ),
        spatial_dispersion( 1 / ( 2 * spatial_dispersion * spatial_dispersion ) ),
        ref( ref ),
        ref_weighted_color_sum( ref.image.Width(), ref.image.Height(), true, 0 ),
        ref_weighted_color_sum_sq( ref.image.Width(), ref.image.Height(), true, 0 )
  {
    assert( window_size > 0 );
    assert( color_dispersion > 0 );
    assert( spatial_dispersion > 0 );

    precalculate();
  }

  PhotometricConsistencyOptimized( const PhotometricConsistencyOptimized<Sampler>& src ) = default;
  PhotometricConsistencyOptimized( PhotometricConsistencyOptimized<Sampler>&& src )      = default;

  PhotometricConsistencyOptimized<Sampler>& operator=( const PhotometricConsistencyOptimized<Sampler>& src ) = default;
  PhotometricConsistencyOptimized<Sampler>& operator=( PhotometricConsistencyOptimized<Sampler>&& src ) = default;

public:
  FloatT operator()( const SrcView& src, const Vector2i& x )
  {
    return ( *this )( src, x, ref.normal( x ), ref.depth( x ) );
  }

  // Version with cache
  FloatT operator()( const SrcView& src, const Vector2i& x, const Normal& normal, const FloatT depth, const PhotometricConistencyCache& cache ) const
  {
    const Vector3 X = ref.unproject( convertToFloat( x ), depth );

    if ( !src.isVisible( src.project( X ) ) )
      return -1;

    const FloatT center_color = Sample( ref.image, x );

    FloatT weights_sum                = 0;
    FloatT weighted_src_color_sum     = 0;
    FloatT weighted_src_ref_color_sum = 0;
    FloatT weighted_src_color_sum_sq  = 0;

    const FloatT weighted_ref_color_sum    = ref_weighted_color_sum( x );
    const FloatT weighted_ref_color_sum_sq = ref_weighted_color_sum_sq( x );

    const FloatT d = X.dot( normal );

    size_t cur_index = 0;
    for ( int window_row = -window_size; window_row <= window_size; window_row++ )
    {
      for ( int window_col = -window_size; window_col <= window_size; window_col++ )
      {
        assert( cur_index < cache._cache_ray.size() );
        const Normal  ray   = cache._cache_ray[ cur_index ]; // ref.ray( convertToFloat( ( x + x_ ).eval() ) );
        const Vector2 x_src = src.project( ( d / ray.dot( normal ) ) * ray );

        // TODO: since this function is called lot's of times, it may be interesting to cache reference image values
        const FloatT ref_color = cache._cache_ref_color[ cur_index ]; // Sample( ref.image, ( x + x_ ).eval() );
        const FloatT src_color = Sample( src.image, x_src, false );

        const FloatT w = weight( window_row, window_col, center_color - ref_color );

        const FloatT wSrcColor = w * src_color;

        weighted_src_color_sum += wSrcColor;
        weighted_src_color_sum_sq += wSrcColor * src_color;
        weighted_src_ref_color_sum += wSrcColor * ref_color;
        weights_sum += w;

        ++cur_index;
      }
    }

    const FloatT inv_sum = 1 / weights_sum;

    weighted_src_color_sum *= inv_sum;
    weighted_src_color_sum_sq *= inv_sum;
    weighted_src_ref_color_sum *= inv_sum;

    const FloatT ref_color_var       = weighted_ref_color_sum_sq - ( weighted_ref_color_sum * weighted_ref_color_sum );
    const FloatT src_color_var       = weighted_src_color_sum_sq - ( weighted_src_color_sum * weighted_src_color_sum );
    const FloatT ref_src_color_covar = weighted_src_ref_color_sum - ( weighted_ref_color_sum * weighted_src_color_sum );

    if ( ref_color_var < 1e-6 )
      return -1;
    if ( src_color_var < 1e-6 )
      return -1;

    FloatT ncc = ref_src_color_covar / std::sqrt( ref_color_var * src_color_var );
    ncc        = std::max( static_cast<FloatT>( -1 ), std::min( static_cast<FloatT>( +1 ), ncc ) );
    return ncc;
  }

  // Version without cache
  FloatT operator()( const SrcView& src, const Vector2i& x, const Normal& normal, const FloatT depth ) const
  {
    const Vector3 X = ref.unproject( convertToFloat( x ), depth );

    if ( !src.isVisible( src.project( X ) ) )
      return -1;

    const FloatT center_color = Sample( ref.image, x );

    FloatT weights_sum                = 0;
    FloatT weighted_src_color_sum     = 0;
    FloatT weighted_src_ref_color_sum = 0;
    FloatT weighted_src_color_sum_sq  = 0;

    const FloatT weighted_ref_color_sum    = ref_weighted_color_sum( x );
    const FloatT weighted_ref_color_sum_sq = ref_weighted_color_sum_sq( x );

    const FloatT d = X.dot( normal );

    for ( int window_row = -window_size; window_row <= window_size; window_row++ )
    {
      for ( int window_col = -window_size; window_col <= window_size; window_col++ )
      {
        const Vector2i x_( window_col, window_row );

        const Normal  ray   = ref.ray( convertToFloat( ( x + x_ ).eval() ) );
        const Vector2 x_src = src.project( ( d / ray.dot( normal ) ) * ray );

        const FloatT ref_color = Sample( ref.image, ( x + x_ ).eval() );
        const FloatT src_color = Sample( src.image, x_src, false );

        const FloatT w = weight( window_row, window_col, center_color - ref_color );

        const FloatT wSrcColor = w * src_color;

        weighted_src_color_sum += wSrcColor;
        weighted_src_color_sum_sq += wSrcColor * src_color;
        weighted_src_ref_color_sum += wSrcColor * ref_color;
        weights_sum += w;
      }
    }

    const FloatT inv_sum = 1 / weights_sum;

    weighted_src_color_sum *= inv_sum;
    weighted_src_color_sum_sq *= inv_sum;
    weighted_src_ref_color_sum *= inv_sum;

    const FloatT ref_color_var       = weighted_ref_color_sum_sq - ( weighted_ref_color_sum * weighted_ref_color_sum );
    const FloatT src_color_var       = weighted_src_color_sum_sq - ( weighted_src_color_sum * weighted_src_color_sum );
    const FloatT ref_src_color_covar = weighted_src_ref_color_sum - ( weighted_ref_color_sum * weighted_src_color_sum );

    if ( ref_color_var < 1e-6 )
      return -1;
    if ( src_color_var < 1e-6 )
      return -1;

    FloatT ncc = ref_src_color_covar / std::sqrt( ref_color_var * src_color_var );
    ncc        = std::max( static_cast<FloatT>( -1 ), std::min( static_cast<FloatT>( +1 ), ncc ) );
    return ncc;
  }

  /**
    * @brief Compute and store in cache values for a given pixel in reference image 
    * 
    * @param x 
    */
  void computeCache( const Vector2i& x, PhotometricConistencyCache& cache ) const
  {
    const size_t cache_size = ( 2 * window_size + 1 ) * ( 2 * window_size + 1 );
    if ( cache_size > cache._cache_ref_color.size() )
    {
      cache._cache_ref_color.resize( cache_size );
      cache._cache_ray.resize( cache_size );
    }

    size_t index = 0;
    for ( int window_row = -window_size; window_row <= window_size; window_row++ )
    {
      for ( int window_col = -window_size; window_col <= window_size; window_col++ )
      {
        const Vector2i x_( window_col, window_row );
        assert( index < cache_size );

        cache._cache_ray[ index ]       = ref.ray( convertToFloat( ( x + x_ ).eval() ) );
        cache._cache_ref_color[ index ] = Sample( ref.image, ( x + x_ ).eval() );
        ++index;
      }
    }
  }

private:
  void precalculate()
  {
    // Precompute weights
    _precomputed_spatial_weight.resize( ( 2 * window_size + 1 ) * ( 2 * window_size + 1 ) );
    for ( int window_row = -window_size; window_row <= window_size; window_row++ )
    {
      for ( int window_col = -window_size; window_col <= window_size; window_col++ )
      {
        const int id_y  = window_row + window_size;
        const int id_x  = window_col + window_size;
        const int index = id_x + ( 2 * window_size + 1 ) * id_y;
        assert( index < _precomputed_spatial_weight.size() );

        _precomputed_spatial_weight[ index ] = std::exp( -( window_row * window_row + window_col * window_col ) * spatial_dispersion );
      }
    }
    _precomputed_color_weight.resize( 512 );
    for ( int color_diff = -255; color_diff <= 255; ++color_diff )
    {
      const FloatT floatDiff = 1.f / 255.f;
      const size_t cur_index = 255 + color_diff;
      assert( cur_index < _precomputed_color_weight.size() );
      _precomputed_color_weight[ 255 + color_diff ] = std::exp( -( floatDiff * floatDiff ) * color_dispersion );
    }

    for ( int row = window_size; row < ( ref.image.Height() - window_size ); row++ )
    {
#pragma omp parallel for
      for ( int col = window_size; col < ( ref.image.Width() - window_size ); col++ )
      {
        const FloatT center_color = Sample( ref.image, Vector2i( col, row ) );

        FloatT weights_sum           = 0;
        FloatT weighted_color_sum    = 0;
        FloatT weighted_color_sum_sq = 0;

        for ( int window_row = -window_size; window_row <= window_size; window_row++ )
        {
          for ( int window_col = -window_size; window_col <= window_size; window_col++ )
          {
            const FloatT color = Sample( ref.image, Vector2i( col + window_col, row + window_row ) );
            const FloatT w     = weight( window_row, window_col, center_color - color );

            const FloatT wColor = w * color;
            weighted_color_sum += wColor;
            weighted_color_sum_sq += wColor * color;
            weights_sum += w;
          }
        }

        ref_weighted_color_sum( row, col )    = weighted_color_sum / weights_sum;
        ref_weighted_color_sum_sq( row, col ) = weighted_color_sum_sq / weights_sum;
      }
    }
  }

  FloatT weight( int row_diff, int col_diff, const FloatT color_diff ) const
  {
    // uses precomputed weigth
    // exp( a + b ) = exp( a ) * exp( b )
    // exp( a ) could be precomputed (row_diff and col_diff are limited to 2 * window size)
    // exp( b ) also color_diff is limited since it corresponds to image pixel values (0;255)
    // for history, the computed value is:
    // return std::exp( -( color_diff * color_diff ) * color_dispersion - ( row_diff * row_diff + col_diff * col_diff ) * spatial_dispersion );
    const int id_y  = row_diff + window_size;
    const int id_x  = col_diff + window_size;
    const int index = id_x + id_y * ( window_size * 2 + 1 ); // TODO : precompute win*2+1

    const float spatial_weight = _precomputed_spatial_weight[ index ];

    const int   i_color_diff = std::min( std::max( (int)( color_diff * 255 ) + 255, 0 ), 510 );
    const float color_weight = _precomputed_color_weight[ i_color_diff ];

    return spatial_weight * color_weight;
  }

  FloatT Sample( const Image<FloatT>& image, const Vector2& x, bool clamp_to_edge ) const
  {
    if ( clamp_to_edge && !image.Contains( x ) )
    {
      Vector2i x_;
      x_( 0 ) = std::min( image.Width() - 1, std::max( 0, convertToInt( x( 0 ) ) ) );
      x_( 1 ) = std::min( image.Height() - 1, std::max( 0, convertToInt( x( 1 ) ) ) );

      return Sample( image, x_ );
    }

    /*
    const int pix_x = std::min( image.Width() - 1, std::max( 0, convertToInt( x( 0 ) ) ) );
    const int pix_y = std::min( image.Height() - 1, std::max( 0, convertToInt( x( 1 ) ) ) );
    */
    return sampler( image, convertToOpenMVG( x.y() ), convertToOpenMVG( x.x() ) );
  }

  FloatT Sample( const Image<FloatT>& image, const Vector2i& x ) const
  {
    return image( x );
  }

private:
  int    window_size;
  FloatT color_dispersion;
  FloatT spatial_dispersion;

  const RefView& ref;

  Image<FloatT> ref_weighted_color_sum;
  Image<FloatT> ref_weighted_color_sum_sq;

  std::vector<FloatT, Eigen::aligned_allocator<FloatT>> _precomputed_spatial_weight;
  std::vector<FloatT, Eigen::aligned_allocator<FloatT>> _precomputed_color_weight;

  openMVG::image::Sampler2d<Sampler> sampler;
};

#endif