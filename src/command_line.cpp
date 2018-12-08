#include "command_line.hpp"

void usage( const std::string& pgm )
{
  std::cout << std::endl;
  std::cout << "Usage: " << pgm << " [options] -i /path/to/sfm_data.bin -o /path/to/output/folder"
            << std::endl
            << std::endl;
  std::cout << "Options: " << std::endl;
  std::cout << "--------------------------------------------------------------------------------------" << std::endl;
  std::cout << " General parameters" << std::endl;
  std::cout << "  -i        Input openMVG output (ie: sfm_data.bin file)" << std::endl;
  std::cout << "  -o        Output project folder" << std::endl;
  std::cout << "  -g 1|0    Use geometric consistency (0: disable, 1: enable)         [ default: 0   ]" << std::endl;
  std::cout << "--------------------------------------------------------------------------------------" << std::endl;
  std::cout << " MVS parameters  " << std::endl;
  std::cout << std::endl;
  std::cout << "  -mc       Number of montecarlo samples                              [ default: 15  ]" << std::endl;
  std::cout << "  -n        Number of iterations                                      [ default: 5   ]" << std::endl;
  std::cout << "  -w        Window size                                               [ default: 5   ]" << std::endl;
  std::cout << "  -pe       Pertubation weight                                        [ default: 1   ]" << std::endl;
  std::cout << "  -pc 1|0   Photometric consistency filtering (0: disable, 1: enable) [ default: 1   ]" << std::endl;
  std::cout << "  -qz       Minimum filtering photometric probability                 [ default: 0.2 ]" << std::endl;
  std::cout << "  -s        Minimum view threshold for filtering                      [ default: 2   ]" << std::endl;
  std::cout << "  -fg 1|0   Geometric consistency filtering (0: disable, 1: enable)   [ default: 1   ]" << std::endl;
  std::cout << "--------------------------------------------------------------------------------------" << std::endl;
  std::cout << " Fusion parameters" << std::endl;
  std::cout << std::endl;
  std::cout << "  -f        Minimum points                                            [ default: 4   ]" << std::endl;
  std::cout << "  -r        Reprojection error threshold                              [ default: 1.3 ]" << std::endl;
  std::cout << std::endl;

  exit( EXIT_FAILURE );
}

/**
 * @brief Get command line arguments for global parameters 
 * 
 * @param input_path      Input SfM file path 
 * @param output_path     Output pwmvs folder path 
 * @param geometric       Use geometric 
 * @param argc            main arguments count 
 * @param argv            main arguments values 
 */
void handle_global_parameters( std::string& input_path, std::string& output_path, bool& geometric, int argc, char** argv )
{
  geometric = false;
  for ( int id_arg = 1; id_arg < argc; id_arg += 2 )
  {
    const std::string cur_arg( argv[ id_arg ] );
    if ( ( id_arg + 1 ) == argc )
    {
      std::cerr << "Missing parameter value for argument : " << cur_arg << std::endl;
      exit( EXIT_FAILURE );
    }

    if ( cur_arg == "-i" )
    {
      input_path = argv[ id_arg + 1 ];
    }
    else if ( cur_arg == "-o" )
    {
      output_path = argv[ id_arg + 1 ];
    }
    else if ( cur_arg == "-g" )
    {
      const bool g = std::stoi( argv[ id_arg + 1 ] );
      geometric    = g;
    }
  }
  if ( input_path.empty() )
  {
    std::cerr << "No input path given" << std::endl;
    usage( argv[ 0 ] );
  }
  if ( output_path.empty() )
  {
    std::cerr << "No output path given" << std::endl;
    usage( argv[ 0 ] );
  }
}

/**
 * @brief Get Command line arguments for PWMVS 
 * 
 * @param pwmvs_controller 
 * @param argc 
 * @param argv 
 */
void handle_pwmvs_arguments( PWMVS::Options& pwmvs_options, int argc, char** argv )
{
  // Default values
  pwmvs_options.monte_carlo_samples            = 15;
  pwmvs_options.num_iterations                 = 5;
  pwmvs_options.window_size                    = 5;
  pwmvs_options.perturbation                   = 1;
  pwmvs_options.filter_photometric_consistency = true;
  pwmvs_options.filter_min_color_similarity    = 0.2;
  pwmvs_options.filter_min_num_consistent      = 2;
  pwmvs_options.filter_geometric_consistency   = true;

  for ( int id_arg = 1; id_arg < argc; id_arg += 2 )
  {
    const std::string cur_arg( argv[ id_arg ] );
    if ( id_arg + 1 == argc )
    {
      std::cerr << "Missing parameter value for argument : " << cur_arg << std::endl;
      exit( EXIT_FAILURE );
    }
    else if ( cur_arg == "-mc" )
    {
      const int n                       = std::stoi( argv[ id_arg + 1 ] );
      pwmvs_options.monte_carlo_samples = n;
    }
    else if ( cur_arg == "-n" )
    {
      const int n                  = std::stoi( argv[ id_arg + 1 ] );
      pwmvs_options.num_iterations = n;
    }
    else if ( cur_arg == "-w" )
    {
      const int n               = std::stoi( argv[ id_arg + 1 ] );
      pwmvs_options.window_size = n;
    }
    else if ( cur_arg == "-pe" )
    {
      const double d             = std::stod( argv[ id_arg + 1 ] );
      pwmvs_options.perturbation = d;
    }
    else if ( cur_arg == "-pc" )
    {
      const bool g                                 = std::stoi( argv[ id_arg + 1 ] );
      pwmvs_options.filter_photometric_consistency = g;
    }
    else if ( cur_arg == "-qz" )
    {
      const double d                            = std::stod( argv[ id_arg + 1 ] );
      pwmvs_options.filter_min_color_similarity = d;
    }
    else if ( cur_arg == "-s" )
    {
      const int s                             = std::stoi( argv[ id_arg + 1 ] );
      pwmvs_options.filter_min_num_consistent = s;
    }
    else if ( cur_arg == "-fg" )
    {
      const int g                                = std::stoi( argv[ id_arg + 1 ] );
      pwmvs_options.filter_geometric_consistency = g;
    }
  }
}

/**
 * @brief Get command line arguments for fusion 
 * 
 * @param fusion_options  Options to set 
 * @param argc            main arguments count 
 * @param argv            main arguments values 
 */
void handle_fusion_arguments( Fusion::Options& fusion_options, int argc, char** argv )
{
  // Default arguments
  fusion_options.min_points             = 4;
  fusion_options.max_reprojection_error = 1.3;

  for ( int id_arg = 1; id_arg < argc; id_arg += 2 )
  {
    const std::string cur_arg( argv[ id_arg ] );
    if ( id_arg + 1 == argc )
    {
      std::cerr << "Missing parameter value for argument : " << cur_arg << std::endl;
      exit( EXIT_FAILURE );
    }

    if ( cur_arg == "-f" )
    {
      const int n               = std::stoi( argv[ id_arg + 1 ] );
      fusion_options.min_points = n;
    }
    else if ( cur_arg == "-r" )
    {
      const double d                        = std::stod( argv[ id_arg + 1 ] );
      fusion_options.max_reprojection_error = d;
    }
  }
}